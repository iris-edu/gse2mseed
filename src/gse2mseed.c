/***************************************************************************
 * gse2mseed.c
 *
 * Simple waveform data conversion from GSE to Mini-SEED.
 *
 * Written by Chad Trabant, IRIS Data Management Center
 *
 * modified 2013.138
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <errno.h>

#include <libmseed.h>

#include "cm6.h"

#define VERSION "1.12"
#define PACKAGE "gse2mseed"

static void packtraces (flag flush);
static int gse2group (char *gsefile, MSTraceGroup *mstg);
static int parameter_proc (int argcount, char **argvec);
static char *getoptval (int argcount, char **argvec, int argopt);
static int readlistfile (char *listfile);
static void addfile (char *filename);
static void record_handler (char *record, int reclen, void *handlerdata);
static void usage (void);

static int   verbose     = 0;
static int   ignorecs    = 0;
static char  bufferall   = 0;
static char *forcenet    = 0;
static char *forceloc    = 0;
static int   packreclen  = -1;
static int   encoding    = -1;
static int   byteorder   = -1;
static char *outputfile  = 0;
static FILE *ofp         = 0;

struct filelink {
  char *filename;
  struct filelink *next;
};

/* A list of input files */
struct filelink *filelist = 0;

static MSTraceGroup *mstg = 0;

static int packedtraces  = 0;
static int packedsamples = 0;
static int packedrecords = 0;


int
main (int argc, char **argv)
{
  struct filelink *flp;
  
  /* Process given parameters (command line and parameter file) */
  if (parameter_proc (argc, argv) < 0)
    return -1;
  
  /* Init MSTraceGroup */
  mstg = mst_initgroup (mstg);
  
  /* Open the output file if specified otherwise stdout */
  if ( outputfile )
    {
      if ( strcmp (outputfile, "-") == 0 )
        {
          ofp = stdout;
        }
      else if ( (ofp = fopen (outputfile, "wb")) == NULL )
        {
          fprintf (stderr, "Cannot open output file: %s (%s)\n",
                   outputfile, strerror(errno));
          return -1;
        }
    }
  
  /* Read input GSE files into MSTraceGroup */
  flp = filelist;
  
  while ( flp != 0 )
    {
      if ( verbose )
	fprintf (stderr, "Reading %s\n", flp->filename);

      gse2group (flp->filename, mstg);
      
      flp = flp->next;
    }
  
  /* Pack any remaining, possibly all data */
  packtraces (1);
  packedtraces += mstg->numtraces;
  
  fprintf (stderr, "Packed %d trace(s) of %d samples into %d records\n",
           packedtraces, packedsamples, packedrecords);
  
  /* Make sure everything is cleaned up */
  mst_freegroup (&mstg);
  
  if ( ofp )
    fclose (ofp);
  
  return 0;
}  /* End of main() */


/***************************************************************************
 * packtraces:
 *
 * Pack all traces in a group.
 *
 * Returns 0 on success, and -1 on failure
 ***************************************************************************/
static void
packtraces (flag flush)
{
  MSTrace *mst;
  MSRecord *msr = NULL;
  int64_t trpackedsamples = 0;
  int64_t trpackedrecords = 0;
  struct blkt_1000_s Blkt1000;
  struct blkt_1001_s Blkt1001;
  
  mst = mstg->traces;
  while ( mst )
    {
      if ( mst->numsamples <= 0 )
        {
          mst = mst->next;
          continue;
        }
      
      /* Initialize MSRecord template for packing */
      msr = msr_init(msr);
      strncpy (msr->network, mst->network, sizeof(msr->network));
      strncpy (msr->station, mst->station, sizeof(msr->station));
      strncpy (msr->location, mst->location, sizeof(msr->location));
      strncpy (msr->channel, mst->channel, sizeof(msr->channel));
      
      /* Add blockettes 1000 & 1001 to template */
      memset (&Blkt1000, 0, sizeof(struct blkt_1000_s));
      msr_addblockette (msr, (char *) &Blkt1000, sizeof(struct blkt_1001_s), 1000, 0);
      memset (&Blkt1001, 0, sizeof(struct blkt_1001_s));
      msr_addblockette (msr, (char *) &Blkt1001, sizeof(struct blkt_1001_s), 1001, 0);
      
      trpackedrecords = mst_pack (mst, &record_handler, 0, packreclen, encoding, byteorder,
                                  &trpackedsamples, flush, verbose-2, msr);
      
      if ( trpackedrecords < 0 )
        {
          fprintf (stderr, "Error packing data\n");
        }
      else
        {
          packedrecords += trpackedrecords;
          packedsamples += trpackedsamples;
        }
      
      mst = mst->next;
    }
  
  msr_free (&msr);
}  /* End of packtraces() */


/***************************************************************************
 * gse2group:
 * Read a GSE file and add data samples to a MSTraceGroup.  As the GSE
 * is read in a MSRecord struct is used as a holder for the input
 * information.
 *
 * Returns 0 on success, and -1 on failure
 ***************************************************************************/
static int
gse2group (char *gsefile, MSTraceGroup *mstg)
{
  FILE *ifp;
  MSRecord *msr = 0;
  char line[1025];
  int linesize;

  char timestr[25];
  char sampstr[10];
  char ratestr[15];
  char chkstr[10];
  
  int expectdata = 0;
  int ochksum, cchksum;
  int blockend = 0;
  int retval = 0;
  int format = 0;  /* 1 = CM6, 2 = INT */
  
  char *cm6buf = 0;
  int cm6bufsize = 0;   /* Number of characters in cm6buf */
  int32_t *intbuf = 0;
  int intbufcount = 0;  /* Number of 32-bit integers in intbuf */
  
  /* Open input file */
  if ( (ifp = fopen (gsefile, "r")) == NULL )
    {
      fprintf (stderr, "Cannot open input file: %s (%s)\n",
	       gsefile, strerror(errno));
      return -1;
    }
  
  /* Open .mseed output file if needed, replacing .gse if present */
  if ( ! ofp )
    {
      char mseedoutputfile[1024];
      int filelen;
      strncpy (mseedoutputfile, gsefile, sizeof(mseedoutputfile)-6 );
      filelen = strlen (mseedoutputfile);
      
      /* Truncate file name if .gse is at the end */
      if ( filelen > 4 )
	if ( (*(mseedoutputfile + filelen - 1) == 'e' || *(mseedoutputfile + filelen - 1) == 'E') &&
	     (*(mseedoutputfile + filelen - 2) == 's' || *(mseedoutputfile + filelen - 2) == 'S') &&
	     (*(mseedoutputfile + filelen - 3) == 'g' || *(mseedoutputfile + filelen - 3) == 'G') &&
	     (*(mseedoutputfile + filelen - 4) == '.') )
	  {
	    *(mseedoutputfile + filelen - 4) = '\0';
	  }

      strcat (mseedoutputfile, ".mseed");
      
      if ( (ofp = fopen (mseedoutputfile, "wb")) == NULL )
        {
          fprintf (stderr, "Cannot open output file: %s (%s)\n",
                   mseedoutputfile, strerror(errno));
          return -1;
        }
    }

  if ( ! (msr = msr_init(msr)) )
    {
      fprintf (stderr, "Cannot initialize MSRecord strcture\n");
      return -1;
    }
  
  while ( fgets (line, 1025, ifp) )
    {
      linesize = strlen (line);
      
      if ( ! strncmp ("WID2", line, 4) && ! expectdata )
	{
	  if ( linesize < 68 ) 
	    {
	      fprintf (stderr, "[%s] WID2 line is too short, only %d characters:\n%s\n",
		       gsefile, linesize, line);
	      retval = -1;
	      break;
	    }
	  
	  /* Extract values from WID2 line and populate the msr holder */
	  strncpy (timestr, line + 5, 23); timestr[23] = '\0';
	  msr->starttime = ms_timestr2hptime (timestr);
	  
	  ms_strncpclean (msr->network, forcenet, 2);
	  ms_strncpclean (msr->station, line + 29, 5);
	  ms_strncpclean (msr->location, line + 39, 2);
	  ms_strncpclean (msr->channel, line + 35, 3);
	  
	  if ( forceloc )
	    {
	      ms_strncpclean (msr->location, forceloc, 2);
	    }
	  
	  if ( ! strncmp (line + 44, "CM6", 3) )
	    {
	      format = 1;
	    }
	  else if ( ! strncmp (line + 44, "INT", 3) )
	    {
	      format = 2;
	    }
	  else
	    {
	      fprintf (stderr, "[%s] %s %s: Only CM6 and INT formatted data are supported, not %.3s\n",
		       gsefile, msr->station, msr->channel, line + 44);
	      retval = -1;
	      break;
	    }
	  
	  ms_strncpclean (sampstr, line + 48, 8);
	  msr->samplecnt = strtol (sampstr, NULL, 10);
	  
	  ms_strncpclean (ratestr, line + 57, 11);
	  msr->samprate = strtod (ratestr, NULL);
	}
      
      else if ( ! strncmp ("STA2", line, 4) && ! expectdata )
	{
	  if ( linesize < 14 )
	    {
	      fprintf (stderr, "[%s] %s %s: STA2 line is too short, only %d characters:\n%s\n",
		       gsefile, msr->station, msr->channel, linesize, line);
	      retval = -1;
	      break;
	    }

	  if ( ! forcenet )
	    {
	      ms_strncpclean (msr->network, line + 5, 2);
	    }
	}
      
      else if ( ! strncmp ("DAT2", line, 4) && ! expectdata )
	{
	  if ( format == 0 )
	    {
	      fprintf (stderr, "[%s] DAT2 line read but data format is not yet known!\n",
		       gsefile);
	      retval = -1;
	      break;
	    }
	  
	  expectdata = 1;
	}
      
      else if ( ! strncmp ("CHK2 ", line, 5) )
	{
	  if ( linesize < 6 )
	    {
	      fprintf (stderr, "[%s] %s %s: CHK2 line is too short, only %d characters:\n%s\n",
		       gsefile, msr->station, msr->channel, linesize, line);
	      retval = -1;
	      break;
	    }
	  
	  /* Test that data was actually expected */
	  if ( ! expectdata )
	    {
	      fprintf (stderr, "[%s] %s %s: CHK2 was found but no DAT2 line indicated the start of data\n",
		       gsefile, msr->station, msr->channel);
	      retval = -1;
	      break;
	    }
	  
	  /* Parse original chksum from the line */
	  ms_strncpclean (chkstr, line + 5, 8);
	  ochksum = strtol (chkstr, NULL, 10);
	  
	  /* Unpack CM6 */
	  if ( format == 1 )
	    {
	      if ( (intbuf = unpackcm6 (cm6buf, cm6bufsize, intbuf, &intbufcount, -1, 2)) == NULL )
		{
		  fprintf (stderr, "[%s] %s %s: Error unpacking CM6 compressed data\n",
			   gsefile, msr->station, msr->channel);
		  retval = -1;
		  break;
		}
	    }
	  
	  if ( msr->samplecnt != intbufcount )
	    {
	      fprintf (stderr, "[%s] %s %s: Unpacked %d of %lld samples!\n",
		       gsefile, msr->station, msr->channel, intbufcount,
		       (long long int) msr->samplecnt);
	      msr->samplecnt = intbufcount;
	    }
	  
	  if ( verbose >= 3 )
	    {
	      int tint;
	      
	      fprintf (stderr, "[%s] %s %s: First 6 samples:\n",
		       gsefile, msr->station, msr->channel);
	      
	      for ( tint = 0; tint < 6 && tint < msr->samplecnt; tint++ )
		{
		  fprintf (stderr, "%10d ", *(intbuf+tint));
		}
	      fprintf (stderr, "\n");
	    }
	  
	  /* Compute chksum and compare */
	  cchksum = gsechksum (intbuf, intbufcount);
	  
	  if ( ochksum != cchksum )
	    {
	      fprintf (stderr, "[%s] %s %s: Calculated chksum does not match chksum from CHK2 line\n",
		       gsefile, msr->station, msr->channel);
	      fprintf (stderr, "Original: %d, Calculated: %d\n", ochksum, cchksum);
	      
	      if ( ! ignorecs )
		{
		  retval = -1;
		  break;
		}
	    }
	  
	  blockend = 1;
	  expectdata = 0;
	  format = 0;
	}
      
      /* Read in data lines */
      else if ( expectdata )
	{
	  int datalinesize = 0;
	  int spanned;
	  char *tptr;
	  
	  datalinesize = linesize;
	  
	  /* Truncate at first newline or carriage return character */
	  if ( (spanned = strcspn (line, "\n\r")) != datalinesize )
	    {
	      line[spanned] = '\0';
	      datalinesize = spanned;
	    }
	  
	  /* Process CM6 data */
	  if ( format == 1 )
	    {
	      /* Truncate at first space character */
	      if ( (tptr = strchr (line, ' ')) )
		{
		  *tptr = '\0';
		  datalinesize = tptr - line;
		}
	      
	      /* Skip empty lines */
	      if ( datalinesize <= 0 )
		{
		  continue;
		}
	      
	      /* Make sure all characters are in the CM6 character table */
	      if ( strspn(line, cm6table) != datalinesize )
		{
		  fprintf (stderr, "[%s] %s %s: Expected a line with CM6 characters but got:\n'%s'\n",
			   gsefile, msr->station, msr->channel, line);
		  retval = -1;
		  break;
		}
	      
	      /* Add new data to CM6 buffer */
	      cm6buf = realloc (cm6buf, cm6bufsize + datalinesize);
	      
	      memcpy ((cm6buf + cm6bufsize), line, datalinesize);
	      
	      cm6bufsize += datalinesize;
	    }
	  
	  /* Process INT data, one to many ASCII integer(s) per line are possible */
	  if ( format == 2 )
	    {
	      tptr = &line[0];
	      
	      /* Skip leading space character(s) */
	      while ( isspace (*tptr) )
		{
		  tptr++;
		  datalinesize--;
		}
	      
	      /* Skip empty lines */
	      if ( datalinesize <= 0 )
		{
		  continue;
		}
	      
	      /* (Re)Allocate intbuf if needed */
	      if ( intbufcount == 0 )
		{
		  intbuf = realloc (intbuf, sizeof(int32_t) * msr->samplecnt);
		}
	      
	      if ( (intbufcount+1) > msr->samplecnt )
		{
		  fprintf (stderr, "[%s] %s %s: More than %lld INT samples found in input file\n",
			   gsefile, msr->station, msr->channel,
			   (long long int) msr->samplecnt);
		  retval = -1;
		  break;
		}
	      
	      while ( *tptr )
		{
		  *(intbuf + intbufcount) = (int32_t) strtol (tptr, NULL, 10);
		  
		  intbufcount++;
		  
		  /* Skip to next space character and then to next non-space character */
		  while ( *tptr && ! isspace (*tptr) )
		    tptr++;
		  while ( *tptr && isspace (*tptr) )
		    tptr++;
		}
	    }
	}
      
      if ( blockend )
	{
	  /* Add data to MSTraceGroup */
	  msr->datasamples = intbuf;
	  msr->numsamples = intbufcount;
	  msr->sampletype = 'i';
	  
	  if ( verbose >= 1 )
	    {
	      fprintf (stderr, "[%s] %lld samps @ %.6f Hz for N: '%s', S: '%s', L: '%s', C: '%s'\n",
		       gsefile, (long long int) msr->numsamples, msr->samprate,
		       msr->network, msr->station, msr->location, msr->channel);
	    }
	  
	  if ( ! mst_addmsrtogroup (mstg, msr, 0, -1.0, -1.0) )
	    {
	      fprintf (stderr, "[%s] Error adding samples to MSTraceGroup\n", gsefile);
	    }
	  
	  /* Unless buffering all files in memory pack any MSTraces now */
          if ( ! bufferall )
            {
              packtraces (1);
              packedtraces += mstg->numtraces;
              mst_initgroup (mstg);
            }

	  /* Cleanup and reset state */
	  msr->datasamples = 0;
	  msr = msr_init (msr);

	  cm6bufsize = 0;
	  intbufcount = 0;
	  ochksum = 0;
	  blockend = 0;
	}
    }

  fclose (ifp);
  
  if ( ofp  && ! outputfile )
    {
      fclose (ofp);
      ofp = 0;
    }

  if ( cm6buf )
    free (cm6buf);
  
  if ( intbuf )
    free (intbuf);
  
  if ( msr )
    msr_free (&msr);
      
  return retval;
}  /* End of gse2group() */


/***************************************************************************
 * parameter_proc:
 * Process the command line parameters.
 *
 * Returns 0 on success, and -1 on failure
 ***************************************************************************/
static int
parameter_proc (int argcount, char **argvec)
{
  int optind;

  /* Process all command line arguments */
  for (optind = 1; optind < argcount; optind++)
    {
      if (strcmp (argvec[optind], "-V") == 0)
	{
	  fprintf (stderr, "%s version: %s\n", PACKAGE, VERSION);
	  exit (0);
	}
      else if (strcmp (argvec[optind], "-h") == 0)
	{
	  usage();
	  exit (0);
	}
      else if (strncmp (argvec[optind], "-v", 2) == 0)
	{
	  verbose += strspn (&argvec[optind][1], "v");
	}
      else if (strcmp (argvec[optind], "-i") == 0)
	{
	  ignorecs = 1;
	}
      else if (strcmp (argvec[optind], "-B") == 0)
	{
	  bufferall = 1;
	}
      else if (strcmp (argvec[optind], "-n") == 0)
	{
	  forcenet = getoptval(argcount, argvec, optind++);
	}
      else if (strcmp (argvec[optind], "-l") == 0)
	{
	  forceloc = getoptval(argcount, argvec, optind++);
	}
      else if (strcmp (argvec[optind], "-r") == 0)
	{
	  packreclen = atoi (getoptval(argcount, argvec, optind++));
	}
      else if (strcmp (argvec[optind], "-e") == 0)
	{
	  encoding = atoi (getoptval(argcount, argvec, optind++));
	}
      else if (strcmp (argvec[optind], "-b") == 0)
	{
	  byteorder = atoi (getoptval(argcount, argvec, optind++));
	}
      else if (strcmp (argvec[optind], "-o") == 0)
	{
	  outputfile = getoptval(argcount, argvec, optind++);
	}
      else if (strncmp (argvec[optind], "-", 1) == 0 &&
	       strlen (argvec[optind]) > 1 )
	{
	  fprintf(stderr, "Unknown option: %s\n", argvec[optind]);
	  exit (1);
	}
      else
	{
	  addfile (argvec[optind]);
	}
    }
  
  /* Make sure an output file is specified if buffering all */
  if ( bufferall && ! outputfile )
    {
      fprintf (stderr, "Need to specify output file with -o if using -B\n");
      exit(1);
    }
  
  /* Make sure input files were specified */
  if ( filelist == 0 )
    {
      fprintf (stderr, "No input files were specified\n\n");
      fprintf (stderr, "%s version %s\n\n", PACKAGE, VERSION);
      fprintf (stderr, "Try %s -h for usage\n", PACKAGE);
      exit (1);
    }
  
  /* Report the program version */
  if ( verbose )
    fprintf (stderr, "%s version: %s\n", PACKAGE, VERSION);
  
  /* Check the input files for any list files, if any are found
   * remove them from the list and add the contained list */
  if ( filelist )
    {
      struct filelink *prevlp, *lp;
      
      prevlp = lp = filelist;
      while ( lp != 0 )
	{
	  if ( *(lp->filename) == '@' )
	    {
	      /* Remove this node from the list */
	      if ( lp == filelist )
		filelist = lp->next;
	      else
		prevlp->next = lp->next;
	      
	      /* Read list file, skip the '@' first character */
	      readlistfile (lp->filename + 1);
	      
	      /* Free memory for this node */
	      free (lp->filename);
	      free (lp);
	    }
	  else
	    {
	      prevlp = lp;
	    }
	  
	  lp = lp->next;
	}
    }
  
  return 0;
}  /* End of parameter_proc() */


/***************************************************************************
 * getoptval:
 * Return the value to a command line option; checking that the value is 
 * itself not an option (starting with '-') and is not past the end of
 * the argument list.
 *
 * argcount: total arguments in argvec
 * argvec: argument list
 * argopt: index of option to process, value is expected to be at argopt+1
 *
 * Returns value on success and exits with error message on failure
 ***************************************************************************/
static char *
getoptval (int argcount, char **argvec, int argopt)
{
  if ( argvec == NULL || argvec[argopt] == NULL ) {
    fprintf (stderr, "getoptval(): NULL option requested\n");
    exit (1);
    return 0;
  }
  
  /* Special case of '-o -' usage */
  if ( (argopt+1) < argcount && strcmp (argvec[argopt], "-o") == 0 )
    if ( strcmp (argvec[argopt+1], "-") == 0 )
      return argvec[argopt+1];
  
  if ( (argopt+1) < argcount && *argvec[argopt+1] != '-' )
    return argvec[argopt+1];
  
  fprintf (stderr, "Option %s requires a value\n", argvec[argopt]);
  exit (1);
  return 0;
}  /* End of getoptval() */


/***************************************************************************
 * readlistfile:
 * Read a list of files from a file and add them to the filelist for
 * input data.
 *
 * Returns the number of file names parsed from the list or -1 on error.
 ***************************************************************************/
static int
readlistfile (char *listfile)
{
  FILE *fp;
  char line[1024];
  char *ptr;
  int  filecnt = 0;
  int  nonspace;

  char filename[1024];
  int  fields;
  
  /* Open the list file */
  if ( (fp = fopen (listfile, "rb")) == NULL )
    {
      if (errno == ENOENT)
        {
          fprintf (stderr, "Could not find list file %s\n", listfile);
          return -1;
        }
      else
        {
          fprintf (stderr, "Error opening list file %s: %s\n",
                   listfile, strerror (errno));
          return -1;
        }
    }
  
  if ( verbose )
    fprintf (stderr, "Reading list of input files from %s\n", listfile);
  
  while ( (fgets (line, sizeof(line), fp)) !=  NULL)
    {
      /* Truncate line at first \r or \n and count non-space characters */
      nonspace = 0;
      ptr = line;
      while ( *ptr )
        {
          if ( *ptr == '\r' || *ptr == '\n' || *ptr == '\0' )
            {
              *ptr = '\0';
              break;
            }
          else if ( *ptr != ' ' )
            {
              nonspace++;
            }
          
          ptr++;
        }
      
      /* Skip empty lines */
      if ( nonspace == 0 )
        continue;
      
      fields = sscanf (line, "%s", filename);
      
      if ( fields != 1 )
	{
	  fprintf (stderr, "Error parsing filename from: %s\n", line);
	  continue;
	}
      
      if ( verbose > 1 )
	fprintf (stderr, "Adding '%s' to input file list\n", filename);
      
      addfile (filename);
      filecnt++;
      
      continue;
    }
  
  fclose (fp);
  
  return filecnt;
}  /* End readlistfile() */


/***************************************************************************
 * addfile:
 *
 * Add file to end of the global file list (filelist).
 ***************************************************************************/
static void
addfile (char *filename)
{
  struct filelink *lastlp, *newlp;
  
  if ( filename == NULL )
    {
      fprintf (stderr, "addfile(): No file name specified\n");
      return;
    }
  
  lastlp = filelist;
  while ( lastlp != 0 )
    {
      if ( lastlp->next == 0 )
        break;
      
      lastlp = lastlp->next;
    }
  
  newlp = (struct filelink *) malloc (sizeof (struct filelink));
  newlp->filename = strdup(filename);
  newlp->next = 0;
  
  if ( lastlp == 0 )
    filelist = newlp;
  else
    lastlp->next = newlp;
  
}  /* End of addfile() */


/***************************************************************************
 * record_handler:
 * Saves passed records to the output file.
 ***************************************************************************/
static void
record_handler (char *record, int reclen, void *handlerdata)
{
  if ( fwrite(record, reclen, 1, ofp) != 1 )
    {
      fprintf (stderr, "Error writing to output file\n");
    }
}  /* End of record_handler() */


/***************************************************************************
 * usage:
 * Print the usage message and exit.
 ***************************************************************************/
static void
usage (void)
{
  fprintf (stderr, "%s version: %s\n\n", PACKAGE, VERSION);
  fprintf (stderr, "Convert GSE2.x/IMS1.0 INT and CM6 waveform data to Mini-SEED.\n");
  fprintf (stderr, "GSE structures recognized: WID2, STA2, DAT2 and CHK2, all other\n");
  fprintf (stderr, "information in the input file(s) is ignored.\n\n");
  fprintf (stderr, "Usage: %s [options] file1 [file2 file3 ...]\n\n", PACKAGE);
  fprintf (stderr,
	   " ## Options ##\n"
	   " -V             Report program version\n"
	   " -h             Show this usage message\n"
	   " -v             Be more verbose, multiple flags can be used\n"
	   " -i             Ignore GSE checksum mismatch, warn but continue\n"
	   " -B             Buffer data before packing, default packs at end of each block\n"
	   " -n netcode     Specify the SEED network code\n"
	   " -l locid       Specify the SEED location ID\n"
	   " -r bytes       Specify record length in bytes for packing, default: 4096\n"
	   " -e encoding    Specify SEED encoding format for packing, default: 11 (Steim2)\n"
	   " -b byteorder   Specify byte order for packing, MSBF: 1 (default), LSBF: 0\n"
	   " -o outfile     Specify the output file, default is <inputfile>.mseed\n"
	   "\n"
	   " file(s)        File(s) of GSE input data\n"
           "                  If a file is prefixed with an '@' it is assumed to contain\n"
           "                  a list of data files to be read, one file per line.\n"
	   "\n"
	   "Supported Mini-SEED encoding formats:\n"
	   " 1  : 16-bit integers (only works if samples can be represented in 16-bits)\n"
	   " 3  : 32-bit integers\n"
	   " 10 : Steim 1 compression\n"
	   " 11 : Steim 2 compression\n"
	   "\n");
}  /* End of usage() */
