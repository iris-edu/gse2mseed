/***************************************************************************
 * gse2mseed.c
 *
 * Simple waveform data conversion from GSE to Mini-SEED.
 *
 * Written by Chad Trabant, IRIS Data Management Center
 *
 * modified 2005.272
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <errno.h>

#include <libmseed.h>

#include "cm6.h"

#define VERSION "1.3"
#define PACKAGE "gse2mseed"

static int gse2group (char *gsefile, TraceGroup *mstg);
static int parameter_proc (int argcount, char **argvec);
static char *getoptval (int argcount, char **argvec, int argopt);
static void addfile (char *filename);
static void record_handler (char *record, int reclen);
static void usage (void);

static int   verbose     = 0;
static int   packreclen  = -1;
static int   encoding    = -1;
static int   byteorder   = -1;
static int   ignorecs    = 0;
static int   auxtoloc    = 0;
static char *forcenet    = 0;
static char *forceloc    = 0;
static char *outputfile  = 0;
static FILE *ofp         = 0;

struct filelink {
  char *filename;
  struct filelink *next;
};

struct filelink *filelist = 0;

int
main (int argc, char **argv)
{
  struct filelink *flp;
  TraceGroup *mstg = 0;

  int packedsamples = 0;
  int packedrecords = 0;
  
  /* Process given parameters (command line and parameter file) */
  if (parameter_proc (argc, argv) < 0)
    return -1;
  
  /* Init TraceGroup */
  mstg = mst_initgroup (mstg);
  
  /* Open the output file if specified otherwise stdout */
  if ( outputfile )
    {
      if ( strcmp (outputfile, "-") == 0 )
        {
          ofp = stdout;
        }
      else if ( (ofp = fopen (outputfile, "w")) == NULL )
        {
          fprintf (stderr, "Cannot open output file: %s (%s)\n",
                   outputfile, strerror(errno));
          return -1;
        }
    }
  else
    {
      ofp = stdout;
    }
  
  /* Read input GSE files into TraceGroup */
  flp = filelist;
  
  while ( flp != 0 )
    {
      if ( verbose )
	fprintf (stderr, "Reading %s\n", flp->filename);

      gse2group (flp->filename, mstg);
      
      flp = flp->next;
    }

  /* Pack data into Mini-SEED records */
  packedrecords = mst_packgroup (mstg, &record_handler, packreclen, encoding, byteorder,
				 &packedsamples, 1, verbose-1, NULL);
  
  if ( packedrecords < 0 )
    {
      fprintf (stderr, "Error packing data\n");
    }
  else
    {
      fprintf (stderr, "Packed %d trace(s) of %d samples into %d records\n",
	       mstg->numtraces, packedsamples, packedrecords);
    }
  
  /* Make sure everything is cleaned up */
  mst_freegroup (&mstg);
  
  if ( ofp )
    fclose (ofp);
  
  return 0;
}  /* End of main() */


/***************************************************************************
 * gse2group:
 * Read a GSE file and add data samples to a TraceGroup.  As the GSE
 * is read in a MSrecord struct is used as a holder for the input
 * information.
 *
 * Returns 0 on success, and -1 on failure
 ***************************************************************************/
static int
gse2group (char *gsefile, TraceGroup *mstg)
{
  FILE *ifp;
  MSrecord *msr = 0;
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
  int cm6bufsize = 0;
  int32_t *intbuf = 0;
  int32_t intbufsize = 0;
  
  /* Open input file */
  if ( (ifp = fopen (gsefile, "r")) == NULL )
    {
      fprintf (stderr, "Cannot open input file: %s (%s)\n",
	       gsefile, strerror(errno));
      return -1;
    }
  
  if ( ! (msr = msr_init(msr)) )
    {
      fprintf (stderr, "Cannot initialize MSrecord strcture\n");
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
	  ms_strncpclean (msr->channel, line + 35, 3);
	  
	  if ( forceloc )
	    {
	      ms_strncpclean (msr->location, forceloc, 2);
	    }	      
	  else if ( auxtoloc )
	    {
	      ms_strncpclean (msr->location, line + 39, 2);
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
      
      else if ( ! strncmp ("CHK2", line, 4) )
	{
	  if ( linesize < 6 )
	    {
	      fprintf (stderr, "[%s] %s %s: CHK2 line is too short, only %d characters:\n%s\n",
		       gsefile, msr->station, msr->channel, linesize, line);
	      retval = -1;
	      break;
	    }
	  
	  /* Parse original chksum from the line */
	  ms_strncpclean (chkstr, line + 5, 8);
	  ochksum = strtol (chkstr, NULL, 10);
	  
	  /* Unpack CM6 */
	  if ( format == 1 )
	    {
	      if ( (intbuf = unpackcm6 (cm6buf, cm6bufsize, intbuf, &intbufsize, -1, 2)) == NULL )
		{
		  fprintf (stderr, "[%s] %s %s: Error unpacking CM6 compressed data\n",
			   gsefile, msr->station, msr->channel);
		  retval = -1;
		  break;
		}
	    }
	  
	  if ( msr->samplecnt != intbufsize )
	    {
	      fprintf (stderr, "[%s] %s %s: Unpacked %d of %d samples!\n",
		       gsefile, msr->station, msr->channel, intbufsize, msr->samplecnt);
	      msr->samplecnt = intbufsize;
	    }
	  
	  if ( verbose >= 3 )
	    {
	      int tint;

	      printf ("[%s] %s %s: First 6 samples:\n",
		      gsefile, msr->station, msr->channel);
	      
	      for ( tint = 0; tint < 6; tint++ )
		{
		  printf ("%10d ", *(intbuf+tint));
		}
	      printf ("\n");
	    }
	  
	  /* Compute chksum and compare */
	  cchksum = gsechksum (intbuf, intbufsize);
	  
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
	  char *tptr;
	  
	  datalinesize = linesize;
	  
	  /* Truncate at first newline character */
	  if ( (tptr = strchr (line, '\n')) )
	    {
	      *tptr = '\0';
	      datalinesize = tptr - line;
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
	  
	  /* Process INT data, one ASCII integer per line is expected */
	  if ( format == 2 )
	    {
	      tptr = &line[0];
	      
	      /* Remove leading space character(s) */
	      while ( *tptr == ' ' )
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
	      if ( intbufsize == 0 )
		{
		  intbuf = realloc (intbuf, sizeof(int32_t) * msr->samplecnt);
		}
	      
	      if ( (intbufsize+1) > msr->samplecnt )
		{
		  fprintf (stderr, "[%s] %s %s: More than %d INT samples found in input file\n",
			   gsefile, msr->station, msr->channel, msr->samplecnt);
		  retval = -1;
		  break;
		}
	      
	      *(intbuf + intbufsize) = (int32_t) strtol (tptr, NULL, 10);
	      
	      intbufsize++;
	    }
	}
      if ( blockend )
	{
	  /* Add data to TraceGroup */
	  msr->datasamples = intbuf;
	  msr->numsamples = intbufsize;
	  msr->sampletype = 'i';
	  
	  if ( verbose >= 1 )
	    {
	      fprintf (stderr, "[%s] %d samps @ %.6f Hz for N: '%s', S: '%s', L: '%s', C: '%s'\n",
		       gsefile, msr->numsamples, msr->samprate,
		       msr->network, msr->station,  msr->location, msr->channel);
	    }
	  
	  if ( ! mst_addmsrtogroup (mstg, msr, -1.0, -1.0) )
	    {
	      fprintf (stderr, "[%s] Error adding samples to TraceGroup\n", gsefile);
	    }
	  
	  /* Cleanup and reset state */
	  msr->datasamples = 0;
	  msr = msr_init (msr);
	  cm6bufsize = 0;
	  intbufsize = 0;
	  ochksum = 0;
	  blockend = 0;
	}
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
      else if (strcmp (argvec[optind], "-n") == 0)
	{
	  forcenet = getoptval(argcount, argvec, optind++);
	}
      else if (strcmp (argvec[optind], "-l") == 0)
	{
	  forceloc = getoptval(argcount, argvec, optind++);
	}
      else if (strcmp (argvec[optind], "-a") == 0)
	{
	  auxtoloc = 1;
	}
      else if (strcmp (argvec[optind], "-i") == 0)
	{
	  ignorecs = 1;
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

  /* Make sure an input files were specified */
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
  }
  
  /* Special case of '-o -' usage */
  if ( (argopt+1) < argcount && strcmp (argvec[argopt], "-o") == 0 )
    if ( strcmp (argvec[argopt+1], "-") == 0 )
      return argvec[argopt+1];
  
  if ( (argopt+1) < argcount && *argvec[argopt+1] != '-' )
    return argvec[argopt+1];
  
  fprintf (stderr, "Option %s requires a value\n", argvec[argopt]);
  exit (1);
}  /* End of getoptval() */


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
record_handler (char *record, int reclen)
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
  fprintf (stderr, "Convert GSE2.x/IMS1.0 CM6 and INT waveform data to Mini-SEED.\n");
  fprintf (stderr, "GSE structures recognized: WID2, STA2, DAT2 and CHK2, all other\n");
  fprintf (stderr, "information in the input file(s) is ignored.\n\n");
  fprintf (stderr, "Usage: %s [options] file1 file2 file3...\n\n", PACKAGE);
  fprintf (stderr,
	   " ## Options ##\n"
	   " -V             Report program version\n"
	   " -h             Show this usage message\n"
	   " -v             Be more verbose, multiple flags can be used\n"
	   " -n netcode     Force the SEED network code\n"
	   " -l loccode     Force the SEED location code\n"
	   " -a             Map GSE auxiliary id code to SEED location id code\n"
	   " -i             Ignore GSE checksum mismatch, warn but continue\n"
	   " -r bytes       Specify record length in bytes for packing, default: 4096\n"
	   " -e encoding    Specify SEED encoding format for packing, default: 11 (Steim2)\n"
	   " -b byteorder   Specify byte order for packing, MSBF: 1 (default), LSBF: 0\n"
	   " -o outfile     Specify the output file, default is stdout.\n"
	   "\n"
	   " file(s)        File(s) of GSE input data\n"
	   "\n"
	   "Supported Mini-SEED encoding formats:\n"
	   " 1  : 16-bit integers (only works if samples can be represented in 16-bits)\n"
	   " 3  : 32-bit integers\n"
	   " 10 : Steim 1 compression\n"
	   " 11 : Steim 2 compression\n"
	   "\n");
}  /* End of usage() */
