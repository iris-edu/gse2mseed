
/*
 * Modified from the original by Chad Trabant, IRIS Data Management Center.
 *
 * 2005.123:
 *   Update packcm6 to use the cm6table defined in cm6.h
 *   Update unpackcm6 to use static variables and only calculate the CM6
 *   decompression table once per session; a minor optimization.
 *
 * 2004.198:
 *   Change unpackcm6 to take Nreq, the number of integers to unpack and
 *   modify it so the number of integers does not need to be known a priori.
 *
 * 2004.182:
 *   Add gsechksum() (chk2() from libsismoutil).
 *
 * 2004.180:
 *   Changed variable names and code style (added lots of spaces) for clarity.
 *   Minor changes for safety and add use comments.
 *
 * Original copyright notice follows.
 *
 * This file originally came from libsismoutil.
 *
 * Copyright (C) 2002, 2003, Sebastien Judenherc <sebastien.judenherc@na.infn.it>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
 * USA
 *
 *
 * Modifed: 2005.123
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>

#include "cm6.h"


/***************************************************************************
 * delta:
 *
 * Calculate single order differences for an array of 32-bit integers
 * in place.
 *
 ***************************************************************************/
static void
delta (int32_t *intbuf, int32_t Nint)
{
  register int32_t i;
  int32_t m1=intbuf[0];
  int32_t tmpi=intbuf[1];

  for (i=1; i<Nint; i++)
    {
      tmpi=intbuf[i];
      intbuf[i]=tmpi-m1;
      m1=tmpi;
    }
}


/***************************************************************************
 * undelta:
 *
 * Undo single order differencing for an array of 32-bit integers in
 * place.
 *
 ***************************************************************************/
static void
undelta (int32_t *intbuf, int32_t Nint)
{
  register int32_t i;
  
  for (i=1; i<Nint; i++)
    {
      intbuf[i]+=intbuf[i-1];
    }
}


/***************************************************************************
 * packcm6:
 *
 * Pack 32-bit integers into CM6 compressed data.
 *
 * intbuf  = input 32-bit integers
 * Nint    = number of input integers to pack
 * cm6buf  = CM6 ASCII buffer (will be appended to, final is returned)
 * Ncm6    = number of characters in cm6buf (will be updated)
 * ndiff   = differencing level (2 is highly recommended)
 *
 * Data in the intbuf array will be modified and should, after this
 * function returns, be considered junk by the calling process.
 *
 * cm6buf should be either allocated by the malloc(3) family of
 * functions or set to NULL by the calling program.  The data in
 * cm6buf will be appeneded to if Ncm6 correctly indicates the size.
 *
 * Returns a pointer to the resulting integers on success and NULL on  error.
 ***************************************************************************/
char *
packcm6 (int32_t *intbuf, int32_t Nint,
	 char *cm6buf, int32_t *Ncm6, int32_t ndiff)
{
  char subcm6[6];
  int32_t iint, nbt;
  int32_t i, j, k;
  
  for (j=0; j < ndiff; j++) delta(intbuf, Nint);
  
  for (i=0; i < Nint; i++)
    {
      iint = (intbuf[i] > 0) ? intbuf[i] : -intbuf[i];
      
      for (nbt=4; ((1<<nbt)<=iint) && (nbt<29); nbt+=5) ;
      
      k = 0;
      nbt -= 4;                  /* first value on 4 bits */

      do
	{
	  subcm6[k++] = ((iint>>nbt) & 31) + 32;
	  nbt -= 5;              /* next on 5 bits */
	} while (nbt >= 0);
      
      subcm6[k-1] -= 32;         /* last does not continue */
      subcm6[0] += (intbuf[i]<0) ? 16 : 0;
      
      cm6buf = (char *) realloc (cm6buf, sizeof(char)*(*Ncm6+k+1));
      if (!cm6buf) return NULL;
      
      for (j=0; j<k; j++) cm6buf[j+*Ncm6] = cm6table[(int32_t)subcm6[j]];
      
      *Ncm6 += k;
    }

  for (j=0; j<ndiff; j++) undelta(intbuf, Nint);

  return (cm6buf);
}


/***************************************************************************
 * unpackcm6:
 *
 * Unpack CM6 compressed data into 32-bit integers.
 *
 * cm6buf = input CM6 ASCII
 * Ncm6   = number of input characters
 * intbuf = 32-bit integer buffer (contents destroyed, final is returned)
 * Nint   = number of output integers unpacked
 * Nreq   = number of output integers to unpack, if negative unpack all
 * ndiff  = differencing level (2 is highly recommended)
 *
 * intbuf should be either allocated by the malloc(3) family of
 * functions or set to NULL by the calling program.
 *
 * It is assumed that if Nreq is positive and intbuf is not NULL that
 * intbuf has enough space to hold at least Nreq 32-bit integers.
 *
 * If the original intbuf is big enough to hold Nreq integers it will
 * not be reallocated.
 *
 * Returns the unpacked integers on success and NULL on error.
 ***************************************************************************/

int32_t *
unpackcm6 (const char *cm6buf, int32_t Ncm6,
	   int32_t *intbuf, int32_t *Nint,
	   int32_t Nreq, int32_t ndiff)
{
  static int8_t tablecalculated = 0;
  static int32_t table[256];

  int32_t cin;
  int32_t val, sign=1, cont=0;
  int32_t i, j;
  int32_t Nintbuf = 0;
  int8_t  alloced = 0;

  /* constants for CM6: */
  int32_t SIGNMASK=16;
  int32_t CONTMASK=32;
  int32_t VAL1MASK=15;
  int32_t VAL2MASK=31;
  int32_t CLBTMASK=127;
  int32_t SHIFTVAL=32;
  
  if ( ! tablecalculated )
    {
      for (j=0,i=0; i < 256; i++)
	{
	  if (i=='+') table[i]=j;
	  else if (i=='-') table[i]=j;
	  else if ((i>='0') && (i<='9')) table[i]=j;
	  else if ((i>='a') && (i<='z')) table[i]=j;
	  else if ((i>='A') && (i<='Z')) table[i]=j;
	  else { table[i]=0; continue; }
	  j++;
	}
      
      tablecalculated = 1;
    }
  
  if ( Nreq > 0 )
    Nintbuf = Nreq;
  else
    Nintbuf = 10;
  
  if ( intbuf == NULL )
    {
      intbuf = (int32_t *) malloc (Nintbuf * sizeof(int32_t));
      if (!intbuf) return NULL;
      alloced = 1;
    }
  
  memset (intbuf, 0, Nintbuf * sizeof(int32_t));
  
  for (j=i=0; j < Ncm6; i++)
    {
      if ( Nreq > 0 && i >= Nreq ) break;
      
      /* If intbuf is full allocate more space for 10 more ints */
      if ( i >= Nintbuf )
	{
	  Nintbuf += 10;
	  
	  intbuf = (int32_t *) realloc (intbuf, Nintbuf * sizeof(int32_t));
	  if (!intbuf) return NULL;
	  alloced = 1;
	  
	  memset (intbuf+i, 0, 10 * sizeof(int32_t));
	}

      cin = cm6buf[j++];
      val = table[cin] & CLBTMASK;
      sign = val & SIGNMASK;
      cont = val & CONTMASK;
      intbuf[i] = val & VAL1MASK;

      while (cont)
	{
	  intbuf[i] *= SHIFTVAL;
	  
	  if ( j >= Ncm6 )
	    {
	      fprintf (stderr, "unpackcm6: truncated data\n");
	      break;
	    }
	  
	  cin = cm6buf[j++];
	  val = table[cin] & CLBTMASK;
	  cont = val & CONTMASK;
	  intbuf[i] += val & VAL2MASK;
	}
      
      if (sign) intbuf[i] = -intbuf[i];
      sign=0;
    }
  
  *Nint = i;
  
  /* Make sure intbuf is just the right size, should always truncate */
  if ( alloced && *Nint != Nintbuf )
    {
      intbuf = (int32_t *) realloc (intbuf, (*Nint) * sizeof(int32_t));
      if (!intbuf) return NULL;
    }
  
  if ( Nreq > 0 && *Nint != Nreq )
    fprintf (stderr, "unpackcm6: %d of %d samples unpacked\n", *Nint, Nintbuf);
  
  while (ndiff--) undelta (intbuf, *Nint);
  
  return (intbuf);
}


/***************************************************************************
 * gsechksum:
 *
 * Calculate a checksum using the algorithm specified by GSETT-3 (used
 * in GSE as the value for CHK2 specifier).
 *
 * intbuf  = input 32-bit integers
 * Nint    = number of input integers
 *
 * Returns the checksum.
 ***************************************************************************/
int32_t
gsechksum (int32_t *intbuf, int32_t Nint)
{
  int32_t i, val;
  int32_t cksum=0L;
  int32_t MODULO=100000000;

  for (i=0; i<Nint; i++)
    {
      if ( abs(intbuf[i]) > MODULO ) val=intbuf[i] - MODULO * (intbuf[i] / MODULO);
      else val=intbuf[i];

      cksum+=val;

      if ( abs(cksum) > MODULO )
	cksum= cksum - MODULO * (cksum/MODULO);
    }

  return (abs(cksum));
}
