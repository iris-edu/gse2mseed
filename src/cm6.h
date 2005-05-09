
#ifndef CM6_H
#define CM6_H

#ifdef __cplusplus
extern "C" {
#endif

#include <inttypes.h>

/* Table of all possible CM6 characters */
static char cm6table[64]= { '+','-','0','1','2','3','4','5','6','7',
			    '8','9','A','B','C','D','E','F','G','H',
			    'I','J','K','L','M','N','O','P','Q','R',
			    'S','T','U','V','W','X','Y','Z','a','b',
			    'c','d','e','f','g','h','i','j','k','l',
			    'm','n','o','p','q','r','s','t','u','v',
			    'w','x','y','z' };
  
char *packcm6 (int32_t *intbuf, int32_t Nint, char *cm6buf,
	       int32_t *Ncm6, int32_t ndiff);
  
int32_t *unpackcm6 (const char *cm6buf, int32_t Ncm6, int32_t *intbuf,
		    int32_t *Nint, int32_t Nreq, int32_t ndiff);

int32_t gsechksum (int32_t *intbuf, int32_t Nint);

#ifdef __cplusplus
}
#endif

#endif
