 /* C source for I1MACH */
 /* Note that some values may need changing -- see the comments below. */
#include <stdio.h>
#include <float.h>
#include <limits.h>
#include <math.h>
 
#include "fpreproc/f77name.h"   /* fortran linkage for i1mach, r1mach, d1mach */
 
#define LOCAL_LONG_MAX32 2147483647
 
int F77NAME(i1mach)(int *i)
{
 
       switch(*i){
         case 1:  return 5;    /* standard input  unit -- may need changing */
         case 2:  return 6;    /* standard output unit -- may need changing */
         case 3:  return 7;    /* standard punch  unit -- may need changing */
         case 4:  return 0;    /* standard error  unit -- may need changing */
#if __CRAY
	 case 5:  return 64;
#else
         case 5:  return 32;   /* bits per integer -- may need changing */
#endif
         case 6:  return 1;    /* Fortran 77 value: 1 character */
                               /*    per character storage unit */
         case 7:  return 2;    /* base for integers -- may need changing */
#if __CRAY
	 case 8:  return 63;
#else
         case 8:  return 31;   /* digits of integer base -- may need changing */
#endif
#if __CRAY
	 case 9:  return LONG_MAX;
#else
         case 9:  return LOCAL_LONG_MAX32;
#endif
         case 10: return FLT_RADIX;
         case 11: return FLT_MANT_DIG;
         case 12: return FLT_MIN_EXP;
         case 13: return FLT_MAX_EXP;
 
         case 14: return DBL_MANT_DIG;
         case 15: return DBL_MIN_EXP;
         case 16: return DBL_MAX_EXP;
         }
 
       fprintf(stderr, "invalid argument: i1mach(%ld)\n", *i);
       exit(1);
       return 0; /* for compilers that complain of missing return values */
       }
 
 /* C source for R1MACH */
 
float F77NAME(r1mach)(int *i)
{
       switch(*i){
         case 1: return FLT_MIN;
         case 2: return FLT_MAX;
         case 3: return FLT_EPSILON/FLT_RADIX;
         case 4: return FLT_EPSILON;
         case 5: return (float)log10((double)FLT_RADIX);
         }
 
       fprintf(stderr, "invalid argument: r1mach(%ld)\n", *i);
       exit(1);
       return 0; /* for compilers that complain of missing return values */
       }
 
 /* Standard C source for D1MACH -- remove the * in column 1 */
 
double F77NAME(d1mach)(int *i)
{
       switch(*i){
         case 1: return DBL_MIN;
         case 2: return DBL_MAX;
         case 3: return DBL_EPSILON/FLT_RADIX;
         case 4: return DBL_EPSILON;
         case 5: return log10((double)FLT_RADIX);
         }
       fprintf(stderr, "invalid argument: d1mach(%ld)\n", *i);
       exit(1); return 0; /* some compilers demand return values */
}
