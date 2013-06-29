#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "fmedian_src.h" /* Takes care of idl_export.h stuff */

/* Quicksort-based (partitioning) k'th element finder
   Based on http://www.mathcs.carleton.edu/courses/course_resources/
   cs227_w96/swanjorg/algorithms.html

   INPUTS:
   A - Array with data (N elements)
   K - The number of the element to be found
   N - The number of elements in A
*/

#define QSRT_k_func2(TYPE)				\
QSRT_k_func2_sig(TYPE)  				\
{							\
  TYPE av;						\
  long i,j,l,r,test;					\
  l = 0;						\
  r = n-1;						\
							\
  while (r > l) {					\
    av = a[r];						\
    i = l - 1;						\
    j = r;						\
							\
    test = 1;						\
    while (test) {					\
      ++i;						\
      while ((a[i] < av) && (i < j)) ++i;		\
      --j;						\
      while ((a[j] > av) && (j > i)) --j;		\
      if (i < j) SWAP2(TYPE,a[i],a[j])			\
	else test = 0;					\
    }							\
							\
    SWAP2(TYPE,a[i],a[r])				\
							\
      /* Terminate if to the k'th ascending point */	\
      if (i == k) return;				\
							\
    /* Select the correct interval */			\
    if (k < i) r = i-1;					\
    else       l = i+1;					\
  }							\
  return;						\
}

QSRT_k_func2(UCHAR)
QSRT_k_func2(IDL_INT)
QSRT_k_func2(IDL_UINT)
QSRT_k_func2(IDL_LONG)
QSRT_k_func2(IDL_ULONG)
QSRT_k_func2(IDL_LONG64)
QSRT_k_func2(IDL_ULONG64)
QSRT_k_func2(float)
QSRT_k_func2(double)

#define FMEDIAN_CORE_FUNC2(TYPE)					\
  FMEDIAN_CORE_FUNC2_SIG(TYPE)	/* See fmedian_src.h */			\
{									\
    /* Half width points */						\
  /* Even sizes (2) rounded down, odd sizes(3) unaffected */		\
  long midx = (xbox - 1)/2;						\
  long midy = (ybox - 1)/2;						\
  									\
  long outx, outy;							\
  									\
  /* printf("Entered FMEDIAN_CORE_FUNC v2\n");	*/			\
  /* If box size given for median filter are out of range, */		\
  /* return original array unmodified.  */				\
  									\
  if ((((xbox > ybox) ? xbox : ybox) < 2) ^ (nx < 1)) {			\
    long i;								\
    for (i=0; i < (nx * ny); ++i) out[i] = in[i];			\
  } else {								\
    									\
    if (find_missing) {							\
      long i;								\
      for (i=0; i < nx*ny; i++)						\
	missing = MIN2(*(in+i),missing);				\
      missing -= 1;							\
    }									\
    									\
    /* Find box with points to use */					\
    /* Careful not to exceed the limits of the input array */		\
    									\
    for (outy=0; outy < ny; ++outy) {	/* Out-array Y loop */		\
      /* Input source box Y top/bottom lim within input boundaries */	\
      long iny_bot = MAX2(outy - midy,0);				\
      long iny_top = MIN2(outy - midy + ybox,ny) - 1;			\
      									\
      for (outx=0; outx < nx; ++outx) {	/* Out-array X loop */		\
	/* Input source box X left/right lim within input boundaries */	\
        long inx_lft = MAX2(outx-midx,0);				\
        long inx_rgt = MIN2(outx - midx + xbox,nx) - 1;			\
									\
        out[outx+nx*outy] = in[outx+nx*outy]; /* Default value */	\
        								\
        /* Store & count points from source box into WORK */		\
        if (always || in[outx+nx*outy]==missing) {			\
	  long np=0; /* Number of non-missing points inside box */	\
	  long iny, inx; /* x/y position, in-array */			\
          for (iny=iny_bot; iny<=iny_top; ++iny) {			\
            for (inx=inx_lft; inx<=inx_rgt; ++inx) {			\
              if (in[inx+nx*iny] != missing) work[np++] = in[inx+nx*iny]; \
	    }								\
	  }								\
	  								\
          /* If zero or one point, leave org (missing or other) */	\
	  if (np==0 || np==1) continue;					\
	  								\
	  /* Special: if two points, use average(!) */			\
	  if ( np == 2) {						\
	    out[outx+nx*outy] = (work[0]+work[1])/2;			\
	    continue;							\
	  }								\
	  /* If more than one point, do sort & pick  midpoint */	\
	  QSRT_k_call2(TYPE,work,np/2,np);				\
	  out[outx+nx*outy] = work[np/2];				\
	}								\
      }									\
    }									\
  }									\
  /* printf("Returning from FMEDIAN_CORE_FUNC2\n"); */			\
}


FMEDIAN_CORE_FUNC2(UCHAR)

FMEDIAN_CORE_FUNC2(IDL_INT)

FMEDIAN_CORE_FUNC2(IDL_UINT)

FMEDIAN_CORE_FUNC2(IDL_LONG)

FMEDIAN_CORE_FUNC2(IDL_ULONG)

FMEDIAN_CORE_FUNC2(IDL_LONG64)

FMEDIAN_CORE_FUNC2(IDL_ULONG64)

FMEDIAN_CORE_FUNC2(float)

FMEDIAN_CORE_FUNC2(double)
