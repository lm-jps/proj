/* code pulled from ana for a standalone despike for radiation hits */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/time.h>
#define ALN2I 1.442695022
#define TINY 1.0e-5
#define ABS(x) ((x)>=0?(x):-(x))
// double systime();  	/* note that systime is defined elsewhere in the ana environment, remove
// 			all refs for standalone (or add func systime) */
 /* 3/17/2010 - interface change to I*4, also work on locations to remove duplicates (some
  pixels are despiked twice when iterated) and do not use output array, modify the input array */
 /* 9/29/2010 - some mods to remove negative values from the data */
 /* 11/20/10 - now want to keep neg values, so remove those mods and make sure
    we don't get too many extra spikes by adjusting absmin (from 1 to 4), we
    already had a test to avoid despiking neg pixels from 9/29/2010 mods */
 /*---------------------------------------------------------------------------*/
int aia_despike(
     float *array, unsigned char *mask, int nx, int ny, float frac, int level, int niter, /* void **outarray, */
     int *badblobs, int sizeofbads,
     int *nspikes, float **oldvalues, float **spikelocs, float **newvalues);
int mask_erode(unsigned char *pbase, unsigned char *qbase, int n, int m);   /* erosion function */
 /*------------------------------------------------------------------------- */
int aia_despike(            /* despike function, modified heavily for AIA */
     float *array, unsigned char *mask, int nx, int ny, float frac, int level, int niter, /* void **outarray, */
     int *badblobs, int sizeofbads,
     int *nspikes, float **oldvalues, float **spikelocs, float **newvalues)
 /* the call is
  int *ptr, *oldvalues, *newvalues;
  int nx, ny, level=7, niter=3, *spikelocs;
  float frac = 0.5;
  unsigned char *mask;
  
  stat = despike(ptr, mask, nx, ny, frac, level, niter,
  &nspikes, &oldvalues, &spikelocs, &newvalues);
  
  note that a NULL passed for mask implies full image (a fake mask is generated internally)
  a 0 for sizeofbads implies no bad pixel correction
  
  note that the nx by ny I*4 array in ptr will be despiked thereby modifying the original,
  the original can be recovered using the offsets in spikelocs and the original values
  in oldvalues, the new values are also available in the newvalues array, these 3 arrays
  are malloc'ed here and the caller is responsible for free'ing them
  
  The bad pixels are represented as bad "blobs" to allow for possible clusters. The badblobs
  array is a null terminated packed list of these, the size of the list must also be
  passed in sizeofbads to avoid any overflows (and allow for a quick check for none),
  the format of the list is [nperim, perim_offsets, npix, pix_offsets]  in series for
  each "blob", nperim is the # of perimeter points (offsets in ptr) to use for averaging
  to obtain a value to substitute for the bad pixels. perim_offsets is the list
  of nperim offsets. npix is the count of bad pixels to be replaced by this average and
  pix_offsets is the npix offsets. Note that all bad pixels are replaced by the same value
  (no gradients or such) computed from the average of the perimeter. The length of each
  in the the list is nperim + npix + 2. The code just steps through.
 */
 {
 int	iq, result_sym, n, m, sum, nc;
 int	lognb2, jj, jc;
 float	absmin = 4.0;
 int	nxc, nyc, offset;
 int	save_niter, ntotal, npix, nslocal, ndups, itercount;
 float	fsum, cfrac, tq, fdif, *lastiterendaddr;
 float  *p, *p1, *p2, *p3, *ptr, *base, *pps, arr[16], *ss;
 float	*spikeold, *oldptr;
 unsigned char	*mptr, *mp, mskv, *mskbase, *eroded, *eroded2;
 float	*spikebufadd, *sptradd;
 float	*spikenew, *newptr;
 float	*spikestarts[20], *spikeends[20];
 //double	t1, t2, t3, t4;
 //t1 = systime();
 lognb2 = (log((double) 16)*ALN2I+TINY);  /* needed for median sort */

 /* always an int array for this routine */
 if (nx < 5 || ny < 5 ) {
 	printf("dimensions must be 5 or greater, nx, ny = %d %d\n", nx,ny);
 	return -1; }
 /* bad pixel correction (if any), uses blob list as noted above */
 if (sizeofbads > 0 && badblobs != 0 ) {
   /* correction done in place, if badblobs not defined correctly, a memory access violation
   could occur */
   int  countdown = sizeofbads;
   int	*pblob, nperim, npix, xq, offset;
   float  acc;
   pblob = badblobs;
   while (1) {
     acc = 0.0;
     nperim = *pblob++;
     if (nperim <= 0 || nperim > countdown) break;
     n = nperim;
     /* this takes the average of the perimeter values, a median might be better */
     while (n--)  { offset = *pblob++;  acc +=  *(array + offset); }
     xq = rintf(acc/(float) nperim);
     npix = *pblob++;
     if (npix <= 0 || npix > countdown) break;
     n = npix;
     while (n--)  { offset = *pblob++;  *(array + offset) = xq; }
     countdown = countdown - nperim - npix - 2;
     if (countdown <= 0) break;
   }
 }

 /* flip image */
 {
    int i;
    float *toprow, *botrow, *tmprow;
    botrow = array;
    toprow = botrow + 4095*4096;
    tmprow = (float *) malloc(sizeof(int)*4096);
    for (i = 0; i < 2048; i++) {
      memcpy(tmprow, toprow, 16384);
      memcpy(toprow, botrow, 16384);
      memcpy(botrow, tmprow, 16384);
      botrow += 4096;
      toprow -= 4096;
    }
    free(tmprow);
 }
 /* if no mask, we actually make one full of 1's, partly for testing */
 eroded = malloc(nx*ny*sizeof(unsigned char));
 if (!eroded) { printf("malloc error in local mask copy\n");  return 1; }
 mp = eroded;
 /* unfortunately, we need to check array for 0x80000000 values and avoid using
 these for despike decisions and we must not modify them. This is done by  zeroing
 the mask for the corresponding points before erosion */
 p = array;
 /* 9/29/2010 - also catch any other negative values while we are looping and zero them */
 if (mask == 0) {
   int	nq = nx * ny;
   while (nq--) {
     /* 9/29/2010 - change here to remove negative values (convert to 0) */
     /* 11/20/10 - try leaving the negative values in */
     //if (*p == 0x80000000) { *mp++ = 0; } else  { *mp++ = 1;  if (*p < 0) *p = 0; }
     if (isnan(*p)) { *mp++ = 0; } else  { *mp++ = 1; }
     p++;
   }
   /* mask gets set to "eroded" further down */
 } else {
   /* have a real mask, need to erode it after checking for 0x80000000's */
   int	nq = nx * ny;
   mp = mask;
   while (nq--) {
     /* 9/29/2010 - change here to remove negative values (convert to 0) */
     /* 11/20/10 - try leaving the negative values in */
     //if (*p == 0x80000000) { *mp++ = 0; } else  { if (*p < 0) *p = 0; }
     if (isnan(*p)) { *mp++ = 0; }
     p++;
   }
   eroded2 = malloc(nx*ny*sizeof(unsigned char));
   if (!eroded2) { printf("malloc error in local mask copy\n");  return 1; }
   mask_erode(mask, eroded2, nx, ny);
   mask_erode(eroded2, eroded, nx, ny);
   /* leaves final result in eroded, so free eroded2 now and redefine mask */
   free(eroded2);
 }
 mask = mp = eroded;
 ptr = base = array;
 /* and the mask */
 mptr = mskbase = mask;
 cfrac = 1.0/(1.0 + frac);	/* 3/2/2010 changed from 1.0 - frac which had less range */
 nc = ntotal = nslocal = 0;
 niter = ABS(niter);
 /* 4/7/2010 - add a check for niter < 0 to avoid runaways */
 if (niter > 20 || niter < 0) { printf("DESPIKE - error, excessive # of iterations = %d\n",
 	niter);  return 1; }

 /* add internal iteration 10/8/98 */
 save_niter = niter;
 spikeold = oldptr = malloc(nx*ny*sizeof(float));  /* that allows every pixel */
 spikenew = newptr = malloc(nx*ny*sizeof(float));  /* that allows every pixel */
 spikebufadd = sptradd = malloc(nx*ny*sizeof(float));  /* that allows every pixel */
 lastiterendaddr = spikebufadd;
 ndups = 0;
 //t2 = systime();
 npix = 0;
 itercount = 0;
 spikestarts[itercount] = spikebufadd;
 /* note that niter = 0 should not do any despike */
 while (niter--) {
   ptr = array + 2*nx;
   /* and the mask */
   mptr = mask + 2*nx;
   /* skip the outer edges */
   m = ny-4;
   while (m--) {
     /* skip the outer edges */
     p = ptr;
     p += 2;
     /* and the mask */
     mp = mptr + 2;
     n = nx-4;
     p2 = p - 1;
     p1 = p2 - nx;
     p3 = p2 + nx;
     while (n--) {
      npix++;
      /* watch out for 0's, lots in test images */
      /* 3/22/2010 - mask check finally added for each position */
      /* 9/29/2010 - avoid processing negative values as well as zeroes */
      if ( (*p > 0)  && *mp ) {
	/* add the 8 around this point */
	tq = (cfrac * (float) *p);
	sum = *p1 + *(p1+1) + *(p1+2) + *p2 + *(p2+2) + *p3 + *(p3+1) + *(p3+2);
	/* note the absmin term, this is to avoid problems with data that has large
	zones of zeroes, we then get too many hits and erosion of edges, mostly a problem
	with artifical data */
	fsum = (float) sum * 0.125;
	/* now the test */
	if ( (fsum < tq)  && ((tq-fsum) > absmin) ) {  /* we have a bad one, zap it */
	  nc++;
	  /* load up sort array and find the desired one */
	  ss = arr;	pps = p - 2*nx -2;
	  *ss++ = *pps++; *ss++ = *pps++; *ss++ = *pps++; *ss++ = *pps++; *ss++ = *pps++;
	  *ss++ = *(p - nx -2);  *ss++ = *(p - nx +2);
	  *ss++ = *(p -2);  *ss++ = *(p +2);
	  *ss++ = *(p + nx -2);  *ss++ = *(p + nx +2);
	  pps = p + 2 *nx - 2;
	  *ss++ = *pps++; *ss++ = *pps++; *ss++ = *pps++; *ss++ = *pps++; *ss++ = *pps++;
	  /* this is basically a median filter (adjustable using level value) to find
	  a "good" value from the surroundings, we do a median rather than an average to
	  try to avoid some effects of spike overflow on the CCD */
	  /* a built in sorter here */
	  { int nn,m,j,i,n=16;
	     int t;
	     m=n;
	     for (nn=1;nn<=lognb2;nn++) {
		m >>= 1;
	       for (j=m;j<n;j++) {
		  i=j-m;
		  t=arr[j];
		  while (i >= 0 && arr[i] > t) {
		    arr[i+m]=arr[i];
		    i -= m;
		  }
		  arr[i+m]=t;
	      }
	    }
	  }
	  /* 3/11/2010 - log into various arrays */
	  *oldptr++ = *p;  offset = p - base;
	  //if (nslocal == 0) printf("sptradd, offset = %d, %d\n", sptradd, offset);
	  *sptradd++ = offset;
	  nslocal++;  /* yet another spike counter */
	  /* now get the indicated one using level, also save into spikenew */
	  *newptr++ = arr[level];
	 } 
       }
       p++;  mp++;  p1++;  p2++;  p3++;

     }
     /* skip the outer edges */
     //p = p + 2;
     ptr = ptr + nx;
     /* and the mask */
     //mp = mp + 2;
     mptr = mptr + nx;
   }
   /* skip the last 2 rows */
   /* end of this iteration, reconfigure for next one, note that for
   a single iteration we don't need out2 */
//   printf("despike got %d spikes in iteration %d out of %d pixels\n", nc, itercount, npix);
//   printf("nslocal = %d\n", nslocal);
   ntotal += nc;  nc = 0;
   npix = 0;
   if (niter < (save_niter - 1)) {
     float	*ip, *jp, *iplast = spikebufadd;
     int	ioff, joff,  delta1, delta2, previter;
     /* check the range spikebufadd to lastiterendaddr against lastiterendaddr+1 to sptradd-1,
     the idea here is to check if any of our new spikes have the same location from a previous
     iteration, this isn't a n*m operation because each set is monotonic, note however that
     we have to scan each iteration set separately (which is annoying) */
     float *jprevstart[20];
     spikestarts[itercount] = spikeends[itercount-1] + 1;
     for (previter=0; previter<itercount; previter++) { jprevstart[previter] = spikestarts[previter]; }
     for (ip=spikestarts[itercount]; ip<sptradd; ip++) {
       ioff = *ip;
       /* now compare with all the previous ones */
       for (previter=0; previter<itercount; previter++) {

	 for (jp=jprevstart[previter]; jp<=spikeends[previter]; jp++) {
           joff = *jp;
	   if (joff == ioff) {
	     /* we want the new instance of this location to be the one used, it needs
	     the earlier old value, it already has the latest new value */
	     delta2 = jp - spikebufadd;
	     delta1 = ip - spikebufadd;
	     *(spikeold + delta1) = *(spikeold + delta2);
	     /* and zap the original location */
	     *jp = -1;
	     ndups++;
	     break;
	   }
	   if (joff > ioff) {
	     /* because both are monotonic, can't be any more now */
	     
	     break;
	   }
	 }
	 /* reset the start to avoid wasting time on check for next jp since both are monotonic */
	 jprevstart[previter] = jp;
       }
     }
     //printf("# of dups = %d\n", ndups);
   }
   spikeends[itercount] = sptradd - 1;
//    printf("setting spikestarts[itercount], spikeends[itercount] = %d, %d\n", spikestarts[itercount],
//    	spikeends[itercount]);
   /* and we now set the new values from this iteration in the input array */
   {
     int ioff, nq;
     float *snew, *ip;;
     ip=spikestarts[itercount];
     snew = spikenew + (ip - spikebufadd);
     while ( ip<sptradd ) {
       ioff = *ip++;
       array[ioff] = *snew++;
     }
   }
   itercount++;
 }
 /* wrap up the spike log */
 *nspikes = nslocal;
// printf("# of spikes logged = %d, dups = %d\n", nslocal, ndups);
 *nspikes = nslocal - ndups;
// printf("# of unique spikes logged = %d\n", *nspikes);
 n = (nslocal - ndups)*sizeof(int);
 /* the duplicates are tagged, copy only the unique ones for export */
 if (n > 0) {
   int	i;
   float	*p1, *p2, *p3, *q1, *q2, *q3;
   *oldvalues = p1 = malloc(n);
   q1 = spikeold;
   *newvalues = p2 = malloc(n);
   q2 = spikenew;
   *spikelocs = p3 = malloc(n);
   q3 = spikebufadd;
   for (i=0;i<nslocal;i++) {
     if (*q3 > 0) { *p1++ = *q1; *p2++ = *q2; *p3++ = *q3; }
     q1++;  q2++;  q3++;
   }
 } else {
   *oldvalues = 0; *spikelocs = 0;  *newvalues = 0;
 }
 /* free up (reverse order) */
 free(spikebufadd); free(spikenew); free(spikeold); free(eroded);
 
// printf("despike got %d spikes in %d iterations\n", ntotal, save_niter);
// t3 = systime();

// printf("despike times, setup %10.6fs, iters %10.6fs, total %10.6fs\n", t2-t1,t3-t2,t3-t1);
 return 0;
 }
 /*------------------------------------------------------------------------- */
 /* a standalone erode function specialized for eroding the AIA masks prior to
 despiking */
 /*------------------------------------------------------------------------- */
int mask_erode(unsigned char *pbase, unsigned char *qbase, int n, int m)   /* erosion function */
 {
 unsigned char *p, *q, *p_endline, *qabove, *qbelow;
 int    iq, mq, type;
 double t1, t2, t3, t4;
// t1 = systime();
 if ( n < 3 || m < 3 ) {
        printf("ERODE: array too small, nx, ny = %d %d\n", n, m);
        return -1; }

 bcopy(pbase, qbase, n*m); 
 p = pbase;
 q = qbase;

 /* the edges and corners are done as special cases */
 /* first corner */
 if (*p == 0) {  /* zap the 3 neighbors */
 *q++;  *q = 0; q = qbase + n;  *q++ = 0; *q = 0;}
 /* starting row */
 p_endline = pbase + n - 2;
 p++;
//  t2 = systime();
 while (1) {
  if (*p == 0) {
  /* got a hit, this means we need to clear 6 pixels in the output image */
  q = (qbase - pbase - 1) + p;   /* this is the q just before the corresponding p */
  qabove = q + n;
  *q++ = 0;             q++;       *q++ = 0;
  *qabove++ = 0;        *qabove++ = 0;  *qabove++ = 0;
  if (p >= p_endline) break;
  p++;
  /* if we continue to get consecutive hits, just need to set next 2 */
  while (*p == 0) {
   if (p >= p_endline) break;
   *q++ = 0;    *qabove++ = 0;
   p++; }
  }
  if (p >= p_endline) break;
  p++;
  }
  p++;
  /* last point in starting row, set independent of previous so we may
  set some pixels twice */
  if (*p == 0) {  
  q = (qbase - pbase - 1) + p;   /* this is the q just before the corresponding p */
  qabove = q + n;
  *q = 0;
  *qabove++ = 0;        *qabove = 0;
  }
  p++;


 /* central area */ 
 /* now the interior rows */
 mq = m - 2;
 while (mq-- > 0) {
 /* start row, not top or bottom */
 /* left hand edge */
 if (*p == 0) {      /* set 6 points */
  q = (qbase - pbase - 1) + 1 + p;        /* q at p */
  qabove = q + n;       qbelow = q - n; 
  *q++;             *q = 0;
  *qabove++ = 0;        *qabove = 0;
  *qbelow++ = 0;        *qbelow = 0;
 }
 p_endline = p + n - 2;
 p++;

 /* done with left edge, now the middle */
 
 while (1) {
  if (*p == 0) {
  /* got a hit, this means we need to clear 8 pixels in the output image */
  q = (qbase - pbase - 1) + p;
  qabove = q + n;       qbelow = q - n; 
  *q++ = 0;             q++;       *q++ = 0;
  *qabove++ = 0;        *qabove++ = 0;  *qabove++ = 0;
  *qbelow++ = 0;        *qbelow++ = 0;  *qbelow++ = 0;
  if (p >= p_endline) break;
  p++;
  /* if we continue to get consecutive hits, just need to clear next 3 */
  while (*p == 0) {
   if (p >= p_endline) break;
   *q++ = 0;    *qabove++ = 0;  *qbelow++ = 0;
   p++; }
  }
  if (p >= p_endline) break;
  p++;
  }
  p++;


 /* the last point in row */
  if (*p == 0) {  
  q = (qbase - pbase - 1) + p;   /* this is the q just before the corresponding p */
  qabove = q + n;       qbelow = q - n; 
  *q = 0;
  *qabove++ = 0;        *qabove++ = 0;
  *qbelow++ = 0;        *qbelow++ = 0;
  }
  p++;

 }

 /* at last, the last row */
 /* left hand edge */
 if (*p == 0) {      /* set 4 points */
  q = (qbase - pbase - 1) + 1 + p;        /* q at p */
  qbelow = q - n;       
  q++;              *q = 0;
  *qbelow++ = 0;         *qbelow = 0;
 }
 p_endline = p + n - 2;
 p++;

 /* now the middle of the last row */
 
 while (1) {
  if (*p == 0) {
  /* got a hit, this means we need to clear 9 pixels in the output image */
  q = (qbase - pbase - 1) + p;
  qbelow = q - n;       
  *q++ = 0;             q++;       *q++ = 0;
  *qbelow++ = 0;        *qbelow++ = 0;  *qbelow++ = 0;
  if (p >= p_endline) break;
  p++;
  /* if we continue to get consecutive hits, just need to set next 3 */
  while (*p == 0) {
  if (p >= p_endline) break;
  *q++ = 0;     *qbelow++ = 0;  p++; }
  }
  if (p >= p_endline) break;
  p++;
  }
  p++;
  /* printf("pbase, p_endline, p = %#x, %#x, %#x\n", pbase, p_endline, p); */
 /* the last point in last row */
  /* printf("pbase, p = %#x, %#x\n", pbase, p); */
  if (*p == 0) {  
  q = (qbase - pbase - 1) + p;   /* this is the q just before the corresponding p */
  qbelow = q - n;       
  *q = 0;
  *qbelow++ = 0;        *qbelow = 0;
  }
//  t3 = systime();
//   printf("AIA erode total time = %10.2f ms, for setup = %10.2f ms, part 2 = %10.2f ms\n",
//     1.E3*(t3-t1), 1.E3*(t2-t1), 1.E3*(t3-t2));
 return 0;
 }
