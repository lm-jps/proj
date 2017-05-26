#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <ctype.h>

/*------------------------------------------------------------------------- */
extern int verbose;
static float lastr1, lastr2, lastyc, lastxc;
float sdisk_xc=2047.5, sdisk_yc=2047.5, sdisk_r, chicircle;
static unsigned char *mask, *mask2, *masks;
static int lastnx, lastny;
static float gsmooth_width = 2.3;	/* limit for gsmooth kernel, good for .005, use 4 for 1.e-7 */
static long *jstarts, *jvalues, *jstarts_s, *jvalues_s, *jstarts2, *jvalues2;
static int num_inann, num_inann2, num_inanns, num_in2use;
static int *in2use, *thann;
static  int  *ixnann, *iynann;
static long *inann, *inanns, *inann2;
static float *annth, *annthr;
int scounts ;
typedef unsigned char byte;
#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define ABS(x) ((x)>=0?(x):-(x))
#define ALN2I 1.442695022
#define TINY 1.0e-5
#define MAXNITER 100
/* arrays showing convergence, length is MAXNITER */
float xcenters[MAXNITER], ycenters[MAXNITER], radiai[MAXNITER];
long	pcounts[MAXNITER], nextcenter;
int max4histogram = 10000;   /* sized here for radius calculations */
void mmq(float *x, int n, float *max, float *min);
void indexf(int n, float ra[],int indx[]);
long getindexofmax(int *x, int n);
int pit(float *qxbase, float *qybase, int nxq, int npow, double *a, double *cfbase, double *fbase );
float parabolicmax5(int *x, int nx);
void sobel_index(float *x, int nx, int ny, long *iq, int nindex, float **grad, float **ang);
void mm_f(float *x, int n, float *max, float *min);
void mm_int(int *x, int n, int *max, int *min);
void indexf(int n, float ra[],int indx[]);
int circle_fit_lsq(float *x, float *y, float *w, int n, float *xc_out, float *yc_out, float *r_out, float *chi_out);
int key_locations(int *array, int n, int *nv, int **counts, int **values, int ***indexptrs);
long hist_max_xvalue(float *x, int n);
int *hist_array(float *x, int n, int *nh, int *hmin);
int hist_max_freq(float *x, int n);
void getmin9(int *p, int ix, int iy, int nx, float *x0, float *y0);
void zeroinonlimb(float *xlimb, float *ylimb, int nin, float dr, float *g);
void  minihough(int nc, float *xlimb, float *ylimb, int niter);
void decomp(double *x, int n, int nd);
void solve(double *a, double *b, int n, int nd);
int pit(float *qxbase, float *qybase, int nxq, int npow, double *a, double *cfbase, double *fbase );
float parabolicmax5(int *x, int nx);
void printarray(float *x, int n);
void printintarray(int *x, int n);
#define EPSILON 1.e-10
double systime();				/* internal systime */
/*------------------------------------------------------------------------------------------------*/
double systime()				/* internal systime */
{
  double t;
  struct timeval tp;
  struct timezone tzp;
  gettimeofday(&tp, &tzp);
  t = (double) tp.tv_sec + .000001* (double) tp.tv_usec;
  return t;
}

float rsun_offset(int wavelength)
{
  switch (wavelength) {
    case  94: return 7.3;
    case 131: return 7.6;
    case 171: return 9.3;
    case 193: return 9.9;
    case 211: return 10.4;
    case 304: return 18.5;
    case 335: return 11.0;
    case 1600: return 1.4;
    case 1700: return 0.2;
    case 4500: return -1.2;
    default: return FP_NAN;
  }
}
/*--------------------------------------------------------------------------*/
void getmin9(int *p, int ix, int iy, int nx, float *x0, float *y0)
 {
 /* modified from the ana internal getmin, computes local minimum for a 3x3 subarray in p using x index
 ix, y index iy, and x dimension nx */
 float	f11, f12, f13, f21, f22, f23, f31, f32, f33;
 float	toler, fx, fy, t, fxx, fyy, fxy;
 int *pp;
 long off;
 /* find the min, p points to a 3x3 array */
 /* assuming ix and iy are not on edges */
 off = (iy-1)*nx;
 pp = p + ix + off -1;
 f11 = (float) *pp++;   f21 = (float) *pp++;	f31 = (float) *pp;
 pp = pp + nx - 2;
 f12 = (float) *pp++;	f22 = (float) *pp++;	f32 = (float) *pp;
 pp = pp + nx - 2;
 f13 = (float) *pp++;	f23 = (float) *pp++;	f33 = (float) *pp;
 
 fx = 0.5 * ( f32 - f12 );	fy = 0.5 * ( f23 - f21 );
 t = 2.* ( f22 );
 fxx =  f32 + f12 - t;	fyy = f23 + f21 - t;
 /* find in which quadrant the minimum lies */
 if (f33 < f11)  { if (f33 < f31) {
	 if (f33 < f13) fxy = f33+f22-f32-f23; else fxy = f23+f12-f22-f13; }
	 else			 {
	 if (f31 < f13) fxy = f32+f21-f31-f22; else fxy = f23+f12-f22-f13; }
  } else   { if (f11 < f31) {
	 if (f11 < f13) fxy = f22+f11-f21-f12; else fxy = f23+f12-f22-f13; }
	 else			 {
	 if (f31 < f13) fxy = f32+f21-f31-f22; else fxy = f23+f12-f22-f13; }
 }
 t = -1./(fxx *fyy - fxy *fxy);
 *x0 = t * (fx * fyy - fy * fxy);
 *y0 = t * (fy * fxx - fx * fxy);
 if ( ABS(*x0) >= 0.75 || ABS(*x0) >= 0.75 ) { *x0 = -fx/fxx; *y0 = -fy/fyy; }
 return;
 }
 /*------------------------------------------------------------------------- */
void zeroinonlimb(float *xlimb, float *ylimb, int nin, float dr, float *g)
 {
 /* suppose we have a good circle and find points near it from our set in inlimb
 and fit these to make a better circle (we hope), use value of g as weight
 the currrent circle fit is assumed in globals sdisk_xc, sdisk_yc, and sdisk_r and get updated */
 int	i, j, n;
 long	*index;
 float  *x, *y, *w;
 float  rdlimb;
 /* we have to select a subset of xlimb, ylimb, g based on distance from circle */
 n = 0;
 index = malloc(nin * sizeof(long));
 i = 0;
 for (j=0;j<nin;j++) {
   rdlimb = hypotf(xlimb[j]-sdisk_xc, ylimb[j]-sdisk_yc) - sdisk_r;
   if ( fabsf(rdlimb) < dr) { index[i] = j; i++; }
 }
 n = i;  /* the # that passed the test */
 x = malloc(n * sizeof(float));   y = malloc(n * sizeof(float));
 w = malloc(n * sizeof(float));
 for (i=0;i<n;i++) {
   j = index[i];
   x[i] = xlimb[j];  y[i] = ylimb[j];  w[i] = g[j];
 } 
 circle_fit_lsq(x, y, w, n, &sdisk_xc, &sdisk_yc, &sdisk_r, &chicircle);
 if (verbose>2) printf("zeroinonlimb: sdisk_xc, sdisk_yc, sdisk_r = %g, %g, %g\n",sdisk_xc, sdisk_yc, sdisk_r);
 free(index);
 free(x);  free(y);  free(w);
 }
 /*------------------------------------------------------------------------- */
#define NXC 11
#define NYC 11
#define NRC 11
void  minihough(int nc, float *xlimb, float *ylimb, int niter)
 {
 /* assumes current best fit in globals, gets all within 1 pixel in r from center
  as our measure over a range of centers and r
  xlimb, ylimb is a subset of indices from the original image that are points we
  consider, usually the results of an edge selection
 */
 float xq, yq, xoffset, yoffset, *rdlimb, dx, dy, dr, xc1, yc1, rx, dxprev, dyprev;
 int mhcount[NXC * NYC], k, ixc, iyc, irc, ic, inx, ix, iy, *hr, nh, hmin, rmax, iq, ndone;
 double	t0, t1, t2, t3, t4, t1s, tstart, tend, dt1, dt2;
 tstart = systime();
 rdlimb = (float *) malloc(nc * sizeof(float));
 dxprev = dyprev = 5.0;
 if (verbose>2) printf("minihough, nc = %d\n", nc);
 if (niter > MAXNITER) {
   if (verbose>1) printf("iteration count request (%d) exceeds MAXNITER (%d)\n", niter, MAXNITER);
   niter = MAXNITER; }
 
 dt1 = dt2 = 0.0;
 for (k=0;k<niter;k++) {
   t0 = systime();
   /* the currrent circle fit is assumed in globals sdisk_xc, sdisk_yc, and sdisk_r and get updated */
   xoffset = sdisk_xc - (float) (NXC/2);
   yoffset = sdisk_yc - (float) (NYC/2);
   ic = 0;
   rmax = inx = 0;
   for (iyc=0;iyc<NYC;iyc++) {
     /* first we vary the center over a 2-D range and look at all the values of r that result */
     dy = (float) iyc + yoffset;
     for (ixc=0;ixc<NXC;ixc++) {
       dx = (float) ixc + xoffset;
       /* compute the nc rdlimb values */
       for (irc=0;irc<nc;irc++) {
         xq = xlimb[irc] - dx;
	 yq = ylimb[irc] - dy;
	 rdlimb[irc] = hypotf(xq, yq) - sdisk_r;
       }
       /* do a histogram and get the max count aka frequency */
       iq = hist_max_freq(rdlimb, nc);
       mhcount[ic] = -iq;
       if (iq > rmax) { rmax = iq;  inx = ic; }
       ic++;
       //printf("ixc, iyc, iq, rmax = %d, %d, %d, %d\n", ixc, iyc, iq, rmax);
     }
   }
   t1 = systime();
   dt1 += t1 - t0;
   
   //printf("inx = %d\n", inx);
   /* inx is index of the max in mhcount, get ix and iy, we use rmax for pcounts */
   ix = inx%NXC;  iy = inx/NXC;
   dx = ix - (float) (NXC/2);  /*offset from previous */
   dy = iy - (float) (NYC/2);
   //printf("ix, iy, dx, dy = %d, %d, %g, %g\n", ix, iy, dx, dy);
   /* check if on the edge and bail out if so */
   if (ix <= 0 || ix >= (NXC-1)) {
      fprintf(stderr, "x center max at edge, ix = %d\n", ix); return;   }
   if (iy <= 0 || iy >= (NYC-1)) {
      fprintf(stderr, "y center max at edge, iy = %d\n", iy); return;   }
   /* get the new center in the image from these */
   /* call getmin9 */
   getmin9(mhcount, ix, iy, NXC, &xc1, &yc1);
   //printf("xc1, yc1 = %g, %g\n", xc1, yc1);
   /* damp */
   dx = (dx + xc1) * 0.5;
   dy = (dy + yc1) * 0.5;
   /* also prevent dx and dy from growing */
//    if (fabsf(dx) > dxprev) { if (dx < 0) dx = -dxprev * 0.75; else dx = dxprev * 0.75; }
//    if (fabsf(dy) > dyprev) { if (dy < 0) dx = -dyprev * 0.75; else dy = dyprev * 0.75; }
   dxprev = fabsf(dx);
   dyprev = fabsf(dy);
   if (verbose>4) printf("dx, dy = %g, %g\n", dx, dy);
   /* and now use the dx, dy to re-compute the nc rdlimb values */
   for (irc=0;irc<nc;irc++) {
     xq = xlimb[irc] - dx - sdisk_xc;
     yq = ylimb[irc] - dy - sdisk_yc;
     rdlimb[irc] = hypotf(xq, yq) - sdisk_r;
   }
   /* do a histogram and get the max count, these are dr's and hmin is generally negative */
   hr = hist_array(rdlimb, nc, &nh, &hmin);
   rx = parabolicmax5(hr, nh) + (float) hmin;
   free(hr);
   dr = rx * 0.5;  /* damp */
   if (verbose>4) printf("Hough dx, dy, dr = %g, %g, %g\n", dx, dy, dr);
   sdisk_xc += dx;
   sdisk_yc += dy;
   sdisk_r += dr;
   xcenters[nextcenter] = sdisk_xc; ycenters[nextcenter] = sdisk_yc;
   radiai[nextcenter] = sdisk_r; pcounts[nextcenter] = rmax;  nextcenter++;
   /* check if we are done */
   if ((dx*dx + dy*dy) < 0.0002 && fabsf(dr) < 0.01) {
     if (verbose>2) printf("close enough, k = %d\n", k);  break;
   }
   t2 = systime();
   dt2 += t2 - t1;
 }
 ndone = nextcenter;
 /* print out the iteration results while debugging */
 if (verbose>3) {
   printf("$xcenters =");
   printarray(xcenters, ndone);

   printf("$ycenters =");
   printarray(ycenters, ndone);

   printf("$radiai =");
   printarray(radiai, ndone);

   printf("$pcounts =");
   for (k=0;k<ndone;k++) {
     if (k%6 == 0) printf("\n");
     printf("%12.5e", (float) pcounts[k]);
   }
   printf("\n");
   printf("minihough times %8.3f %8.3f\n", dt1, dt2);
   tend = systime();
   printf("total minihough time %8.3f\n", tend - tstart);
 }
 free(rdlimb);
 }
 /*------------------------------------------------------------------------- */
void printarray(float *x, int n)
 {
 /* print out array similar to ana for checks */
 int k;
 printf("%d count array\n", n);
 for (k=0;k<n;k++) {
   if (k%6 == 0) printf("\n");
   printf("%13.5e", x[k]);
 }
 printf("\n");
 }
 /*------------------------------------------------------------------------- */
void printintarray(int *x, int n)
 {
 /* print out array similar to ana for checks */
 int k;
 printf("%d count array\n", n);
 for (k=0;k<n;k++) {
   if (k%8 == 0) printf("\n");
   printf("%10d", x[k]);
 }
 printf("\n");
 }
 /*------------------------------------------------------------------------- */
int circle_fit_lsq(float *x, float *y, float *w, int n, float *xc_out, float *yc_out, float *r_out, float *chi_out)
 {
 /* 2016 - modified from the ana internal function */
 /* 4/11/96 - adapted from the C++ code in the Gismo package which in turn
 was taken from something that was taken from PREPMORT file by CVERT
 on 1 Mar 1991 which is after the method of chernov and ososkov in  computer
 physics communications 33 (1984) 329-333. */
 /* return convention same as Gismo version:
 	0	OK
	1	< 3 points (but not checked here)
	2	no convergence, returns non-iterative solution
	3	singular determinant
	4	rsquared < 0 (imaginary radius!)
 */

 double xxyy,
        sx,sy,sw,swinv,                     /*sums*/
        sxx8,syy8,sxy8,sxxy8,syxy8,sxxyy8,  /*more sums*/
        xt,xtp,top,bot,
        x8,y8,w8,x08,y08,rsq,chisq8,
        f,g,h,p,q,t,                        /*see the paper*/
        a0,b0,c0,d0,det,
        gamma, gamma0,g0inv;
 
 float	*xp, *yp, *wp;
 int i, nptw, iflag = 0;
 int	j, k, nn;
 
 /*  To improve accuracy, we do all calculation in a shifted coordinate*/
 /*  system where the weighted sum of the xi and the yi are both zero  */
 /*  this seems to work even in single precision, but I leave it in    */
 /*  double just for clean living.                                     */
  
 /*  Form the sums, which are used to shift the coordinate system      */
 
 sx = 0.0; sy = 0.0; sw = 0.0; nptw = 0;

 nn = n;
 xp = x;	yp = y;		wp = w;
 sx = sy = sw = 0.0;
 while (nn--) { w8 = *wp++;  sw += w8;  sx += w8 * *xp++;  sy += w8 * *yp++; }

 swinv = 1./sw;	
 sx = sx*swinv; sy = sy*swinv;		/* normalize */
 
 /*  These are the sums needed in the calculation                      */
 
 sxx8 = syy8 = sxy8 = sxxy8 = syxy8 = sxxyy8 = 0.0;
 nn = n;
 xp = x;	yp = y;		wp = w;
 while (nn--) {
   w8 = *wp++;   x8 = *xp++ - sx;	y8 = *yp++ - sy;
   xxyy    =  x8*x8 + y8*y8;		sxx8   +=  x8*x8*w8;
   syy8   +=  y8*y8*w8;			sxy8   +=  x8*y8*w8;
   sxxy8  +=  x8*xxyy*w8;		syxy8  +=  y8*xxyy*w8;
   sxxyy8 +=  xxyy*xxyy*w8;
 }

 f = swinv*(3.*sxx8 + syy8);
 g = swinv*(sxx8 + 3.*syy8);
 h = 2.*swinv*sxy8;
 p = swinv*sxxy8;
 q = swinv*syxy8;
 t = swinv*sxxyy8;

 /*  The initial guess gamma comes from the non-iterative fast 'fit' */
 
 gamma0 = swinv*(sxx8 + syy8);
 g0inv  = 1./gamma0;
 
 a0 = -4.;
 b0 = (f*g - t - h*h)*g0inv*g0inv;
 c0 = (t*(f + g) - 2.*(p*p + q*q))*g0inv*g0inv*g0inv;
 d0 = (t*(h*h - f*g) + 2.*(p*p*g + q*q*f) - 4.*p*q*h)*g0inv*g0inv*g0inv*g0inv;
 
 /*  Here we get the root by newton's method                         */
 
 xt = 1.;
 for(i=  0; i < 10; i++)
   {top = xt*(xt*(xt*(xt + a0) + b0) + c0) + d0;
    bot = xt*(xt*(4.*xt + 3.*a0) + 2.*b0) + c0;
    xtp = xt -top/bot;
    if(fabs(xtp-xt) < EPSILON) break;
    xt = xtp;
   }
 if (i >= 10) { iflag = 2;	xtp = 1.0; } /* no convergence */
 /* for the no convergence case, we use the original gamma0 above */
 gamma = gamma0*xtp;
 det = h*h - (f-gamma)*(g-gamma);
 
 /*   Check on singularity of determinant */
 
 if(fabs(det) <= 1.0e-16) {return 3;}
 
 x08 = (q*h - p*(g-gamma))/det;
 y08 = (p*h - q*(f-gamma))/det;
 
 rsq = x08*x08 + y08*y08 + gamma;
 if (rsq < 0.0) {  return 4;  }		/* a problem, r^2 < 0 */
 *r_out  = (float) sqrt (rsq);
 
 /*  Shift the coordinate system back */
 
 *xc_out = (float) (x08 + sx);
 *yc_out = (float) (y08 + sy);
 
 /*  Calculate pseudo-chi-squared/dof */
 
 if(n == 3) {
    *chi_out = 0.0;
 }
 else {
  chisq8 = 0.0;
  nn = n;	wp = w;		xp = x;		yp = y;
  while (nn--) {
   w8 = *wp++;   x8 = *xp++ - sx;	y8 = *yp++ - sy;
   chisq8 += 
     w8*((x08+sx-x8)*(x08+sx-x8) + (y08+sy-y8)*(y08+sy-y8) - rsq)*
     ((x08+sx-x8)*(x08+sx-x8) + (y08+sy-y8)*(y08+sy-y8) -rsq);
	    }
  *chi_out = (float) (chisq8/(4.0*rsq)/(n-3));
  }
 return iflag;
 }
 /*------------------------------------------------------------------------- */
 void indexf(int n, float ra[],int indx[])
 /* a heap sorter for F*4, returns sorted index */
 {
  int l,j,ir,i,indxt;
  float q;
 
  for (i=0;i<n;i++) indx[i] = i;
  l = (n/2);
  ir = n-1;
  for (;;) {
   if (l > 0)
    q=ra[(indxt=indx[--l])];
   else {
    q=ra[(indxt=indx[ir])];
    indx[ir]=indx[0];
    if (--ir == 0) {indx[0]=indxt;return; }
   }
   i = l;
   j = l + l + 1;
   while (j <= ir) {
    if (j < ir && ra[indx[j]] < ra[indx[j+1]]) j++;
    if (q < ra[indx[j]]) {
     indx[i] = indx[j];
     j += (i=j) + 1;
    }
    else j = ir + 1;
   }
   indx[i] =indxt;
  }
 }
 /*------------------------------------------------------------------------- */
int key_locations(int *array, int n, int *nv, int **counts, int **values, int ***indexptrs)
 /* locate integer keys, nv is the returned # of unique values */
 /* for integer arrays only, finds unique values and all their locations, it is similar
 in some ways to where_unique but also returns all the indices for each unique value
 (key) in a sorted indices array */
 {
 int	iq, i, nd, j, range, type, m, k, mq;
 int	min, max;
 int	nuniq;
 int	*vscratch, *p, *q, val_sym, cnt_sym, indexsymarr;
 int	*countp, **qp, **qpbase, **qpstart, *valueoffsets, *pv;
 /* always need the range */
 mm_int(array, n, &max, &min);
 range = max - min + 1;
 /* need some scratch to keep the temporary values and indices */
 vscratch = (int *) malloc(range * sizeof (int));
 valueoffsets = (int *) malloc(range * sizeof (int));
 bzero( (void *) vscratch, range * sizeof (int));
 bzero( (void *) valueoffsets, range * sizeof (int));
 nuniq = 0;
 /* first a pass to get the count of each unique value and the # of uniques */
 q = array;
 for (j=0;j<n;j++) {
   i = *q++ - min;
   p = vscratch + i;
   if ( *p == 0 ) { nuniq++; }
   *p += 1;
 }
 
 /* and now set up output arrays */
 //printf("unique count = %d\n", nuniq);
 *nv = nuniq;
 /* the scratch array vscratch is a histogram of the unique values, may be sparse */
 /* create arrays with just the unique values and counts */
 p = *values = (int *) malloc(nuniq * sizeof (int));
 q = *counts = (int *) malloc(nuniq * sizeof (int));
 /* the indexptrs term is an array of pointers, each pointing to an array of indices */
 *indexptrs = qpbase = (int **) malloc(nuniq * sizeof (int *));/* there are nuniq ptr to index sets */
 qp = qpstart = (int **) malloc(nuniq * sizeof (int *));  /* for an incrementing set of scratch pointers */
 i = 0;
 pv = valueoffsets;  /* valueoffsets is used later to help populate the offset arrays */
 for (j=0;j<range;j++) {
   int cq = *(vscratch + j);
   if ( cq ) {
     *p++ = j + min;	/* this is the value */
     *q++ = cq;		/* this is the count */
     /* also allocate the arrays for the indices, each is count long */
     *qpbase++ = *qp++ = (int *) malloc(cq * sizeof (int));
     *pv++ = i;  i++;
     if (i > nuniq) { printf("internal error 1 in key_locations\n");  return -1; }
   } else *pv++ = -1;  /* we need a pv for each value within "range" */
 }

 /* now we have to scan the original array again to collect all the locations
 and populate the index arrays, we need valueoffsets to connect things */
 q = array;
 for (j=0;j<n;j++) {
   i = *q++ - min;  /* value in original array */
   /* get the index in the output value array, from a scratch array that covers full range */
   k = *(valueoffsets + i);
   if (k < 0) printf("bad k = %d\n", k);
   if (k >= nuniq) printf("bad k = %d\n", k);
   *(qpstart[k])++ = j;     
 }

 free(vscratch);
 free(qpstart);
 free(valueoffsets);
 return 1;
 }
 /*------------------------------------------------------------------------- */
void mm_f(float *x, int n, float *max, float *min)
 {
 float	rmin, rmax;
 /* just return the min and max of a floating array */
 if (n <= 0 ) { printf("bad count in mm_f, n = %d\n",n); return; }
 rmax = rmin = *x++;	n--;
 while (n--) { if (*x > rmax) rmax = *x; if (*x < rmin) rmin = *x; x++; }
 *min = rmin;	*max = rmax;
 return;
 }
 /*------------------------------------------------------------------------- */
void mm_int(int *x, int n, int *max, int *min)
 {
 int	rmin, rmax;
 /* just return the min and max of a I*4 array */
 if (n <= 0 ) { printf("bad count in mm_int, n = %d\n",n); return; }
 rmax = rmin = *x++;	n--;
 while (n--) { if (*x > rmax) rmax = *x; if (*x < rmin) rmin = *x; x++; }
 *min = rmin;	*max = rmax;
 return;
 }
 /*------------------------------------------------------------------------- */
long getindexofmax(int *x, int n)
  /* x is a 1-D F*4 vector, find the index of the nmax */
 {
 int	rmax;
 long inx, i;
 if (n <= 0 ) { printf("bad count in getindexofmax, n = %d\n",n); return -1; }
 rmax = *x++;	n--;  inx = i = 0;
 while (n--) { i++; if (*x > rmax) { rmax = *x; inx = i; } x++; }
 return inx;
 }
 /*------------------------------------------------------------------------- */
 /* these histogram functions convert float values to rounded integers */
 /*------------------------------------------------------------------------- */
long hist_max_xvalue(float *x, int n)
 /* generates a histogram of the n values in x after casting into rounded integers
 and returns the most populus x value, corrected for min */
 {
 int *hist, xq;
 float xmax, xmin;
 long max, min, imax, i, iq, range;
 mm_f(x, n, &xmax, &xmin);
 max = lrintf(xmax);
 min = lrintf(xmin);
 range = max - min + 1;
 if (range > max4histogram) { printf("values in array > %d\n", max4histogram);  return 0; }
 //printf("max value in hist_max_xvalue %ld, range = %ld\n", max, range);
 hist = (int *) malloc(sizeof(int) * range);
 bzero((void *) hist, sizeof(int) * range);
 for (i=0;i<n;i++) {
   iq = lrintf(x[i]) - min;
   hist[iq] += 1;
 }
 /* now scan for the location of the max */
 imax = 0;
 xq = 0;
 for (i=0;i<range;i++) {
   if (hist[i] > xq) { imax = i;  xq = hist[i]; }
 }
 //printf("max hist value and location %d, %ld\n", xq, imax);
 free(hist);
 imax = imax + min;
 return imax;
 }
 /*------------------------------------------------------------------------- */
int hist_max_freq(float *x, int n)
 /* generates a histogram of the n values in x after casting into rounded integers
 and returns the value of the max entry */
 {
 int *hist, xq;
 float xmax, xmin;
 long max, min, imax, i, iq, range;
 mm_f(x, n, &xmax, &xmin);
 max = lrintf(xmax);
 min = lrintf(xmin);
 range = max - min + 1;
 if (range > max4histogram) { printf("values in array > %d\n", max4histogram);  return 0; }
 //printf("n %d, max value in hist_max_freq %ld, range = %ld\n", n, max, range);
 hist = (int *) malloc(sizeof(int) * range);
 bzero((void *) hist, sizeof(int) * range);
 for (i=0;i<n;i++) {
   iq = lrintf(x[i]) - min;
   hist[iq] += 1;
 }
 /* now scan for the location of the max */
 imax = 0;
 xq = 0;
 for (i=0;i<range;i++) {
   if (hist[i] > xq) { imax = i;  xq = hist[i]; }
 }
 //printf("max hist value and location %d, %ld\n", xq, imax);
 free(hist);
 return xq;
 }
 /*------------------------------------------------------------------------- */
int *hist_array(float *x, int n, int *nh, int *hmin)
 /* generates a histogram of the n values in x after casting into rounded integers
 and returns the histogram array which must be free'ed by calling routine */
 {
 int	*hist, xq;
 float	xmax, xmin;
 long	max, min, i, iq, range;
 mm_f(x, n, &xmax, &xmin);
 max = lrintf(xmax);
 min = lrintf(xmin);
 *nh = range = max - min + 1;
 *hmin = (int) min;
 if (range > max4histogram) { printf("values in array > %d\n", max4histogram);  return 0; }
 //printf("max value in hist_max_xvalue %ld, range = %ld\n", max, range);
 hist = (int *) malloc(sizeof(int) * range);
 bzero((void *) hist, sizeof(int) * range);
 for (i=0;i<n;i++) {
   iq = lrintf(x[i]) - min;
   hist[iq] += 1;
 }
 return hist;
 }
/*--------------------------------------------------------------------------*/
void decomp(double *x, int n, int nd)
/* translated from fortran anadecomp.for which may be easier to follow since
	it uses subscripts rather than the pointers used here */
/* no pivoting in this version !, so diagonals must be != 0 */
{
  register	double	sum1, sum2, *p1, *p2, *p3, *p4;
  double	*qd, * q2, *q1, div, *p;
  int	nq, mq, lq, k;
  p1 = x;
  div = 1.0 / *p1; nq = n-1; while (nq--) { p1 += nd; *p1 *= div; }
  nq = n-1;	k = 1;		qd = x + 1 + nd; /*qd is diagonal ptr */
  while (nq--) {						/*outer loop */
	  /*diagonal term */
	  sum1 = *qd;  mq = k;  p1 = p2 = qd;
	  while (mq--) { p1 -= nd; p2--; sum1 -= (*p1) * (*p2); }
	  *qd = sum1;
	  q1 = q2 = qd; lq = n - k - 1;
	  while (lq--) { q1 += nd;  q2 += 1;  mq = k;
	   sum1 = *q1;	sum2 = *q2;
	   p1 = p4 = qd; p2 = q1; p3 = q2;
	   while (mq--) { p1 -= nd; p2 -= 1; p3 -= nd; p4 -= 1; /*inner */
	   sum1 -= (*p1) * (*p2);  sum2 -= (*p3) * (*p4); }     /*end of */
	  *q2 = sum2;  *q1 = sum1 / *qd;
	  }					/*end of middle loop */
  k++; qd += (nd + 1); }				/*end of outer loop */
}
/*------------------------------------------------------------------------- */
void solve(double *a, double *b, int n, int nd)
{
  register	double	sum, *p1, *p2;
  double	*qd, *q1, div;
  int	nq, mq, lq, k;
  /*printf("in d_solve\n");*/
  p1 = a;	q1 = b; *q1 = *q1 / *p1; q1++;
  nq = n-1;	k = 1;		qd = a + 1 + nd; /*qd is diagonal ptr */
  while (nq--) { sum = 0.0;  mq = k;
	  p1 = qd;	p2 = q1;
	  while (mq--) { p1 -= nd; p2--; sum += (*p1) * (*p2); }
	  *q1 = (*q1 - sum) / *qd;
  k++; qd += (nd + 1); q1++; }	
  nq = n-1;	k = n - 1;	qd = a + nd*nd -2;  q1 = b + n - 2;
  while (nq--) { sum = 0.0;  mq = n - k;
	  p1 = qd;	p2 = q1;
	  while (mq--) { p2++; sum += (*p1) * (*p2); p1 += nd; }
	  *q1 = *q1 - sum;
  k--; qd -= (nd + 1); q1--; }	
}
/*--------------------------------------------------------------------------*/
int pit(float *qxbase, float *qybase, int nxq, int npow, double *a, double *cfbase, double *fbase )
 /* based on ana version, but simplified to handle just F*4 input */
 {
 double	sum, *f, *cf;
 float *qx, *qy;
 int	i, n, ib, ie, k;
 cf = cfbase;	f = fbase;
 qx = qxbase; qy = qybase;

 //for (k=0;k<nxq;k++) printf("k, x(k), y(k) %d %f %f\n",k,qxbase[k],qybase[k]);

 n = nxq;  while (n--) *f++ = *qx++;  /* note that a F*8 copy is made of qxbase */
 *a = nxq;
 for (k=0;k < (2*npow);k++) { sum = 0.0; f = fbase; n = nxq; qx = qxbase;
	 while (n--) {sum += *f;  *f = (*f) * (*qx++); f++; }
	 ib = (k+1-npow) > 0 ? (k+1-npow) : 0;		/*max of */
	 ie = (k+1) > (npow) ? (npow) : (k+1);		/*min of */
	 for (i=ib;i<=ie;i++) {
	 *( a + i + (npow+1)*(ie+ib-i) ) = sum;
	 }
 }
							 /* now the rhs */
 f = fbase; n = nxq;   while (n--) *f++ = *qy++;  /* note that a F*8 copy is made of qybase */
 for (k=0;k <= (npow);k++) { sum = 0.0; f = fbase; n = nxq; qx = qxbase;
	 while (n--) {sum += *f;  *f = (*f) * (*qx++); f++; }
	 *cf++ = sum; }

 
 decomp(a, npow+1, npow+1);
 solve(a, cfbase, npow+1, npow+1);
 
 /*for (k=0;k<=npow;k++) printf("k, cf(k) %d %f\n",k,*(cfbase+k));*/
 return 1;
 }
 /*------------------------------------------------------------------------- */
float parabolicmax5(int *x, int nx)
 /* x is a 1-D eveny spaced vector (I*4 here), do a parabolic fit for 5 points around max and return
 max position (in xq) including fractional value */
 {
 long i, i1, i2;
 float xq, xpit[5], ypit[5];
 float *qybase;
 double *cfbase, *fbase, *a;
 int npow = 2, *ixp;
 /* first first the max position */
 i = getindexofmax(x, nx);
 if (nx < 5) {
    printf("parabolicmax5 input too small for fit, using max index only, no interpolation\n");
    xq = (float) i;
    return xq;
 } else {
    if ((i+2) >= nx) i = nx - 3;
    if (i < 2) i = 2;
    i1 = i - 2;  /* start of subarray to pass to pit */
    ixp = x + i1;
    for (i=0;i<5;i++) xpit[i] = (float) i;	/* the x array for call to cana_pit */
    for (i=0;i<5;i++) ypit[i] = (float) *ixp++;	/* the y array for call to cana_pit */
    cfbase = (double *) malloc((npow+1) * sizeof(double));
    fbase  = (double *) malloc( nx * sizeof(double));
    a = (double *) malloc((npow+1)*(npow+1) * sizeof(double));
    pit( xpit, ypit, 5, 2, a, cfbase, fbase);
    /* the coefficients are in cfbase of length 3 */
    if (cfbase[2] == 0.0) {
      printf("parabolicmax5 - singularity, using max index only, no interpolation\n");
      xq = (float) i;
      return xq;
    }
    xq = - cfbase[1]*.5/cfbase[2];
    xq = xq + i1;
    free(a);  free(fbase);	free(cfbase);
    return xq;
 }
 }
 /*------------------------------------------------------------------------- */
void sobel_index(float *x, int nx, int ny, long *iq, int nindex, float **grad, float **ang)
  /* x is an array (nx, ny), we want to process only the values with indices iq, these
   are offsets in x from 0 and must not be on an edge! nindex is iq count, returned values
   for the gradient and angle in grad and ang, each also nmask values
  */
 {
 float gx, gy, xq;
 float *pbot, *ptop, *pt2, *pt3;
 int	jj;
 long	in;
 pt2 = *grad = malloc(sizeof(float) * nindex);
 pt3 = *ang  = malloc(sizeof(float) * nindex);
 //printf("pt2, pt3 = %p, %p\n", pt2, pt3);
 for (jj=0; jj<nindex; jj++) {
   in = iq[jj];  /* this is index in x */
   //if (jj <20) printf("in = %ld\n", in);

   pbot = x + in - nx -1; /* surrounding points on diagonals */
   ptop = x + in + nx +1;
   //if (jj <20) printf("pbot = %p, ptop = %p\n", pbot, ptop);
   gx = gy = *ptop-- - *pbot++;
   gy = gy + 2 * (*ptop-- - *pbot++);
   xq = *pbot - *ptop;
   gy = gy - xq;
   gx = gx + xq + 2 * (*(pbot + nx) - *(ptop - nx));
   //if (jj <20) printf("grad = %12.4e, ang = %12.4e\n", 0.125 * hypotf(gx, gy),atan2f(gy, gx));
   *pt2++ = 0.125 * hypotf(gx, gy);
   *pt3++ = atan2f(gy, gx);
 }
 return;
 }
 /*------------------------------------------------------------------------- */
// void sobel_mask(float *x, int nx, int ny, long *iq, byte *maskbase, float *grad, float *ang)
//    /* x is an array (nx, ny), we want to process only the values with indices iq, these
//     are offset in x from 0 and must not be on an edge! nindex is iq count, returned values
//     for the gradient and angle in grad and ang, each also nmask values
//    */
//   {
//   float gx, gy, xq;
//   float *pbot, *ptop, *pt2, *pt3;
//   int	jj, in;
//   byte *maskptr, *maskptr2;
//   pt2 = grad = malloc(sizeof(float) * nindex);
//   pt3 = ang  = malloc(sizeof(float) * nindex);
// 
// 
// 
// 
//   for (jj=0; jj<nindex; jj++) {
//     in = iq[jj];  /* this is index in x */
//     
//     pbot = x + in - nx -1; /* surrounding points on diagonals */
//     ptop = x + in + nx +1;
//     gx = gy = *ptop-- - *pbot++;
//     gy = gy + 2 * (*ptop-- - *pbot++);
//     xq = *pbot - *ptop;
//     gy = gy - xq;
//     gx = gx + xq + 2 * (*(pbot + nx) - *(ptop - nx));
//     *pt2++ = 0.125 * hypotf(gx, gy);
//     *pt3++ = atan2f(gy, gx);
//   }
//   return;
//   }
 /*------------------------------------------------------------------------- */
void gsmooth_index(float *x, float fwhm, int nx, int ny, long *iq, int nindex, float *result);
void gsmooth_index(float *x, float fwhm, int nx, int ny, long *iq, int nindex, float *result)
   /* x is an array (nx, ny), we want to smooth only the values with indices iq, these
    are offsets in x from 0 and must not be on an edge! nindex is iq count, returned values
    for the gradient and angle in grad an d ang, each also nmask values
   */
  {
  int	jj, in, i, i1, i2, n2, ng, n, ns;
  float *pt1, *pt2, *pt3, *gkern, width, sum, gsum, sumg, wq, xq;
  width = 0.6005612 * fwhm;	/* convert fwhm to gaussian width*/
  pt2 = result;
  /* check for trivial cases, just return original */
  if ( width < 0.25 ) {
    float *pt = x;
    for (jj=0; jj<nindex; jj++) {  *pt2++ = *pt++; }
    return;
  }
  /* set up the gaussian convolution kernel */
  n2 = MIN( (int) (gsmooth_width *  width), nx/2 - 2);
  n2 = MAX( n2, 1);
  ng = 2 * n2 +1;    /* ng is odd and can't be larger than nx */
  gkern = (float *) malloc( (ng) * sizeof(float));
  n = ng;   wq = (float) (-n2);  pt3 = gkern;  gsum = 0.0;
  while (n--)
    { sum = wq/width; wq++; xq = exp(-sum*sum); *pt3++ = xq; gsum += xq; } 
  for (jj=0; jj<nindex; jj++) {
    in = iq[jj];  /* this is index in image x */
    i = in % nx;  /* the i value */
    i2 = i + n2;
    i1 = i - n2;
    sum = 0;
    /* the two special cases (near edges) should not be common but we have to check */
    if (i1 < 0) {
      /* a left edge case */
      pt3 = gkern - i1;
      ns = ng + i1;
      pt1 = x + in - n2 - i1;
      sumg = sum = 0.0;
      while (ns--) { sumg  += *pt3;  sum += *pt1++ * *pt3++; }
      *pt2++ = sum / sumg;
    } else if (i2 >= nx) {
      /* a right edge case */
      pt3 = gkern;
      pt1 = x + in - n2;
      ns = ng -i2 + nx - 1;  /* ng - (i2 - (nx-1)) */
      sumg = sum = 0.0;
      while (ns--) { sumg  += *pt3;  sum += *pt1++ * *pt3++; }
      *pt2++ = sum / sumg;
    } else {
      /* if no edges, we can use the full gsum to normalize in this case */
      ns = ng;  pt1 = x + in - n2;  pt3 = gkern;
      while (ns--) { sum += *pt1++ * *pt3++; }
      *pt2++ = sum / gsum;
    }
  }
  free(gkern);
  return;
  }
  /*------------------------------------------------------------------------- */
void gsmoothy_index(float *x, float fwhm, int nx, int ny, long *iq, int nindex, float *result);
void gsmoothy_index(float *x, float fwhm, int nx, int ny, long *iq, int nindex, float *result)
   /* the y direction smooth, x is an array (nx, ny), we want to smooth only the values with 
   indices iq, these are offset in x from 0 and must not be on an edge! nindex is iq count, 
   returned values for the gradient and angle in grad an d ang, each also nmask values
   */
  {
  int	jj, in, j, j1, j2, n2, ng, n, ns;
  float *pt1, *pt2, *pt3, *gkern, width, sum, gsum, sumg, wq, xq;
  width = 0.6005612 * fwhm;	/* convert fwhm to gaussian width*/
  pt2 = result;
  /* check for trivial cases, just return original */
  if ( width < 0.25 ) {
    float *pt = x;
    for (jj=0; jj<nindex; jj++) {  *pt2++ = *pt++; }
  }
  /* set up the gaussian convolution kernel */
  n2 = MIN( (int) (gsmooth_width *  width), ny/2 - 2);
  n2 = MAX( n2, 1);
  ng = 2 * n2 +1;    /* ng is odd and can't be larger than ny */
  gkern = (float *) malloc( (ng) * sizeof(float));
  n = ng;   wq = (float) (-n2);  pt3 = gkern;  gsum = 0.0;
  while (n--)
    { sum = wq/width; wq++; xq = exp(-sum*sum); *pt3++ = xq; gsum += xq; } 
  for (jj=0; jj<nindex; jj++) {
    in = iq[jj];  /* this is index in image x */
    j = in / nx;  /* the j value, we are doing y direction here */
    j2 = j + n2;
    j1 = j - n2;
    sum = 0;
    /* the two special cases (near edges) should not be common but we have to check */
    if (j1 < 0) {
      /* a top edge case */
      pt3 = gkern - j1;
      ns = ng + j1;
      pt1 = x + in - j * nx;
      sumg = sum = 0.0;
      while (ns--) { sumg  += *pt3;  sum += *pt1 * *pt3++;  pt1 += nx; }
      *pt2++ = sum / sumg;
    } else if (j2 >= ny) {
      /* a bottom edge case */
      pt3 = gkern;
      pt1 = x + in - n2 * nx;
      ns = ng -j2 + ny - 1;  /* ng - (j2 - (ny-1)) */
      sumg = sum = 0.0;
      while (ns--) { sumg  += *pt3;  sum += *pt1 * *pt3++;   pt1 += nx; }
      *pt2++ = sum / sumg;
    } else {
      /* if no edges, we can use the full gsum to normalize in this case */
      ns = ng;  pt1 = x + in - n2 * nx;  pt3 = gkern;
      while (ns--) { sum += *pt1 * *pt3++;   pt1 += nx; }
      *pt2++ = sum / gsum;
    }
  }
  free(gkern);
  return;
  }
/*------------------------------------------------------------------------- */
void fill_image(float *x, int nx, int ny, long *iq, int nindex, float *fill);
void fill_image(float *x, int nx, int ny, long *iq, int nindex, float *fill)
  /* populate a (nx,ny) image in x with values of fill at indices iq, both iq and fill
  have nindex values. Often x would be zeroed ahead of time but depends on use */
  {
  int	jj;
  float *ptf = fill;
  long *pt;
  pt = iq;
  for (jj=0; jj<nindex; jj++) {
    x[*pt++] = *ptf++;  
  }  
  return;
  }
/*------------------------------------------------------------------------- */
void fill_indices(float *x, int nx, int ny, long *iq, int nindex);
void fill_indices(float *x, int nx, int ny, long *iq, int nindex)
  /* populate a (nx,ny) image in x with values of 1 at indices iq, works like fill_image but
  always fills with 1, mostly for testing */
  {
  int	jj;
  long *pt = iq;
  for (jj=0; jj<nindex; jj++) {
    x[*pt++] = 1.0;  
  }  
  return;
  }
/*------------------------------------------------------------------------- */
int limbcompute(float *x, int nx, int ny, float xcguess, float ycguess, float rguess, float rrange,
  int limbmode, int useprevflag, float fwhm)   /* stand alone limb fitter */
  {
  /* where x is the raw image dimensioned (nx,ny), rguess is approx size of r in pixels,
  rrange is +/- range (in pixels) to use for our annulus */
  /* limbmode is used for difference approaches, 0 for EUV, 1 for 1600/1700, 2 for 4500,
  3 for 304, 4 for HMI */
  /* useprevflag can be set to 1 for the second and subsequent images in a series to use the previous
  fit vales rather than xcguess,,ycguess, rguess */
  /* fwhm is the amoount of smoothing to use for the gradient calculations, 3 is generally good but a
  higher value (like 6) might be necessary for low count, noisy images */
  /* note that result is saved in globals sdisk_xc, sdisk_yc, sdisk_r */
  float r1, r2, r1s, r2s, da;
  //float fwhm = 3.0;
  float *grad, *ang, *fscratch;
  float *result;
  int  *in2;
  int margin, nxny, nxd2, nyd2, nv, num_in2, n2use;
  int *values, *counts, **in4annulus;
  float rmax, limbtolerance;
  float *xlimb, *ylimb, *xlimbptr, *ylimbptr;
  int i, j, k, njs, njs2, njs_s, jprev, jprev2, jprev_s, nq;
  int makeannflag;
  int debug = 1;
  double	t0, t1, t2, t3, t4, t1s, tstart, tend;
  tstart = systime();
  nxny = nx*ny;
  nxd2 = nx/2;  nyd2 = ny/2;
  /* malloc the masks if first call or nx, ny changed */
  if (nx != lastnx || ny != lastny) {
    free(inann);  free(inanns);  free(inann2);
    inann = malloc(sizeof(long) * nxny);
    inanns = malloc(sizeof(long) * nxny);
    inann2 = malloc(sizeof(long) * nxny);
    lastnx = nx;  lastny = ny;
  } 
  /* need a new annulus coordinate set ? */
  makeannflag = 1;
  r1 = rguess - rrange;
  r2 = rguess + rrange;
  if (r1 == lastr1 && r2 == lastr2 ) {
    /* for the center, just have to be close, try 10 */
    if ( (fabsf(lastxc-xcguess) <= 10) && fabsf(lastyc-ycguess) <= 10) makeannflag = 0;
  }
  if (makeannflag) {
    if (verbose > 1) {
      printf("r1, lastr1, r2, lastr2 = %g %g %g %g\n", r1, lastr1, r2, lastr2);
      printf("lastxc, xcguess, lastyc, ycguess = %g %g %g %g\n", lastxc, xcguess, lastyc, ycguess);
    }
  }
  if (makeannflag) {
    float r1uv, r2uv;
    /* use a larger range for the smooth mask because we do a sobel inside afterwards */
    r1s = r1 - 7;   r2s = r2 + 7;
    /* and a special one for UV, used for the double sobel */
    r1uv = r1 - 1;  r2uv = r2 + 1;
    lastr1 = r1;
    lastr2 = r2;
    lastxc = xcguess;
    lastyc = ycguess;
    /* now square the r's for the tests below */
    r1 = r1 * r1;  r2 = r2 * r2;
    r1s = r1s * r1s;  r2s = r2s * r2s;
    r1uv = r1uv * r1uv;  r2uv = r2uv * r2uv;
    /* generate the annulus coordinate arrays and the masks */
    num_inann = num_inann2 = num_inanns = 0;
    margin = 25;
    rmax = 2100.;  /* to exclude corners, only really needed when center is way off as in flats */
    //njs = njs_s = njs2 = 0;
    jprev = jprev_s = jprev2 = -1;
    for (j=margin; j<(ny-margin); j++) {
	float rys = powf ( j - ycguess, 2);
	for (i=margin; i<(nx-margin); i++) {
	  float rxs = powf( i - xcguess, 2);
	  float rs = rys + rxs;
	  /* there are 3 concentric masks, nested */
	  if (rs > r1s && rs < r2s) {
	    /* also check for CCD corners */
	    float rccd = hypotf((float) (i-nxd2), (float) (j-nyd2));
	    if (rccd < rmax) {
	      int iq = i + nx * j;
	      //if (j != jprev_s) { jstarts_s[njs_s] = iq;  jvalues_s[njs_s] = jprev_s = j;  njs_s++; }
	      /* this is the larger mask for the smooth, called masks in ana script */
	      inanns[num_inanns++] = iq;
	      if (rs > r1uv && rs < r2uv) {
		//if (j != jprev2) { jstarts2[njs2] = iq;  jvalues2[njs2] = jprev2 = j;  njs2++; }
		/* this is used for the first UV double sobel, called mask2 in ana script */
		inann2[num_inann2++] = iq;
		if (rs > r1 && rs < r2) {
	          //if (j != jprev) { jstarts[njs] = iq;  jvalues[njs] = jprev = j;  njs++; }
	          /* the sobel mask for EUV, and second sobel for UV, called mask in ana script */ 
		  inann[num_inann++] = iq;
		}
	      }
	  }
	  }
	}
    }
    /* get the angles, only for inann */
    /* these can be carried over so global statics */
    free(ixnann);
    free(iynann);
    free(annth);
    free(annthr);
    free(thann);
    free(in2use);
    ixnann = malloc(sizeof(int) * num_inann);
    iynann = malloc(sizeof(int) * num_inann);
    annth = malloc(sizeof(float) * num_inann);
    annthr = malloc(sizeof(float) * num_inann);
    thann = malloc(sizeof(int) * num_inann);
    in2use = malloc(sizeof(int) * num_inann);
    for (k=0;k<num_inann;k++) {
      float xq;
      long in = inann[k];  /* this is index in image x */
      i = in % nx;  /* the i value */
      j = in / nx;
      ixnann[k] = i;
      iynann[k] = j;
      xq = 5.729578e+01 * atan2f((float) j - ycguess,(float) i - xcguess);
      /* thann is a set for sectoring into 1 degree slots */
      thann[k] = (int) lroundf(xq);
      //if (thann[k] == -180) thann[k] = 180;  /* otherwise we get 361 values */
      if (thann[k] < 0) thann[k] = thann[k] + 360;
      annth[k] = xq;
      /* annthr is a reversed version used for 4500 */
      annthr[k] = annth[k] + 180.;
      if (annthr[k] > 360.) annthr[k] = annthr[k] - 360.;
      //if (k < 20) printf("k, inann[k], thann[k], annth[k], annthr[k] = %d, %ld, %d, %13.5e, %13.5e\n",
      //  k, inann[k], thann[k], annth[k], annthr[k]);
   }
    /* check for sectors with low populations, we are looking for parts of the annulus that are chopped
    off by an edge. Normally each angle zone will have about the same count except when this happens. This
    doesn't remove many indices, just the ones that are partially occulted since the fully occulted ones
    don't appear anyhow. We expect a broad max in population counts (all full sectors will have close to
    the same #), so just eliminate all with < max(counts)*.75 for a trial */
    key_locations(thann, num_inann, &nv, &counts, &values, &in4annulus);
    /* in4annulus are indices in thann */
    /* check out */
//     for (j=0;j<nv;j++) {
//       printf("j, values[j], counts[j] = %d, %d, %d\n", j, values[j], counts[j]);
//     }
    {
      int max, min, nq, *p1, *p2, nc;
      mm_int(counts, nv, &max, &min);
      nq = (int) (0.75 * (float) max);
      p1 = in2use;  /* working space, already allocated above */
      nc = 0;  /* will be final count */
      for (j=0;j<nv;j++) {
	/* check each zones population */
	if (counts[j] >= nq) {
	  /* meets the critera, load this zone */
	  p2 = in4annulus[j];
	  for (i=0;i<counts[j];i++) {
            *p1++ = *p2++;
	    nc++;
	  }
	  free(in4annulus[j]);
	}
      }
      /* in2use are indices in the thann, innan, etc arrays */
      //printf("nc, num_inann = %d, %d\n", nc, num_inann);
      if (nc > num_inann) { printf("internal error, nc, num_inann = %d, %d\n", nc, num_inann); return 1; }
   //    /* in2use now packed with nc indices from inann, reload the other arrays
   //       note that in2use must be monotonic */
   //    p1 = annth;
   //    for (j=0;j<nc;j++) {  *p1++ = annth[in2use[j]]; }

      num_in2use = nc;
    }
    free(counts);
    free(values);

  }  else if (verbose > 1) { printf("re-using last annulus\n"); }
  t1 = systime();
  if (verbose > 0) printf("time for setup %8.3f\n", t1-tstart);



  fscratch = (float *) malloc(sizeof(float) * nxny);  /* don't forget to free */
  bzero((void *) fscratch, nxny);
  /* limbmode is used for different approaches, 0 for EUV, 1 for 1600/1700, 2 for 4500, 3 for 304, 4 for HMI */
  result = malloc(sizeof(float) * num_inanns);
  if (limbmode == 1) {
    /* the UV case uses an extra mask and does a double Sobel */
    gsmooth_index(x, fwhm, nx, ny, inanns, num_inanns, result);
    /* insert this result for the gsmoothy */
    fill_image(fscratch, nx, ny, inanns, num_inanns, result);
    gsmoothy_index(fscratch, fwhm, nx, ny, inanns, num_inanns, result);
    fill_image(fscratch, nx, ny, inanns, num_inanns, result);
    /* and the sobel step, uses mask 2 for first sobel */
    sobel_index(fscratch, nx, ny, inann2, num_inann2, &grad, &ang);
    /* note that presently, grad and ang are allocated in sobel_index and must be freed */
    /* need another sobel here, expand grad into fscratch */
    fill_image(fscratch, nx, ny, inann2, num_inann2, grad);
    free(grad);   free(ang);
    sobel_index(fscratch, nx, ny, inann, num_inann, &grad, &ang);
    /* the final grad/ang are free'ed later */
  } else {
    /* the non-UV case */
    gsmooth_index(x, fwhm, nx, ny, inanns, num_inanns, result);
    /* insert this result for the gsmoothy */
    fill_image(fscratch, nx, ny, inanns, num_inanns, result);
    //fill_image(xs1, nx, ny, inanns, num_inanns, result);  /* for debug */
    gsmoothy_index(fscratch, fwhm, nx, ny, inanns, num_inanns, result);
    fill_image(fscratch, nx, ny, inanns, num_inanns, result);
    //fill_image(ys1, nx, ny, inanns, num_inanns, result);  /* for debug */
    /* and the sobel step, uses a smaller mask */
    sobel_index(fscratch, nx, ny, inann, num_inann, &grad, &ang);
    //fill_image(gradout, nx, ny, inann, num_inann, grad);
    //fill_image(angout, nx, ny, inann, num_inann, ang);
    /* note that presently, grad and ang are allocated in sobel_index and must be freed */
  }
 free(fscratch);
 free(result);
 t2 = systime();
 if (verbose > 0) printf("time for image processing %8.3f\n", t2-t1);
 /* only inann and in2use matter now, done with inanns and inann2 until next image */
 //printf("theta (aka ang here)\n");
//  for (k=0;k<20;k++) printf("k, inann[k], ang[k], grad[k] = %d, %ld %13.5e, %13.5e\n", k, inann[k], (float) 57.29578 * ang[k], grad[k]);
//  for (k=0;k<20;k++) {
//    int iq = in2use[k];
//    printf("k, iq, inann[iq], ang[iq], grad[iq] = %d, %d, %ld %13.5e, %13.5e\n", k, iq, inann[iq], (float) 57.29578 * ang[iq], grad[iq]);
//  }

 /* mark all within 5 degrees of expected orientation */
 if (limbmode == 3) da = 5.0; else da = 5.0;   /* originally tried different da's for different modes */
 /* ang is in radians and represents the inann set of indices from original image, annth and annthr are
 in degrees */
 {
   float *fp1, *fp2, xq, gq;
   int *p1, *p2, iq, jq, nc, nmatch = 0;
   if (limbmode == 2 || limbmode == 3 || limbmode == 4) {
      fp1 = annthr;   /* use reverse angles */
   } else {
      fp1 = annth;    /* use regular angles */
   }
   in2 = p1 = malloc(sizeof(int) * num_in2use);
   nc = 0;
   if (verbose > 2) printf("num_in2use %d\n", num_in2use);
   for (k=0;k<num_in2use;k++) {
     iq = in2use[k];  /* index for ang and annth (or annthr) */
     xq = fabsf( fp1[iq] - (float) 57.29578 * ang[iq]);
     if ((xq <= da || xq >= (360-da)) && (grad[iq] > 0.0)) { *p1++ = iq; nc++;  /* loading in2 */ 
//        if (nmatch < 20) {
//          printf("a match at k = %d, ang[iq], grad[iq] = %13.5e, %13.5e\n", k, ang[iq], grad[iq]);
// 	 nmatch++;
//        }
     }
   }
   num_in2 = nc;
   if (verbose > 2)
     printf("n within 5 degrees, num_inann = %d, %d\n", num_in2, num_inann);
 }
 /* note that in2 are still indices in the full num_inann set */
 {
   int *theta3, *p1, *p2, *indx, *infaint, *svalues;
   int keepfaint = 300, nfaint, nzones, iq, jq, nc;
   float *xint, *fp1, xtot, *sxint;
   //printf("num_in2 = %d\n", num_in2);
   theta3 = p1 = malloc(sizeof(int) * num_in2);  /* checked free */
   for (k=0;k<num_in2;k++) {
     *p1++ = thann[in2[k]];  /* loading theta3 */
     //if (k <20) printf("k, in2[k], inann[in2[k]], theta3 = %d %d %ld %d\n", k, in2[k], inann[in2[k]], theta3[k]);
   }
   /* run key_locations on theta3, note re-use of other arguments */
   //printf("size of theta3 %d\n", num_in2);
   key_locations(theta3, num_in2, &nzones, &counts, &values, &in4annulus);
   /* the nzones pointers are in in4annulus, these are indices in theta3 */
   if (verbose > 2) printf("nzones = %d\n", nzones);
   indx = malloc(sizeof(int) * nzones);  /* checked free */
   /* next we want to eliminate bright zones on the limb if EUV (limbmode 0 or 3) */
   /* get mean intensity in each zone */
   if (limbmode == 0 || limbmode == 3) {
     xint = fp1 = malloc(sizeof(float) * nzones);  /* checked free */
     for (j=0;j<nzones;j++) {
       xtot = 0.0;
       nc = counts[j];
       p2 = in4annulus[j];
       for (i=0;i<nc;i++) {
	 iq = *p2++;  	/* iq is an index in theta3 */
	 jq = in2[iq];	/* jq is index in thann and inann */
	 xtot += (float) x[inann[jq]];
       }
       *fp1++ = xtot/(float) nc;
       //printf("j, nc, value = %d %d %d\n", j, nc, values[j]);
     }
   /* now keep either all or part of the sorted indx, these are indices in the 360 zones */
     /* now sort the mean intensities in xint */
     indexf(nzones, xint, indx);  /* sort using the mean intensities */
     /* keep a fraction equal to keepfaint/360 for limbmode 0 or 3 */
     nfaint = (int) ((float) keepfaint * (float) nzones/360.);
     free(xint);   /* cleanup, just need the indices now */
   } else {
     /* all are used for limbmode 1, 2, 4 */
     nfaint = nzones;
     for (j=0;j<nzones;j++) { indx[j] = j; }
   }
   //printf("keeping %d\n", nfaint);

   /* note we malloc'ed nzones but will only use nfaint of them */
   /* set up an index array for the ones we want, need a new array called in2use from the ana code */


   if (useprevflag) {
     /* set to last values */
     xcenters[0] = sdisk_xc; ycenters[0] = sdisk_yc; radiai[0] = sdisk_r; pcounts[0] = 0;
     nextcenter = 1;
   } else {
   /* this rather large section is just used to get a stable estimate of the limb from
   a guess that may be fairly far off though still in the annulus, it is normally used just
   for the first image in a time series, after that we normally use the previous fit */
   /* loop over initial angular zones and find a point in each, we are
   only interested in the in3 coords */
   float *ixlimb, *iylimb, *rsec, *weights, rq, thq, *rdlimb;
   long irmax, *inx4ann;
   int jfaint, *indxrdlimb;
   ixlimb = malloc(sizeof(float) * nfaint);  /* checked free */
   iylimb = malloc(sizeof(float) * nfaint);  /* checked free */
   jfaint = 0;  /* use jfaint to load ixlimb in case there are zones with nc =0 */
   for (j=0;j<nfaint;j++) {
     k = indx[j];
      //printf("j, k, counts = %d, %d, %d\n", j, k, counts[k]);
     nc = counts[k];
     p2 = in4annulus[k];  /* indices in theta3 */
     if (nc >0) {
       if (limbmode == 1 || limbmode == 3) {  /* this is the uv and 304 cases */
	 rsec = malloc(sizeof(float) * nc);  /* checked free */
	 for (i=0;i<nc;i++) {
	   float ix, iy;
	   iq = *p2++;  	/* iq is an index in theta3 */
	   jq = in2[iq];	/* jq is index in thann and inann */
	   ix = ixnann[jq];
	   iy = iynann[jq];
	   rsec[i] = 0.5 * hypotf( ix - xcguess, iy - ycguess);
	   if (limbmode == 3) rsec[i] = rsec[i] * .25; /* this squeezes the distribution for 304 */
	 }
	 /* now get the r with the most entries */
	 irmax = hist_max_xvalue(rsec, nc);
         /* irmax is the integer rounded most common radius */
	 rq = 2.0 * (float) irmax;  /* have to use for ixlimb and iylimb */
	 if (limbmode == 3) rq = rq *4.0;
	 /* we still have k, use to get nominal angle for this zone, this is in values from the the last
	 key_locations */
	 thq = values[k] * 0.01745329;
	 ixlimb[jfaint] = rq *cosf(thq) + xcguess;
	 iylimb[jfaint] = rq *sinf(thq) + ycguess;
	 free(rsec);
       } else {
	 /* the limbmode = 0, 2, and 4 cases, find the max grad in each zone */
	 float gmax;
	 int jqgmax;
	 gmax = -1.0E20;  jqgmax = 0;
	 for (i=0;i<nc;i++) {
	   iq = *p2++;  	/* iq is an index in theta3 */
	   jq = in2[iq];	/* jq is index in thann and inann and grad */
           if ( grad[jq] > gmax) { gmax = grad[jq]; jqgmax = jq; }
	 }
	 ixlimb[jfaint] = ixnann[jqgmax];
	 iylimb[jfaint] = iynann[jqgmax];
	 //printf("j, jfaint, jqgmax, ixlimb[j], iylimb[j], gmax = %d, %d, %d, %g, %g, %g\n", j, jfaint, jqgmax, ixlimb[j], iylimb[j], gmax);
       }
       jfaint++;
     }
   }
   /* for all the cases we now have nfaint entries in ixlimb and iylimb */
   if (nfaint != jfaint) printf("nfaint, jfaint do not match: %d %d\n", nfaint, jfaint);
   nfaint = jfaint;  /* usually the same, maybe always? */
   /* now we use the ixlimb, iylimb pairs */
   weights = malloc(sizeof(float) * nfaint);  /* weights, just all 1.0 here */ /* checked free */
   for (j=0;j<nfaint;j++) { weights[j] = 1.0; }
   circle_fit_lsq(ixlimb, iylimb, weights, nfaint, &sdisk_xc, &sdisk_yc, &sdisk_r, &chicircle);
   if (verbose > 0)
     printf("sdisk_xc, sdisk_yc, sdisk_r = %g, %g, %g\n",sdisk_xc, sdisk_yc, sdisk_r);
   
   xcenters[0] = sdisk_xc; ycenters[0] = sdisk_yc; radiai[0] = sdisk_r; pcounts[0] = nfaint;
   nextcenter = 1;

   /* the ones inside the circle are the ones most likely to be bad, so eliminate them first */
   /* we compute them all and sort them because we only want to eliminate the worse 25% or all insiders
      more than 7 pixels from circle, whichever is smaller */
   rdlimb = malloc(sizeof(float) * nfaint); /* checked free */
   for (j=0;j<nfaint;j++) {
     rdlimb[j] = hypotf(ixlimb[j]-sdisk_xc, iylimb[j]-sdisk_yc) - sdisk_r;
   }
   indxrdlimb = malloc(sizeof(int) * nfaint); /* checked free */
   indexf(nfaint, rdlimb, indxrdlimb);
   /* allocate for the reduced set, use nfaint so oversized */
   {
     float *ixlimb2, *iylimb2;
     int n2reject, nreject, n2use;
     n2reject = nfaint/4;
     ixlimb2 = malloc(sizeof(float) * nfaint);
     iylimb2 = malloc(sizeof(float) * nfaint);
     n2use = 0;  
     nreject = 0;
     for (j=0;j<nfaint;j++) {
       i = indxrdlimb[j];
       if ( (rdlimb[i] > -7.0 ) || (nreject > n2reject) ) {
         /* keep this one */
	 ixlimb2[n2use] = ixlimb[i];    iylimb2[n2use] = iylimb[i];   n2use++;
       } else {
         /* a reject */
	 nreject++;  /* when nreject exceeds n2reject, we will keep all */
       }
       //printf("j, n2use, rdlimb[i] = %d, %d, %g\n",j, n2use, rdlimb[i]);
     }
     /* n2use is the new count for a circle fit, it is smaller than nfaint so reuse the weights */
     if (verbose > 1) printf("# of values sent to circle_fit %d\n", n2use);
     circle_fit_lsq(ixlimb2, iylimb2, weights, n2use, &sdisk_xc, &sdisk_yc, &sdisk_r, &chicircle);
     xcenters[nextcenter] = sdisk_xc; ycenters[nextcenter] = sdisk_yc;
     radiai[nextcenter] = sdisk_r; pcounts[nextcenter] = n2use;  nextcenter++;
     if (verbose > 0)
       printf("sdisk_xc, sdisk_yc, sdisk_r = %g, %g, %g\n",sdisk_xc, sdisk_yc, sdisk_r);
     free(ixlimb2);
     free(iylimb2);
   }
   free(weights);  /* OK to free this now, we used it twice */
   free(ixlimb);
   free(iylimb);
   free(indxrdlimb);
   free(rdlimb);
  
 }
 /* we have either the previous center or a preliminary computation, now zero in */

 if (limbmode != 0 && limbmode != 3) {  /* not EUV or 304, gets 1, 2, and 4 */
   float *g3, *g3ptr;
   /* these cases will use more circle fitting but we want the larger set of points so collect */
   if (verbose) printf("in UV or 304 zone, nfaint = %d\n", nfaint);
   nc = 0;  /* get total needed */
   for (j=0;j<nfaint;j++) {
     k = indx[j];
     nc += counts[k];
   }
   n2use = nc;
   xlimb = xlimbptr = malloc(sizeof(float) * n2use); /* checked free */
   ylimb = ylimbptr = malloc(sizeof(float) * n2use); /* checked free */
   /* need a set of gradients that match the xlimb, ylimb */
   g3 = g3ptr = malloc(sizeof(float) * n2use); /* checked free */
   /* load them */
   for (j=0;j<nfaint;j++) {
     k = indx[j];
     nc = counts[k];
     p2 = in4annulus[k];
     if (nc >0) {
	 for (i=0;i<nc;i++) {
	   //printf("i = %d\n", i);
	   iq = p2[i];  	/* iq is an index in theta3 */
	   //printf("iq = %d\n", iq);
	   jq = in2[iq];	/* jq is index in thann and inann and grad */
	   //printf("jq = %d\n", jq);
	   *xlimbptr++ = ixnann[jq];
	   *ylimbptr++ = iynann[jq];
	   *g3ptr++ = grad[jq];
	 }
     }
   }
   /* first a coarse zeroinonlimb for all 3 of these cases */
   {
     float drs[] = {50., 35., 25.,15, 10., 5};
     int ndr = sizeof(drs)/sizeof(float);
     // printf("ndr = %d,\n",ndr);
     for (k=0;k<ndr;k++) {
       zeroinonlimb(xlimb, ylimb, n2use, drs[k], g3);
       xcenters[nextcenter] = sdisk_xc; ycenters[nextcenter] = sdisk_yc;
       radiai[nextcenter] = sdisk_r; pcounts[nextcenter] = nc;  nextcenter++;
     }
   }

   /* for cases 2 or 4, we will finish up with just some more of the zeroinonlimb fits */
   if (limbmode == 2 || limbmode == 4) {
     float drs[] = {3., 2., 1., 1., 1., 1., 1.};
     int ndr = sizeof(drs)/sizeof(float);
     // printf("ndr = %d,\n",ndr);
     for (k=0;k<ndr;k++) {
       zeroinonlimb(xlimb, ylimb, n2use, drs[k], g3);
       xcenters[nextcenter] = sdisk_xc; ycenters[nextcenter] = sdisk_yc;
       radiai[nextcenter] = sdisk_r; pcounts[nextcenter] = nc;  nextcenter++;
     }
   }

   free(xlimb);  free(ylimb);  free(g3);
 }
 
 /* the Hough fits are done for cases 0, 1, and 3 (we are done fitting for 2 and 4) */
 if (limbmode != 2 && limbmode != 4) {
   /*  collect the set of coordinates that match a dr condition and a local threshold in g  */
   /* assume that we are within 10 pixels to reduce candidate points but make it 15 for 304 */
   if (limbmode == 3) limbtolerance = 15; else limbtolerance = 10;
   /* we will replace the set of in4annulus pointer arrays with reduced sets that meet our conditions */

   n2use = 0;  /* get total */
   for (j=0;j<nfaint;j++) {  /* loops through the 1 deg angular zones, but not always all */
     int ncnew;
     k = indx[j];
     nc = counts[k];
     ncnew = 0;
     p2 = in4annulus[k];  /* indices in theta3 */
     if (nc >0) {
	 for (i=0;i<nc;i++) {
	   float ix, iy, rq, drq;
	   iq = p2[i];  	/* iq is an index in theta3 */
	   jq = in2[iq];	/* jq is index in thann and inann and grad */
	   ix = ixnann[jq];
	   iy = iynann[jq];

	   rq = hypotf( ix - sdisk_xc, iy - sdisk_yc);
	   drq = fabsf(rq - sdisk_r);
	   if (drq < limbtolerance) {
	     /* use this one, it meets the conditions */
	     p2[ncnew] = iq;   ncnew++;
	   }
	 }
     }
     counts[k] = ncnew;  /* note that this could be 0 */
     n2use += ncnew;
   }
   if (verbose > 2) printf("new total after limbtolerance %d\n", n2use);

   /* now have a new set in counts and in4annulus that are within limbtolerance, now work on g similarly */

   n2use = 0;  /* get total */
   for (j=0;j<nfaint;j++) {  /* loops through the 1 deg angular zones, but not always all */
     float gmean;
     int ncnew;
     k = indx[j];  /* the intensity sorting */
     nc = counts[k];
     p2 = in4annulus[k];
     gmean = 0.0;
     ncnew = 0;
     if (nc >0) {
	 for (i=0;i<nc;i++) {
	   iq = p2[i];  	/* iq is an index in theta3 */
	   jq = in2[iq];	/* jq is index in thann and inann and grad */
	   /* within this sector, compute the mean g first */
	   gmean = gmean + grad[jq];
	 }
	 gmean = gmean/nc;
	 /* cycle through again to select using gmean */
	 gmean = gmean * 0.9;	/* our threshold */
	 for (i=0;i<nc;i++) {
	   iq = p2[i];  	/* iq is an index in theta3 */
	   jq = in2[iq];	/* jq is index in thann and inann and grad */
	   if (grad[jq] > gmean) { p2[ncnew] = iq;   ncnew++; }
	 }
     }
     //printf("j, values[k], nc, ncnew %d %d %d %d\n", j, values[k], nc, ncnew);
     counts[k] = ncnew;  /* note that this could be 0 */
     n2use += ncnew;
   }
   if (verbose > 2) printf("new total after g check %d\n", n2use);

   /* we need arrays of the remaining ix and iy values */
   if (verbose > 2) printf("next new total (n2use) = %d,\n",n2use);
   xlimb = xlimbptr = malloc(sizeof(float) * n2use);
   ylimb = ylimbptr = malloc(sizeof(float) * n2use);
   /* load them */
   for (j=0;j<nfaint;j++) {
     k = indx[j];  /* the intensity sorting */
     nc = counts[k];
     p2 = in4annulus[k];
     if (nc >0) {
	 for (i=0;i<nc;i++) {
	   iq = p2[i];  	/* iq is an index in theta3 */
	   jq = in2[iq];	/* jq is index in thann and inann and grad */
	   *xlimbptr++ = ixnann[jq];
	   *ylimbptr++ = iynann[jq];
	 }
     }
     free(p2);
   }
   free(in2);  free(counts);   free(values);

   /* so we have nc coordinate pairs, do the Hough iterations */
   t3 = systime();
   if (verbose > 0) printf("time for pre fit %8.3f\n", t3-t2);
   minihough(n2use, xlimb, ylimb, 30);
   free(xlimb);  free(ylimb);
   /* result */
   if (verbose > 0)
     printf("final sdisk_xc, sdisk_yc, sdisk_r = %g, %g, %g\n",sdisk_xc, sdisk_yc, sdisk_r);
 }
 t4 = systime();
 if (verbose > 0) printf("time for Hough %8.3f\n", t4-t3);
 free(theta3);
 free(indx);
 }
 free(grad);  free(ang);
 tend = systime();
 if (verbose > 0) printf("total time %8.3f\n", tend-tstart);
 return 0;  /* success return */   
 }
