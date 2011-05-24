/* code pulled from ana for a standalone image rotate/magnify */
#define MAXREGRIDSIZE 10000
#include "defs.h"
#include <stdio.h>
#include <stdlib.h>
#include "ana_structures.h"
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/time.h>
#define ABS(x) ((x)>=0?(x):-(x))
#define	MIN(a,b) (((a)<(b))?(a):(b))
#define	MAX(a,b) (((a)>(b))?(a):(b))
#define BI_CUBIC_SMOOTH 4
#define BI_CUBIC 3
 static int  ana_type_size[5]={1,2,4,4,8};
 union	types_ptr { byte *b; short *w; int *l; float *f; double *d;};		
 static	int	nm1, mm1, nm2, mm2, n, data_type, stretchmark_flag;
 static	float	fnm1, fnm5, fmm1, fmm5, xl, yl;
 static	union	types_ptr base, out;
 static	int	regridtypeflag = 0;
 static int regrid(void *, int , int ,float *, float *, int , int ,int , int ,void **, int *, int *);
 static void bicubic_f(), bicubic_fc();
 static int	resample_type = BI_CUBIC_SMOOTH;
 int image_magrotate(void *,int,int,int,float,float,float,float,void **,int *,int *,int,int);
 /*------------------------------------------------------------------------- */
int image_magrotate(
     /* rotate about center, magnify (rescale) and then offset the input image
     this is usually the preferred order for aligning to another image, order
     can be varied by constructing the final grid accordingly */
     /* theta is assumed in degrees */
     void *array, int nin, int min, int data_type_input,
     float theta, float mag, float dx, float dy,
     void **outarray, int *nx, int *ny,
     int regridtype_input, int stretchmark_flag_input
     )
 {
 float cx[4], cy[4], xc, yc, cq, sq, cs, magr;
 int i, stat, m;
 /* set some local "globals" that are used by several routines */
 regridtypeflag = regridtype_input;
 stretchmark_flag = stretchmark_flag_input;
 data_type = data_type_input;
 n = nin;  m = min;
 /* setup the x and y grid */
 cx[0] = 0.0; cx[1] = (float) n; cx[2] = 0.0; cx[3] = (float) n;
 cy[0] = 0.0; cy[2] = (float) m; cy[1] = 0.0; cy[3] = (float) m;
 xc = 0.5 * (float) (n-1);
 yc = 0.5 * (float) (m-1);
 magr = 1./mag;
 /* center before rotating grid and apply rescale (note - this assumes we
 want the result array to have the same dimensions as the original, appropiate
 only for small rescales) */
 for (i=0;i<4;i++) { cx[i] = (cx[i] - xc)*magr; cy[i] = (cy[i] - yc)*magr; }
 theta = theta * 0.01745329252;   /* remove if theta already in radians */
 cq = cos(theta);   sq = sin(theta);
 /* rotate and uncenter grid */
 for (i=0;i<4;i++) {
   cs = cx[i]; cx[i] =    cs*cq + cy[i]*sq + xc;
               cy[i] = cy[i]*cq -   cs *sq + yc;
 }
 /* and a final offset */
 for (i=0;i<4;i++) { cx[i] = cx[i] - dx;   cy[i] = cy[i] - dy; }
 /* note that nx and ny should be same as n and m for this case */

 stat = regrid(array, n, m, cx, cy, 2, 2, n, m, outarray, nx, ny);
 
 if (stat) { printf("error in regrid\n"); *outarray = 0; nx = ny = 0; return -1; }
 return 0;
 }
 /*------------------------------------------------------------------------- */
int regrid(
     void *array, int n, int m,
     float *xgbase, float *ygbase, int ng, int mg,
     int ns, int ms,
     void **outarray, int *nxout, int *nyout
     )
 /* also uses several globals shared with bicubic_f and bicubic_fc */
 {
 int	iq, nx, ny, ngrun;
 int	iprun, jrun, jprun, ig, ic, jc, jq, result_sym;
 int	i, j, ind;
 float	fn, fm, yrun, ax, bx, cx, dx, ay, by, cy, dy, xq, beta, xinc, yinc;
 float	xl0, yl0;
 struct	ahead	*h;
 union	types_ptr jpbase, jbase, ipbase, bb;
 void *bq;
 base.l = (int *) array;
 fn = (float) n;	fm = (float) m;
 ngrun = ng;
 /* generate the output array */
 ng--;	mg--;
 nx = ng * ns;	ny = mg * ms;
 *nxout = nx;  *nyout = ny;
 if ( nx > MAXREGRIDSIZE || ny > MAXREGRIDSIZE ) {
	 printf("result array in REGRID would be %d by %d\n",nx,ny);
	 printf("which exceeds current MAXREGRIDSIZE (%d)\n",MAXREGRIDSIZE);
	 return -1; }
 /* creat the result array, the same type as input (data_type) */
 i = ana_type_size[data_type];
 bq = malloc(nx * ny * i);
 
 *outarray = bq;
 jpbase.l = (int *) *outarray;
 yrun = 1.0/ (float) ms;
 /* various increments, in bytes! */
 iprun = ns * i;	jrun = ng * iprun;	jprun = ms * jrun;
 ind = 0;
 /* before any looping, branch on the type of regrid we are doing */
 switch (regridtypeflag) {
 case 0:	/* nearest neighbor regrid */
 /* start 4 level loop */
 while (mg--) 	{					/* outer mg loop */
   ipbase.b = jpbase.b;	jpbase.b = ipbase.b + jprun;
   ig = ng;	j = ind;
   /* ig loop */
   while (ig--) {
				   /* get the 8 grid values for this cell */
     i = j;
     ax = xgbase[i];		ay = ygbase[i];		i++;
     bx = xgbase[i] - ax;	by = ygbase[i] - ay;	i += ngrun;
     dx = xgbase[i] - ax;	dy = ygbase[i] - ay;	i--;
     cx = xgbase[i] - ax;	cy = ygbase[i] - ay;
     dx = dx - bx - cx;	dy = dy - by - cy;	j++;
     xq = 1.0 /(float) ns;
     bx *= xq;	by *= xq;	dx *= xq * yrun;	dy *= xq * yrun;
     cx *= yrun;	cy *= yrun;
     /* setup for jc loop */
     jbase.b = ipbase.b;	ipbase.b = ipbase.b + iprun;
     beta = 0.0;	jc = ms;
     while (jc--) {
       /* setup for ic loop */
       out.b = jbase.b;	jbase.b = jbase.b + jrun;
       xl = ax + beta * cx + 0.5;	yl = ay + beta * cy + 0.5;
       xinc = (bx + beta * dx);	yinc = (by + beta * dy);
       ic = ns;
       /* type switch for inner loop */
       switch (data_type) {
	 case 0:
	   while (ic--) {
	   if ( xl < 0 || xl >= fn || yl < 0 || yl >= fm) *out.b++ = 0;
	   /* else { iq = xl;	jq = yl;  iq += n*jq; *out.b++ = *(base.b + iq); }*/
	   else { *out.b++ = *(base.b + (int) xl + n * (int) yl); }
	   xl += xinc;	yl += yinc;	} break;
	 case 1:
	   while (ic--) {
	   if ( xl < 0 || xl >= fn || yl < 0 || yl >= fm) *out.w++ = 0;
	   else { *out.w++ = *(base.w + (int) xl + n * (int) yl); }
	   xl += xinc;	yl += yinc;	} break;
	 case 2:
	   while (ic--) {
	   if ( xl < 0 || xl >= fn || yl < 0 || yl >= fm) *out.l++ = 0;
	   else { *out.l++ = *(base.l + (int) xl + n * (int) yl); }
	   xl += xinc;	yl += yinc;	} break;
	 case 3:
	   while (ic--) {
	   if ( xl < 0 || xl >= fn || yl < 0 || yl >= fm) *out.f++ = 0;
	   else { *out.f++ = *(base.f + (int) xl + n * (int) yl); }
	   xl += xinc;	yl += yinc;	} break;
	 case 4:
	   while (ic--) {
	   if ( xl < 0 || xl >= fn || yl < 0 || yl >= fm) *out.d++ = 0;
	   else { *out.d++ = *(base.d + (int) xl + n * (int) yl); }
	   xl += xinc;	yl += yinc;	} break;
       }   /* end of switch on type for inner loop */
	 beta++;
     }
   }
     ind += ngrun;
 }
 break;		/* end of nearest neighbor regrid case */
 
 case 1:
  {	/* bicubic with or without stretchmarks case */
   mm1 = m-1;  nm1 = n-1;  nm2 = n-2; mm2 = m-2;
   fnm1 = fn - 1.0;	fnm5 = fn -0.5;
   fmm1 = fm - 1.0;	fmm5 = fm -0.5;
			  /* start 4 level loop */
  do 	{						/* outer mg loop */
    ipbase.b = jpbase.b;	jpbase.b = ipbase.b + jprun;
    ig = ng;	j = ind;
							   /* ig loop */
    do	{
      /* get the 8 grid values for this cell */
      i = j;
      ax = xgbase[i];		ay = ygbase[i];		i++;
      bx = xgbase[i] - ax;	by = ygbase[i] - ay;	i += ngrun;
      dx = xgbase[i] - ax;	dy = ygbase[i] - ay;	i--;
      cx = xgbase[i] - ax;	cy = ygbase[i] - ay;
      dx = dx - bx - cx;	dy = dy - by - cy;	j++;
      xq = 1.0 /(float) ns;
      bx *= xq;	by *= xq;	dx *= xq * yrun;	dy *= xq * yrun;
      cx *= yrun;	cy *= yrun;
						      /* setup for jc loop */
      jbase.b = ipbase.b;	ipbase.b = ipbase.b + iprun;
      beta = 0.0;	jc = ms;
      /* LS notes that adding increments to xl and yl for the inner loop below
      results in roundoff errors for F*4 computations, see notes further below */
      while (jc--) {
	/* setup for ic loop */
	/* some optimizer trouble here on Alpha, re-arrange */
	yinc = (by + beta * dy);
	//ic = ns;
	out.b = jbase.b;	jbase.b = jbase.b + jrun;
	/* xl0 and yl0 added as bases for position computation */
	xl = xl0 = ax+beta*cx;  yl = yl0 = ay+beta*cy;
	xinc = (bx + beta * dx);
	for (ic=0; ic < ns; ic++) {
 	       xl = xl0 + ic*xinc;	yl = yl0 + ic*yinc;
 	       switch (resample_type) {
		 case BI_CUBIC_SMOOTH:	bicubic_f();  break;
		 case BI_CUBIC:		bicubic_fc(); break;
	       }
	}
	beta++;
      }
    } while (--ig > 0);
    ind += ngrun;
  } while (--mg > 0);
 } break;		/* end of bicubic with stretchmarks case */
 
 default:  { printf("illegal regridtypeflag = %d in regrid\n", regridtypeflag);	return -1; }
 }
 
 return 0;
 }
 /*------------------------------------------------------------------------- */
void bicubic_f()	/* internal routine for single pixel */
 {
 /* used by all (most?) routines that do bi-cubic interpolations, runs
 a little slower than the originals, perhaps because of the overhead
 in the call or other adjustments made */
 /*    MODIFIED 1-26-84 TO USE LOWER NOISE CUBIC INTERPOLATION FOUND IN
      S.K. PARK & R.A. SCHOWENGERDT, COMP. VIS. & IM. PROC., VOL. 23,
      P. 258, 1983:  USES THEIR FORMULA WITH ALPHA = -0.5
@Article{park83image,
  author       = "S. K. Park and R. A. Schowengerdt",
  title	       = "Image Reconstruction by Parametric Cubic
                  Convolution",
  journal      = "Computer Vision, Graphics, and Image Processing",
  year	       = 1983,
  volume       = 23,
  pages	       = "258-272"
  }
*/
 int	i1, i2, i3, i4, j1, j2, j3, j4, iq;
 float	c1, c2, c3, c4, b1, b2, b3, b4, dx0, dx1, dx2, dx3, dx4, xq, yq;
 union	types_ptr bb;
 /* the location is in xl, yl; base is the pointer to array; out is
 pointer to output; both are unions */
 i2 = (int) xl;		j2 = (int) yl;
 if ( i2 >= 1 && i2 < nm2 ) {		/* normal interior */
	 dx0 = xl - i2; i1 = i2 - 1; i2 = 1; i3 = 2; i4 = 3;	 }
	 else {				/* edge cases */
	 /* check if stretchmarks required */
	 if (stretchmark_flag == 0) {
	  if ( xl < -0.5 || xl > fnm5) {
	   switch (data_type) {
	   case 0: *out.b++ = 0; break;
	   case 1: *out.w++ = 0; break;
	   case 2: *out.l++ = 0; break;
	   case 3: *out.f++ = 0; break;
	   case 4: *out.d++ = 0; break;
	   }
	  return;
	  }
	 }
	 i2 = MIN(i2, nm1);	i2 = MAX( i2, 0);
	 xq = MIN(xl, fnm1);	xq = MAX(xq, 0);
	 dx0 = xq - i2;
	 i1 = MIN(i2-1, nm1);	i1 = MAX( i1, 0);
	 i3 = MIN(i2+1, nm1);	/* i3 = MAX( i3, 0); */
	 i4 = MIN(i2+2, nm1);	/* i4 = MAX( i4, 0); */
	 i2 = i2 - i1;	i3 = i3 - i1;	i4 = i4 - i1;
	 }
 dx1 = 1.0 - dx0;  dx2 = -dx0 * 0.5;  dx3 = dx0 * dx2;
 dx4 = 3. * dx0 * dx3;
 c1 = dx2*dx1*dx1; c2 = 1.-dx4+5.*dx3; c3 = dx4-(dx2+4.*dx3); c4 = dx3*dx1;
 if ( j2 >= 1 && j2 < mm2 ) {		/* normal interior */
	 j1 = j2 - 1; dx0 = yl - j2; j3 = n; j2 = n; j4 = n;	}
	 else {				/* edge cases */
	 /* check if stretchmarks required */
	 if (stretchmark_flag == 0) {
	  if ( yl < -0.5 || yl > fmm5) {
	   switch (data_type) {
	   case 0: *out.b++ = 0; break;
	   case 1: *out.w++ = 0; break;
	   case 2: *out.l++ = 0; break;
	   case 3: *out.f++ = 0; break;
	   case 4: *out.d++ = 0; break;
	   }
	  return;
	  }
	 }	 j2 = MIN(j2, mm1);	j2 = MAX( j2, 0);
	 xq = MIN(yl, fmm1);	xq = MAX(xq, 0);
	 dx0 = xq - j2;
	 j1 = MIN(j2-1, mm1);	j1 = MAX( j1, 0);
	 j3 = MIN(j2+1, mm1);	/* j3 = MAX( j3, 0); */
	 j4 = MIN(j2+2, mm1);	/* j4 = MAX( j4, 0); */
	 j4 = (j4 - j3) * n;  j3 = (j3 - j2) * n;  j2 = (j2 - j1) *n;
	 }
 dx1 = 1.0 - dx0;  dx2 = -dx0 * 0.5;  dx3 = dx0 * dx2;
 dx4 = 3. * dx0 * dx3;
 b1 = dx2*dx1*dx1; b2 = 1.-dx4+5.*dx3; b3 = dx4-(dx2+4.*dx3); b4 = dx3*dx1;
 /* low corner offset */
 iq = i1 + j1 * n;
 switch (data_type) {
 case 0:
   bb.b = base.b+iq;
   xq = b1*(c1 * *(bb.b) + c2 * *(bb.b+i2)+ c3 * *(bb.b+i3) + c4 * *(bb.b+i4));
   bb.b += j2;
   xq += b2*(c1 * *(bb.b) + c2 * *(bb.b+i2)+ c3 * *(bb.b+i3) + c4 * *(bb.b+i4));
   bb.b += j3;
   xq += b3*(c1 * *(bb.b) + c2 * *(bb.b+i2)+ c3 * *(bb.b+i3) + c4 * *(bb.b+i4));
   bb.b += j4;
   xq += b4*(c1 * *(bb.b) + c2 * *(bb.b+i2)+ c3 * *(bb.b+i3) + c4 * *(bb.b+i4));
   /* byte arrays need to be range restricted, too many simple minds out there */
   xq = MAX( 0, MIN( 255, xq));
   /* also we need to round rather than truncate, taking that extra care */
   *out.b++ = rint(xq); break;
 case 1:
   bb.w = base.w+iq;
   xq = b1*(c1 * *(bb.w) + c2 * *(bb.w+i2)+ c3 * *(bb.w+i3) + c4 * *(bb.w+i4));
   bb.w += j2;
   xq += b2*(c1 * *(bb.w) + c2 * *(bb.w+i2)+ c3 * *(bb.w+i3) + c4 * *(bb.w+i4));
   bb.w += j3;
   xq += b3*(c1 * *(bb.w) + c2 * *(bb.w+i2)+ c3 * *(bb.w+i3) + c4 * *(bb.w+i4));
   bb.w += j4;
   xq += b4*(c1 * *(bb.w) + c2 * *(bb.w+i2)+ c3 * *(bb.w+i3) + c4 * *(bb.w+i4));
   /* also we need to round rather than truncate, taking that extra care */
   *out.w++ = rint(xq); break;
 case 2:
   bb.l = base.l+iq;
   xq = b1*(c1 * *(bb.l) + c2 * *(bb.l+i2)+ c3 * *(bb.l+i3) + c4 * *(bb.l+i4));
   bb.l += j2;
   xq += b2*(c1 * *(bb.l) + c2 * *(bb.l+i2)+ c3 * *(bb.l+i3) + c4 * *(bb.l+i4));
   bb.l += j3;
   xq += b3*(c1 * *(bb.l) + c2 * *(bb.l+i2)+ c3 * *(bb.l+i3) + c4 * *(bb.l+i4));
   bb.l += j4;
   xq += b4*(c1 * *(bb.l) + c2 * *(bb.l+i2)+ c3 * *(bb.l+i3) + c4 * *(bb.l+i4));
   /* also we need to round rather than truncate, taking that extra care */
   *out.l++ = rint(xq); break;
 case 3:
   //printf("float case, i1,i2,i3,i4,j1,j2,j3,j4 = %d %d %d %d, %d %d %d %d\n",i1,i2,i3,i4,j1,j2,j3,j4);
   bb.f = base.f+iq;
   xq = b1*(c1 * *(bb.f) + c2 * *(bb.f+i2)+ c3 * *(bb.f+i3) + c4 * *(bb.f+i4));
   bb.f += j2;
   xq += b2*(c1 * *(bb.f) + c2 * *(bb.f+i2)+ c3 * *(bb.f+i3) + c4 * *(bb.f+i4));
   bb.f += j3;
   xq += b3*(c1 * *(bb.f) + c2 * *(bb.f+i2)+ c3 * *(bb.f+i3) + c4 * *(bb.f+i4));
   bb.f += j4;
   xq += b4*(c1 * *(bb.f) + c2 * *(bb.f+i2)+ c3 * *(bb.f+i3) + c4 * *(bb.f+i4));
   //printf("xq = %f\n", xq);
   *out.f++ = xq; break;
 case 4:
   bb.d = base.d+iq;
   xq = b1*(c1 * *(bb.d) + c2 * *(bb.d+i2)+ c3 * *(bb.d+i3) + c4 * *(bb.d+i4));
   bb.d += j2;
   xq += b2*(c1 * *(bb.d) + c2 * *(bb.d+i2)+ c3 * *(bb.d+i3) + c4 * *(bb.d+i4));
   bb.d += j3;
   xq += b3*(c1 * *(bb.d) + c2 * *(bb.d+i2)+ c3 * *(bb.d+i3) + c4 * *(bb.d+i4));
   bb.d += j4;
   xq += b4*(c1 * *(bb.d) + c2 * *(bb.d+i2)+ c3 * *(bb.d+i3) + c4 * *(bb.d+i4));
   *out.d++ = xq; break;
 }
 return;
 }
 /*------------------------------------------------------------------------- */
void bicubic_fc()	/* internal routine for single pixel */
 {
 /* used by all (most?) routines that do bi-cubic interpolations, runs
 a little slower than the originals, perhaps because of the overhead
 in the call or other adjustments made */
 int	i1, i2, i3, i4, j1, j2, j3, j4, iq;
 float	c1, c2, c3, c4, b1, b2, b3, b4, dx0, dx1, dx2, dx3, dx4, xq, yq;
 union	types_ptr bb;
 /* the location is in xl, yl; base is the pointer to array; out is
 pointer to output; both are unions */
 i2 = (int) xl;		j2 = (int) yl;
 if ( i2 >= 1 && i2 < nm2 ) {		/* normal interior */
	 dx0 = xl - i2; i1 = i2 - 1; i2 = 1; i3 = 2; i4 = 3;	 }
	 else {				/* edge cases */
	 /* check if stretchmarks required */
	 if (stretchmark_flag == 0) {
	  if ( xl < -0.5 || xl > fnm5) {
	   switch (data_type) {
	   case 0: *out.b++ = 0; break;
	   case 1: *out.w++ = 0; break;
	   case 2: *out.l++ = 0; break;
	   case 3: *out.f++ = 0; break;
	   case 4: *out.d++ = 0; break;
	   }
	  return;
	  }
	 }
	 i2 = MIN(i2, nm1);	i2 = MAX( i2, 0);
	 xq = MIN(xl, fnm1);	xq = MAX(xq, 0);
	 dx0 = xq - i2;
 /*	printf("dx0 = %f\n",dx0);*/
	 i1 = MIN(i2-1, nm1);	i1 = MAX( i1, 0);
	 i3 = MIN(i2+1, nm1);	/* i3 = MAX( i3, 0); */
	 i4 = MIN(i2+2, nm1);	/* i4 = MAX( i4, 0); */
	 i2 = i2 - i1;	i3 = i3 - i1;	i4 = i4 - i1;
	 }

 dx1 = 1.0 - dx0;  dx4 = -dx0*dx1;  c1 = dx4 * dx1;  c4 = dx0 * dx4;
 dx2 = dx0 * dx0; dx3 = dx0 * dx2; c2 = 1.-2.*dx2+dx3; c3 = dx0*(1.0+dx0-dx2);
 if ( j2 >= 1 && j2 < mm2 ) {		/* normal interior */
	 j1 = j2 - 1; dx0 = yl - j2; j3 = n; j2 = n; j4 = n;	}
	 else {				/* edge cases */
	 /* check if stretchmarks required */
	 if (stretchmark_flag == 0) {
	  if ( yl < -0.5 || yl > fmm5) {
	   switch (data_type) {
	   case 0: *out.b++ = 0; break;
	   case 1: *out.w++ = 0; break;
	   case 2: *out.l++ = 0; break;
	   case 3: *out.f++ = 0; break;
	   case 4: *out.d++ = 0; break;
	   }
	  return;
	  }
	 }	 j2 = MIN(j2, mm1);	j2 = MAX( j2, 0);
	 xq = MIN(yl, fmm1);	xq = MAX(xq, 0);
	 dx0 = xq - j2;
	 j1 = MIN(j2-1, mm1);	j1 = MAX( j1, 0);
	 j3 = MIN(j2+1, mm1);	/* j3 = MAX( j3, 0); */
	 j4 = MIN(j2+2, mm1);	/* j4 = MAX( j4, 0); */
	 j4 = (j4 - j3) * n;  j3 = (j3 - j2) * n;  j2 = (j2 - j1) *n;
	 }

 dx1 = 1.0 - dx0;  dx4 = -dx0*dx1;  b1 = dx4 * dx1;  b4 = dx0 * dx4;
 dx2 = dx0 * dx0; dx3 = dx0 * dx2; b2 = 1.-2.*dx2+dx3; b3 = dx0*(1.0+dx0-dx2);
 /* low corner offset */
 iq = i1 + j1 * n;
 switch (data_type) {
 case 0:
   bb.b = base.b+iq;
   xq = b1*(c1 * *(bb.b) + c2 * *(bb.b+i2)+ c3 * *(bb.b+i3) + c4 * *(bb.b+i4));
   bb.b += j2;
   xq += b2*(c1 * *(bb.b) + c2 * *(bb.b+i2)+ c3 * *(bb.b+i3) + c4 * *(bb.b+i4));
   bb.b += j3;
   xq += b3*(c1 * *(bb.b) + c2 * *(bb.b+i2)+ c3 * *(bb.b+i3) + c4 * *(bb.b+i4));
   bb.b += j4;
   xq += b4*(c1 * *(bb.b) + c2 * *(bb.b+i2)+ c3 * *(bb.b+i3) + c4 * *(bb.b+i4));
   /* byte arrays need to be range restricted, too many simple minds out there */
   xq = MAX( 0, MIN( 255, xq));
   /* also we need to round rather than truncate, taking that extra care */
   *out.b++ = rint(xq); break;
 case 1:
   bb.w = base.w+iq;
   xq = b1*(c1 * *(bb.w) + c2 * *(bb.w+i2)+ c3 * *(bb.w+i3) + c4 * *(bb.w+i4));
   bb.w += j2;
   xq += b2*(c1 * *(bb.w) + c2 * *(bb.w+i2)+ c3 * *(bb.w+i3) + c4 * *(bb.w+i4));
   bb.w += j3;
   xq += b3*(c1 * *(bb.w) + c2 * *(bb.w+i2)+ c3 * *(bb.w+i3) + c4 * *(bb.w+i4));
   bb.w += j4;
   xq += b4*(c1 * *(bb.w) + c2 * *(bb.w+i2)+ c3 * *(bb.w+i3) + c4 * *(bb.w+i4));
   /* also we need to round rather than truncate, taking that extra care */
   *out.w++ = rint(xq); break;
 case 2:
   bb.l = base.l+iq;
   xq = b1*(c1 * *(bb.l) + c2 * *(bb.l+i2)+ c3 * *(bb.l+i3) + c4 * *(bb.l+i4));
   bb.l += j2;
   xq += b2*(c1 * *(bb.l) + c2 * *(bb.l+i2)+ c3 * *(bb.l+i3) + c4 * *(bb.l+i4));
   bb.l += j3;
   xq += b3*(c1 * *(bb.l) + c2 * *(bb.l+i2)+ c3 * *(bb.l+i3) + c4 * *(bb.l+i4));
   bb.l += j4;
   xq += b4*(c1 * *(bb.l) + c2 * *(bb.l+i2)+ c3 * *(bb.l+i3) + c4 * *(bb.l+i4));
   /* also we need to round rather than truncate, taking that extra care */
   *out.l++ = rint(xq); break;
 case 3:
   bb.f = base.f+iq;
   xq = b1*(c1 * *(bb.f) + c2 * *(bb.f+i2)+ c3 * *(bb.f+i3) + c4 * *(bb.f+i4));
   bb.f += j2;
   xq += b2*(c1 * *(bb.f) + c2 * *(bb.f+i2)+ c3 * *(bb.f+i3) + c4 * *(bb.f+i4));
   bb.f += j3;
   xq += b3*(c1 * *(bb.f) + c2 * *(bb.f+i2)+ c3 * *(bb.f+i3) + c4 * *(bb.f+i4));
   bb.f += j4;
   xq += b4*(c1 * *(bb.f) + c2 * *(bb.f+i2)+ c3 * *(bb.f+i3) + c4 * *(bb.f+i4));
   *out.f++ = xq; break;
 case 4:
   bb.d = base.d+iq;
   xq = b1*(c1 * *(bb.d) + c2 * *(bb.d+i2)+ c3 * *(bb.d+i3) + c4 * *(bb.d+i4));
   bb.d += j2;
   xq += b2*(c1 * *(bb.d) + c2 * *(bb.d+i2)+ c3 * *(bb.d+i3) + c4 * *(bb.d+i4));
   bb.d += j3;
   xq += b3*(c1 * *(bb.d) + c2 * *(bb.d+i2)+ c3 * *(bb.d+i3) + c4 * *(bb.d+i4));
   bb.d += j4;
   xq += b4*(c1 * *(bb.d) + c2 * *(bb.d+i2)+ c3 * *(bb.d+i3) + c4 * *(bb.d+i4));
   *out.d++ = xq; break;
 }
 return;
 }
 /*------------------------------------------------------------------------- */
