/*
 *  rdutil.c
 *
 *  library of utility functions for dealing with ring-diagram fit files
 *
 *    autoweed_vel	automatically weed velocities for inversions
 *    read_fit		reads n, l, f, and d_f from a fit file
 *    read_fit_v	reads n, l, f, d_f, u_i and d_ui from a fit file
 *
 *  Bugs:
 *    read_fit() and read_fit_v() assume that the frequency error estimate
 *	is in column 5 of the fit table; for fit files produced with versions
 *	of rdfitc below 1.2, this is not the case
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

	/*  some or all of the following declarations may be unnecessary  */
/*
int interp_vel(int *n, double *l, double *f, double *ef, double *ux,
    double *eux, double *uy, double *euy, int npts);
int interp_freq(int *n, double *l, double *f, double *ef, int npts,
    int *ni, int *li, double *fi, double *efi, int *npts_out);
int freqdif(int *n1, int *n2, int *l1, double *l2, double *f1, double *f2,
    double *ef1, double *ef2, int npts1, int npts2, int *n, int *l, double *f,
    double *df, double *edf, int *npts);
int autoweed(int *l, int *n, double *f, double *df, double *edf, int *msk,
    int npts);

int nearst(double xb, double x[], int ntab);
int divdif(double xb, double x[], double f[], int *nuse, int ntab,
    double fb[], double aeps, double *dfb, double *ddfb);
int svdevl(int n, int m, double *u, double *v, double sigma[], int lu,
    int lv, double b[], double reps);
int svd(int n, int m, double *a, double *v, double sigma[], int la, int lv);
int llsq(int n, int m, int k, double *x, int ix, double f[], double ef[],
    double a[], double *u, double *v, int iu, int iv, double sigma[],
    double y[], void (* phi) (int , double * , double * ), double reps,
    double *chisq);
void line(int m, double *x, double *y);
*/
/*
 * The following have been commented out for dependencies and because they are
 *    not used by rdfitc (which only uses autoweed_vel())
 *	interp_vel
 *	interp_freq
 *	freqdif
 *	autoweed
 */
int autoweed_vel (int *n, double *l, double *ux, double *uy, int *mask,
    int npts) {
/*
 *  Returns 'mask' of length npts: mask[i] is False for rejected modes.
 *    Ux and Uy are weeded together. Values are rejected when their difference
 *    from the mean exceeds 5 (tol?) std. dev., both the mean and std dev
 *    being calculated over the central 3/5 of the l range for the order
 */
  const int maxn = 8;
  const double tol = 5.0;

  int i, j, offset;
  int n_num[maxn];
  double num;
  double sumx, sumxx, sumy, sumyy, meanx, meany, stdx, stdy;
  double range, uxmin, uxmax, uymin, uymax;
  double lmin,lmax;
	  
  for (i=0; i<maxn; i++) n_num[i] = 0;
  for (i=0; i<npts; i++)
    if (n[i] < maxn) n_num[n[i]]++;

  offset=0;
  for (i=0; i<maxn; i++) {
    num = 0.0;
    sumx = 0.0;
    sumy = 0.0;
    sumxx = 0.0;
    sumyy = 0.0;
							    /*  get l limits  */
    lmin = lmax = l[offset];
    for (j=offset; j < n_num[i] + offset; j++) {
      if (l[j] < lmin) lmin = l[j];
      if (l[j] > lmax) lmax = l[j];
    }
    range = 0.2 * (lmax - lmin);
    lmin += 2.0 * range;
    lmax -= 2.0 * range;
    for (j=offset; j < n_num[i] + offset; j++) {
      if ((l[j] <= lmax) && (l[j] >= lmin)) {
	sumx += ux[j];
	sumxx += ux[j]*ux[j];
	sumy += uy[j];
	sumyy += uy[j]*uy[j];
	num++;
      }
    }
    meanx = sumx / num;
    meany = sumy / num;
    stdx = num*sumxx - sumx*sumx;
    stdy = num*sumyy - sumy*sumy;
    stdx /= (num * (num - 1));
    stdy /= (num * (num - 1));
    stdx = sqrt (stdx);
    stdy = sqrt (stdy);
    uxmin = meanx - 5.0 * stdx;
    uxmax = meanx + 5.0 * stdx;
    uymin = meany - 5.0 * stdy;
    uymax = meany + 5.0 * stdy;
    
    for (j=offset; j<n_num[i]+offset; j++) {
      if (ux[j] <= uxmax && ux[j] >= uxmin &&
	  uy[j] <= uymax && uy[j] >= uymin) mask[j]=1;
      else mask[j]=0;
    }
    offset +=n_num[i];
  }
  return 0;
}

int read_fit (FILE *fpt, int *npts, int **n, double **l, double **f,
    double **ef) {
  int i, nlines;
  char buffer[8192];
  if (ferror (fpt)) return -1;
  if (feof (fpt)) rewind (fpt);

  nlines = 0;
  while (!feof(fpt)) {
    fgets (buffer, 8192, fpt);
    if (buffer[0] != '#' && !feof (fpt)) nlines++;
  }
  if (nlines == 0) return 1;
  
  *n = (int *)malloc (nlines * sizeof (int));
  *l = (double *)malloc (nlines * sizeof (double));
  *f = (double *)malloc (nlines * sizeof (double));
  *ef = (double *)malloc (nlines * sizeof (double));
  
  rewind (fpt);
  
  for (i=0; i<nlines; i++) {
    fgets (buffer, 8192, fpt);
    while (buffer[0] == '#') fgets( buffer, 8192, fpt);
    sscanf (buffer, "%i %lf %*f %lf %lf", &(*n)[i], &(*l)[i], &(*f)[i], 
        &(*ef)[i]);
  } 
  *npts = nlines;
  
  return 0;
}

int read_fit_v (FILE *fpt, int *npts, int **n, double **l, double **f,
    double **ef, double **ux, double **eux, double **uy, double **euy) {
  int i, nlines;
  char buffer[8192];
  if (ferror (fpt)) return -1;
  if (feof (fpt)) rewind (fpt);
  
  nlines = 0;
  while (!feof(fpt)) {
    fgets (buffer, 8192, fpt);
    if (buffer[0] != '#' && !feof (fpt)) nlines++;
  }
  if (nlines == 0) return 1;
  
  *n = (int *)malloc (nlines * sizeof (int));
  *l = (double *)malloc (nlines * sizeof (double));
  *f = (double *)malloc (nlines * sizeof (double));
  *ef = (double *)malloc (nlines * sizeof (double));
  *ux = (double *)malloc (nlines * sizeof (double));
  *eux = (double* )malloc (nlines * sizeof (double));
  *uy = (double *)malloc (nlines * sizeof (double));
  *euy = (double *)malloc (nlines * sizeof (double));
  
  rewind (fpt);
  
  for (i=0; i<nlines; i++) {
    fgets (buffer, 8192, fpt);
    while (buffer[0] == '#') fgets (buffer, 8192, fpt);
    sscanf (buffer, "%i %lf %*f %lf %lf %lf %lf %lf %lf", &(*n)[i], &(*l)[i],
        &(*f)[i], &(*ef)[i], &(*ux)[i], &(*eux)[i], &(*uy)[i], &(*euy)[i]);
  } 
  *npts = nlines;
  
  return 0;
}
/*
int interp_vel(int *n, double *l, double *f, double *ef, double *ux, 
      double *eux, double *uy, double *euy, int npts)	{
*/
  /*
     Function - interp
     Takes power spectra fit parameters, and interpolates them to integer l.
     Data are returned in the input arrays - **data is overwritten!**
  */
/*
  int i, j, nuse, ierr, offset=0;
  int n_num[13];
  double fb[10], ll;
  double dfb, ddfb, aeps=1e-7;
  double *inp_l, *inp_f, *inp_ef, *inp_ux, *inp_uy, *inp_eux, *inp_euy;
  FILE * outBug;
  
  // count number of frequencies for each n
  // note that it is assumed that all frequencies for each 
  //	n are grouped together in the arrays
  for(i=0; i<13; i++) n_num[i]=0;
  for(i=0; i < npts; i++)	{
	  n_num[n[i]] ++;
  }
  
  inp_l=(double*) malloc(npts*sizeof(double));
  inp_f=(double*) malloc(npts*sizeof(double));
  inp_ef=(double*) malloc(npts*sizeof(double));
  inp_ux=(double*) malloc(npts*sizeof(double));
  inp_uy=(double*) malloc(npts*sizeof(double));
  inp_eux=(double*) malloc(npts*sizeof(double));
  inp_euy=(double*) malloc(npts*sizeof(double));
  // loop over n
//  outBug=fopen("tesbug","w");
  for(i=0; i < 10; i++)	{
    for(j=0;j<npts;j++)	{
      inp_l[j]=0.0;
      inp_f[j]=0.0;
      inp_ef[j]=0.0;
      inp_ux[j]=0.0;
      inp_uy[j]=0.0;
      inp_eux[j]=0.0;
      inp_euy[j]=0.0;
    }
    for(j=0;j<n_num[i];j++)	{
      inp_l[j]=l[j+offset];
      inp_f[j]=f[j+offset];
      inp_ef[j]=ef[j+offset];
      inp_ux[j]=ux[j+offset];
      inp_uy[j]=uy[j+offset];
      inp_eux[j]=eux[j+offset];
      inp_euy[j]=euy[j+offset];
    }
  for(j=0;j<n_num[i];j++)	{
    ll=(int)(inp_l[j]+0.5);
    nuse=3;
    ierr=divdif(ll, inp_l, inp_f, &nuse, 
	n_num[i], fb, aeps, &dfb, &ddfb);
    f[j+offset]=fb[nuse];
    nuse=3;
    ierr=divdif(ll, inp_l, inp_ef, &nuse,
	n_num[i], fb, aeps, &dfb, &ddfb);
    ef[j+offset]=fb[nuse];
    nuse=3;
    ierr=divdif(ll, inp_l, inp_ux, &nuse, 
	n_num[i], fb, aeps, &dfb, &ddfb);
    ux[j+offset]=fb[nuse];
    nuse=3;
    ierr=divdif(ll, inp_l, inp_uy, &nuse, 
	n_num[i], fb, aeps, &dfb, &ddfb);
    uy[j+offset]=fb[nuse];
    nuse=3;
    ierr=divdif(ll, inp_l, inp_eux, &nuse, 
	n_num[i], fb, aeps, &dfb, &ddfb);
    eux[j+offset]=fb[nuse];
    nuse=3;
    ierr=divdif(ll, inp_l, inp_euy, &nuse, 
	n_num[i], fb, aeps, &dfb, &ddfb);
    euy[j+offset]=fb[nuse];
    l[j+offset]=ll;
      }
    offset+=n_num[i];
  }
  
  return 0;
}
*/
/*
int interp_freq(int *n, double *l, double *f, double *ef, int npts, 
    int *ni, int *li, double *fi, double *efi, int *npts_out)	{
  int i, j, nn, ll, nmin, nmax, this_npts, nuse;
  double *this_l, *this_f, *this_ef, all, fb[5];
  double reps, dfb, ddfb, ff, error;
  
  this_l = (double*) malloc(npts*sizeof(double));
  this_f = (double*) malloc(npts*sizeof(double));
  this_ef = (double*) malloc(npts*sizeof(double));

  *npts_out = 0;
  nmin = 100;
  nmax = 0;
  for(i=0; i<npts; i++)	{
    if(n[i] < nmin) nmin = n[i];
    if(n[i] > nmax) nmax = n[i];
  }
  for(nn=nmin; nn<=nmax; nn++)	{
    j = 0;
    for(i=0; i<npts; i++)	{
      if(n[i]==nn)	{
	this_l[j] = l[i];
	this_f[j] = f[i];
	this_ef[j] = ef[i];
	j++;
      }
    }
    this_npts = j;
    j = *npts_out;
    //printf("%i %i\n",this_npts, *npts_out);
    for(i=0; i<this_npts; i++)	{
      ll = (int) this_l[i] + 0.5;
      all = (double) ll;	// cast back to double for function argument
      nuse = 3;
      reps = 1.0e-7;
    //printf("got here\n");
      divdif(all, this_l, this_f, &nuse, this_npts, fb, reps, &dfb, &ddfb);
      ff = fb[nuse];
      nuse = 3;
      divdif(all, this_l, this_ef, &nuse, this_npts, fb, reps, &dfb, &ddfb);
      error = fb[nuse];
      //printf("%i %i %f %f\n",nn,ll,ff,error);
      if(error>0)	{
	ni[j] = nn;
	li[j] = ll;
	fi[j] = ff;
	efi[j] = error;
	j++;
      }
    }
    //printf("%i\n",j);
    *npts_out = j;
  }

  return 0;
}
*/
/*
int freqdif(int *n1, int *n2, int *l1, double *l2, double *f1, double *f2, 
    double *ef1, double *ef2, int npts1, int npts2, int *n, int *l, double *f,
    double *df, double *edf, int *npts)	{
  int i, j, nn, ll, *index, nmin, nmax, lmin, lmax, this_npts, ind;
  int nuse;
  double *this_l, *this_f, *this_ef, all, fb[5], dfb, ddfb, reps;
  double ff, error;

  *npts = 0;

  this_l = (double*) malloc(npts2*sizeof(double));
  this_f = (double*) malloc(npts2*sizeof(double));
  this_ef = (double*) malloc(npts2*sizeof(double));

  nmin = n1[0];
  nmax = n1[0];
//  lmin = (int) l2[0]+0.5;
//  lmax = (int) l2[0]+0.5;
  for(i=1; i<npts1; i++)	{
    if(n1[i] < nmin) nmin = n1[i];
    if(n1[i] > nmax) nmax = n1[i];
//    if(l2[i] < lmin) lmin = (int) l2[i]+0.5;
//    if(l2[i] > lmax) lmax = (int) l2[i]+0.5;
  }

//  index = (int*) malloc((nmax-nmin)*(lmax-lmin)*sizeof(int));
//  for(i=0; i<(nmax-nmin)*(lmax-lmin); i++) index[i] = -1;
//  for(i=0; i<npts2; i++)	{
//    ll = (int) l2[i]+0.5;
//    index[n2[i]*(lmax-lmin)+ll-lmin] = i;
//  }

  for(nn=nmin; nn<=nmax; nn++)	{
    j = 0;
    lmin = 10000;
    lmax = 0;
    for(i=0; i<npts2; i++)	{
      if(n2[i] == nn)	{
	this_l[j] = l2[i];
	this_f[j] = f2[i];
	this_ef[j] = ef2[i];
	if(lmin > l2[i]) lmin = (int) l2[i]+0.5;
	if(lmax < l2[i]) lmax = (int) l2[i]+0.5;
	j++;
      }
    }
    this_npts = j;
    //printf("this_npts = %i\n", this_npts);
    j = *npts;
    for(i=0; i<npts1; i++)	{
//      ind = index[n1[i]*(lmax-lmin)+l1[i]-lmin];
//      if(ind>=0)	{
      if(n1[i] == nn && l1[i] >= lmin && l1[i] <= lmax)	{
	//ll = (int) this_l[i] + 0.5;
	ll = l1[i];
	all = (double) ll;	// cast back to double for function argument
	nuse = 3;
	reps = 1.0e-7;
	divdif(all, this_l, this_f, &nuse, this_npts, fb, reps, &dfb, &ddfb);
	ff = fb[nuse];
	nuse = 3;
	divdif(all, this_l, this_ef, &nuse, this_npts, fb, reps, &dfb, &ddfb);
	error = fb[nuse];
	if(error>0.0)	{
	  n[j] = nn;
	  l[j] = ll;
	  f[j] = f1[i];
	  df[j] = f1[i] - ff;
	  edf[j] = sqrt(ef1[i]*ef1[i]+error*error);
	  //printf("%i %i  %i  %f  %f  %f\n",j, n[j], l[j], f[j], df[j], edf[j]);
	  j++;
	}
      }
    }
    *npts = j;
    //printf("%i\n",*npts);
  }

  return 0;
}
*/
void line(int m, double *x, double *y)	{
  y[0] = x[0];
  y[1] = 1.0;
  return;
}
/*
int autoweed(int *l, int *n, double *f, double *df, double *edf, int *msk,
    int npts)	{
  int avglen = 4;
  double tol = 4.0;
  double delnu = 50.0;
  int i, j, nn, this_npts, llim, ulim, offset;
  int data_range[2], *this_l;
  double *this_f, *this_df, *this_edf, *slopes, *a, *u, *y;
  double v[4], sigma[2], chisq;
  double mean, std;

  this_l = (int*) malloc(npts*sizeof(int));
  this_f = (double*) malloc(npts*sizeof(double));
  this_df = (double*) malloc(npts*sizeof(double));
  this_edf = (double*) malloc(npts*sizeof(double));
  slopes = (double*) malloc(npts*sizeof(double));
  a = (double*) malloc(avglen*sizeof(double));
  u = (double*) malloc(2*avglen*sizeof(double));
  y = (double*) malloc(avglen*sizeof(double));
  offset = 0;
  for(i=0; i<npts; i++) msk[i] = 0;
  for(nn=0; nn<7; nn++)	{
    j = 0;
    mean = 0.0;
    std = 0.0;
    for(i=0; i<npts; i++)	{
      if(n[i] == nn)	{
	this_l[j] = l[i];
	this_f[j] = f[i];
	this_df[j] = df[i];
	this_edf[j] = edf[i];
	j++;
      }
    }
    this_npts = j;
    if(this_npts > avglen + 5)	{
      data_range[0] = (int) 2*this_npts/10.;
      data_range[1] = (int) 8*this_npts/10.;
      for(i=0; i<this_npts-avglen; i++)	{
	llsq(avglen, 2, 1, &this_f[i], 1, &this_df[i], &this_edf[i], a, u, 
	    v, 2, 2, sigma, y, *line, 1.0e-6, &chisq);
	slopes[i] = a[0];
	//printf("%i  %f\n",nn,a[0]);
	mean += a[0];
      }
      mean /= (double) this_npts;
      for(i=0; i<this_npts-avglen; i++) std += (slopes[i]-mean)*(slopes[i]-mean);
      std /= (double) this_npts;
      std = sqrt(std);
//      printf("%i %i %i\n",data_range[0], data_range[1], this_npts/2);
      llim = data_range[0];
      ulim = data_range[1];
      while(llim > 0)	{
	if( fabs(slopes[llim]-mean)/std > tol) break;
	llim--;
      }
      while(ulim<this_npts-avglen-1)	{
	if( fabs(slopes[ulim]-mean)/std > tol) break;
	ulim++;
      }
      i = this_npts/2;
      while(i>0)	{
	if(this_f[i]-this_f[i-1] > delnu) break;
	i--;
      }
      if(i > llim) llim = i;
      i = this_npts/2;
      while(i<this_npts)	{
	if(this_f[i]-this_f[i-1] > delnu) break;
	i++;
      }
      if(i < ulim) ulim = i;
      printf("%i %i %i %i\n",nn,this_npts,llim,ulim);
      for(i=offset+llim; i<offset+ulim; i++) msk[i] = 1;//i<this_npts+offset-ulim; i++) msk[i] = 1;
      offset += this_npts;
    }
  }

  return 0;
}
*/

int nearst(double xb, double x[], int ntab)	{
  int low, igh, mid;

  low=0; igh=ntab-1;
  if((xb < x[low]) != (xb < x[igh]) ) {

/*	If the point is within the range of table, then locate it by bisection */

    while(igh-low > 1) {
      mid=(low+igh)/2;
      if((xb < x[mid]) == (xb < x[low])) low=mid;
      else igh=mid;
    }
  }

  if(fabs(xb-x[low]) < fabs(xb-x[igh])) return low;
  else return igh;
}

int divdif(double xb, double x[], double f[], int *nuse, int ntab,
    double fb[], double aeps, double *dfb, double *ddfb)	{
  int i,j,k,next,in,ip,nit,ier, nmax=10;
  double err,px,dpx,ddpx,xn[11],xd[11];

/*	Find the nearest point */

  next=nearst(xb,x,ntab);
  fb[1]=f[next];
  xd[1]=f[next];
  xn[1]=x[next];
  ier=0;
  px=1.0;

/*	Initialisation for the derivatives */

  *dfb=0.0; *ddfb=0.0;
  dpx=0.0; ddpx=0.0;

/*	Points between IN and IP are used for interpolation */

  ip=next; in=next;

/*	Maximum number of points to be used for interpolation */
  nit=*nuse; if(nmax<nit) nit=nmax; if(ntab<nit) nit=ntab;
  if(*nuse>nmax || *nuse>ntab) ier=22;
  if(*nuse<1) {
    ier=21;
    nit=6; if(nmax<nit) nit=nmax; if(ntab<nit) nit=ntab;
  }
  *nuse=1;
	  
/*	Calculate successive interpolation polynomial */
  for(j=2; j<=nit; ++j) {
/*	Choose the next nearest point to XB */
    if(in<=0 ) {
      ip=ip+1; next=ip;
    }
    else if(ip >= ntab-1) {
      in=in-1; next=in;
    }
    else if(fabs(xb-x[ip+1]) < fabs(xb-x[in-1]) ) {
      ip=ip+1; next=ip;
    }
    else {
      in=in-1; next=in;
    }

/*	Calculating the divided differences */
    xd[j]=f[next];
    xn[j]=x[next];
    for(k=j-1; k>=1; --k) xd[k]=(xd[k+1]-xd[k])/(xn[j]-xn[k]);

/*	Calculating the derivatives */
    ddpx=ddpx*(xb-xn[j-1])+2.*dpx;
    dpx=dpx*(xb-xn[j-1])+px;
    *dfb = *dfb+dpx*xd[1];
    *ddfb = *ddfb+ddpx*xd[1];

    px=px*(xb-xn[j-1]);
    err=xd[1]*px;
    fb[j]=fb[j-1]+err;
    *nuse=j;

    if(fabs(err) < aeps) return ier;
  }
  return 23;
}

int svd(int n, int m, double *a, double *v, double sigma[], int la, int lv)

{
	int i,j,k,l,itr,ier, itmax=30;
	double f, g, h, rmax, s, r1, r2, c, x, y, z, aeps, eps=1.e-16;
//	double *e;
	double e[2];

	if(n>m || n<=0 || m<=0 || n>la || n>lv) return 105;
	ier=0;

/*	Reduction to Bidiagonal form using Householder transformations */
	g=0.0; rmax=0.0;
//	e=(double *) calloc(n, sizeof(double));

	for(i=0; i<n; ++i) {
/*	Off-diagonal elements of bidiagonal form  */
		e[i]=g;
		s=0.0;
		for(j=i; j<m; ++j) s=s+a[i+j*la]*a[i+j*la];
		if(s <= 0.0) {
/*	transformation not required */
			g=0.0;
		}
		else {
			f= a[i+i*la];
			g=sqrt(s);
			if(f>=0.0) g=-g;
			h=f*g-s;
			a[i+i*la] = f-g;

			for(j=i+1; j<n; ++j) {
				s=0.0;
				for(k=i; k<m; ++k) s=s+a[i+k*la]*a[j+k*la];
				f=s/h;
				for(k=i; k<m; ++k) a[j+k*la]= a[j+k*la]+f*a[i+k*la];
			}
		}

/*	Diagonal elements of bidiagonal form  */
		sigma[i]=g;
		s=0.0;
		for(j=i+1; j<n; ++j) s=s+a[j+i*la]*a[j+i*la];

		if(s<= 0.0) g=0.0;
		else {
			f= a[i*la+(i+1)];
			g=sqrt(s);
			if(f>= 0.0) g=-g;
			h=f*g-s;
			a[i*la+(i+1)]=f-g;
			for(j=i+1; j<n; ++j) e[j]=a[j+i*la]/h;

			for(j=i+1; j<m; ++j) {
				s=0.0;
				for(k=i+1; k<n; ++k) s=s+a[k+j*la]*a[k+i*la];
				for(k=i+1; k<n; ++k) a[k+j*la] = a[k+j*la]+s*e[k];
			}
		}
		r1=fabs(sigma[i])+fabs(e[i]);
		if(r1 > rmax) rmax=r1;
	}

/*	Accumulation of right hand transformation in array V */
	for(i=n-1; i>=0; --i) {
		if(g != 0.0) {
			h=g*a[i*la+(i+1)];
			for(j=i+1; j<n; ++j) v[i+j*lv]=a[j+i*la]/h;

			for(j=i+1; j<n; ++j) {
				s=0.0;
				for(k=i+1; k<n; ++k) s=s+a[k+i*la]*v[j+k*lv];
				for(k=i+1; k<n; ++k) v[j+k*lv]=v[j+k*lv]+s*v[i+k*lv];
			}
		}

		for(j=i+1; j<n; ++j) {
			v[j+i*lv]=0.0; v[i+j*lv]=0.0;
		}
		v[i+i*lv]=1;
		g= e[i];
	}

/*	Accumulation of left hand transformation overwritten on matrix A */
	for(i=n-1; i>=0; --i) {
		g=sigma[i];
		for(j=i+1; j<n; ++j) a[j+i*la]=0.0;
		if(g != 0.0) {
			h=g*a[i+i*la];

			for(j=i+1; j<n; ++j) {
				s=0.0;
				for(k=i+1; k<m; ++k) s=s+a[i+k*la]*a[j+k*la];
				f=s/h;
				for(k=i; k<m; ++k) a[j+k*la]=a[j+k*la]+f*a[i+k*la];
			}

			for(j=i; j<m; ++j) a[i+j*la]=a[i+j*la]/g;
		}
		else {
			for(j=i; j<m; ++j) a[i+j*la]=0.0;
		}
		a[i+i*la] = a[i+i*la]+1;
	}

/*	Diagonalisation of the bidiagonal form */
	aeps=eps*rmax;
/*	Loop over the singular values */
	for(k=n-1; k>=0; --k) {
/*	The QR transformation */
		for(itr=1; itr<=itmax; ++itr) {

/*	Test for splitting */
			for(l=k; l>=0; --l) {
				if(fabs(e[l]) < aeps) goto split;
				if(fabs(sigma[l-1]) < aeps) break;
			}

/*	cancellation of E[L] if L>1  */
			c=0.0; s=1.0;
			for(i=l; i<=k; ++i) {
				f=s*e[i];
				e[i] = c*e[i];
				if(fabs(f) < aeps) goto split;
				g=sigma[i];
				sigma[i]=sqrt(f*f+g*g);
				c=g/sigma[i];
				s=-f/sigma[i];

				for(j=0; j<m; ++j) {
					r1= a[j*la+(l-1)];
					r2= a[i+j*la];
					a[j*la+(l-1)]=r1*c+r2*s;
					a[i+j*la]=c*r2-s*r1;
				}
			}

split:			z=sigma[k];
			if(l == k) {
/*	QR iteration has converged */
				if(z < 0.0) {
					sigma[k] = -z;
					for(j=0; j<n; ++j) v[k+j*lv]=-v[k+j*lv];
				}
				break;
			}

			if(itr==itmax) {ier=12; break;}
					
/*	calculating shift from bottom 2x2 minor */
			x=sigma[l];
			y=sigma[k-1];
			g=e[k-1];
			h=e[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.*h*y);
			g=sqrt(1.+f*f);
			if(f < 0.0) g=-g;
			f=((x-z)*(x+z)+h*(y/(f+g)-h))/x;

/*	next QR transformation */
			c=1.0; s=1.0;
/*	Given's rotation  */
			for(i=l+1; i<=k; ++i) {
				g=e[i];
				y=sigma[i];
				h=s*g;
				g=c*g;
				e[i-1]=sqrt(f*f+h*h);
				c=f/e[i-1];
				s=h/e[i-1];
				f=c*x+s*g;
				g=c*g-s*x;
				h=s*y;
				y=c*y;

				for(j=0; j<n; ++j) {
					x=v[j*lv+(i-1)];
					z=v[i+j*lv];
					v[j*lv+(i-1)]=c*x+s*z;
					v[i+j*lv]=c*z-s*x;
				}

				sigma[i-1]=sqrt(f*f+h*h);
				if(sigma[i-1] != 0.0) {
					c=f/sigma[i-1];
					s=h/sigma[i-1];
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for(j=0; j<m; ++j) {
					y= a[j*la+(i-1)];
					z= a[i+j*la];
					a[j*la+(i-1)] = c*y+s*z;
					a[i+j*la] = c*z-s*y;
				}
			}

			e[l]=0.0;
			e[k]=f;
			sigma[k]=x;
		}
	}
//	free(e);
	return ier;
}

int svdevl(int n, int m, double *u, double *v, double sigma[], int lu,
	int lv, double b[], double reps)

{
	int i,j;
	double smax, aeps, s;
	double *wk;

/*	Finding the largest singular value */
	smax=0.0;
	for(i=0; i<n; ++i)
		if(sigma[i] > smax) smax=sigma[i];

	aeps=smax*reps;
	wk=(double *)calloc(n, sizeof(double));
	for(i=0; i<n; ++i) {
		s=0.0;
/*	Only SIGMA[I] > AEPS contribute to the solution */
		if(sigma[i] > aeps) {
			for(j=0; j<m; ++j) s=s+b[j]*u[i+j*lu];
			s=s/sigma[i];
		}
		wk[i]=s;
	}

	for(i=0; i<n; ++i) {
		s=0.0;
		for(j=0; j<n; ++j) s=s+v[j+i*lv]*wk[j];
		b[i]=s;
	}
	free(wk);
	return 0;
}

int llsq(int n, int m, int k, double *x, int ix, double f[], double ef[],
	double a[], double *u, double *v, int iu, int iv, double sigma[],
	double y[], void (* phi) (int , double * , double * ), double reps,
	double *chisq)

{
	int i,j,ier;
	double s1;
//	double *wk;
	double wk[2];

	if(m>n || m<=0 || n<=0 || k>ix) return 606;
	
//	wk=(double *) calloc(m, sizeof(double));
/*	Setting up the design matrix and the RHS */
	for(i=0; i<n; ++i) {
		if(ef[i]<=0.0) {/*free(wk); */return 607;}
		a[i]=f[i]/ef[i];
		phi(m,&x[i*ix],wk);
		for(j=0; j<m; ++j) u[j+i*iu]=wk[j]/ef[i];
	}

	ier=svd(m,n,u,v,sigma,iu,iv);
	if(ier>100) {/*free(wk);*/ return ier;}

/*	Calculate the least squares solution */
	ier=svdevl(m,n,u,v,sigma,iu,iv,a,reps);
 
/*	Computing the \chi^2 from fitted coefficients */
	*chisq=0.0;
	for(i=0; i<n; ++i) {
		phi(m,&x[i*ix],wk);
		s1=0.0;
		for(j=0; j<m; ++j) s1=s1+a[j]*wk[j];
		y[i]=s1;
		s1=(f[i]-y[i])/ef[i];
		*chisq=(*chisq)+s1*s1;
	}

//	free(wk);
	return 0;
}
