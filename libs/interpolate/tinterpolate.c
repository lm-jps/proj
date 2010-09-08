// Based on coeffj.c
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <mkl.h>
#include "tinterpolate.h"
#include "drms_defs.h"

#define minval(x,y) (((x) < (y)) ? (x) : (y))
#define maxval(x,y) (((x) > (y)) ? (x) : (y))

int tinterpolate(
  int nsample, // Number of input times
  double *tsample, // Input times
  double tint, // Target time
  int nconst, // Number of polynomial terms exactly reproduced
  float **images, // Pointer array to input images
  unsigned char **masks, // Pointer array to input masks. 0 is good, 1 is missing
  float *image_out, // Interpolated image
  int nx, // Number of points in dimension adjacent in memory
  int ny, // Number of points in dimension not adjacent in memory
  int nlead, // Leading dimension of arrays. nlead>=nx
  int method, // Interpolation method
  char **filenamep, // Pointer to name of file to read covariance from.
                   // Set to actual file used if method > 0.
  float fillval // Value to use if not enough points present
)
{
  const unsigned char maskgood=0; // Value of good entries in masks
  int malign=32;
  int ncor=5000;
  //double dta=1.0; // Time between samples of input autocorrelation
  int i,j,isample;
  FILE *fileptr;
  double ti,dt,dtx;
  double *acor; // Input autocorrelation
  double *a0,*a,*rh0,*coeffd;
  double *a1b,bta1b,*a1r,bta1r,c;
  int idt;
  int info;
  char uplo[] = "U";
  int ione = 1;
  float *coeff;
  int nmasks/*,smask*/,imask,imask1,ngood;
  int *wgood;
  int *smasks,*ngoods;
  float *coeffs,*cp;
  float *sums,*ip;
  double t0,t1,t2;
  unsigned char *mp;
  //int ib, nblock;
  char *filename;
  double sum;
  double *xc,*b,*rhc,*bta1b1,*help;
  int iconst,jconst;

  if ((nsample < 1) || (nsample > 20)) {
// Code breaks at 31 or 32 due to the use of summed masks
// Also, memory usage and efficiency drops due to all masks calculated for >20
    printf("tinterpolate: must have 1<=nsample<=20\n");
    return 1;
  }
  t0=dsecnd();
  switch (method) {
  case 0:
    filename=strdup(*filenamep);
    break;
  case 1:
    filename=strdup(DEFS_MKPATH("/data/acort_hr12.txt"));
    break;
  default:
    printf("Unknown method in init_fill.\n");
    return 1;
    break;
  }
  fileptr = fopen (filename,"r");
  if (fileptr==NULL) {
    printf("File not found in tinterpolate.\n");
    return 1;
  }
  if (method != 0) {
    if (filenamep != NULL) {
      *filenamep=strdup(filename);
    }
  }

// Read autocorrelation
  acor=(double *)(MKL_malloc(ncor*sizeof(double),malign));
//fread ((char*)(acor),sizeof(double),ncor,fileptr);
  for (i=0;i<ncor;i++) fscanf(fileptr,"%lf",acor+i);
  fclose(fileptr);
  t0=dsecnd()-t0;

  t1=dsecnd();
  a0=(double *)(MKL_malloc(nsample*nsample*sizeof(double),malign));
  a=(double *)(MKL_malloc(nsample*nsample*sizeof(double),malign));
  rh0=(double *)(MKL_malloc(nsample*sizeof(double),malign));
  coeffd=(double *)(MKL_malloc(nsample*sizeof(double),malign));
  coeff=(float *)(MKL_malloc(nsample*sizeof(float),malign));
  a1b=(double *)(MKL_malloc(nconst*nsample*sizeof(double),malign));
  a1r=(double *)(MKL_malloc(nsample*sizeof(double),malign));
  xc=(double *)(MKL_malloc(nsample*sizeof(double),malign));
  b=(double *)(MKL_malloc(nconst*nsample*sizeof(double),malign));
  rhc=(double *)(MKL_malloc(nconst*sizeof(double),malign));
  help=(double *)(MKL_malloc(nconst*sizeof(double),malign));
  bta1b1=(double *)(MKL_malloc(nconst*nconst*sizeof(double),malign));

// nmasks=2^nsample
  nmasks=1;
  for (i=0;i<nsample;i++) nmasks=2*nmasks;
  coeffs=(float *)(MKL_malloc(nsample*nmasks*sizeof(float),malign));
  wgood=(int *)(MKL_malloc(nsample*sizeof(int),malign));

  for (i=0;i<nsample;i++) {
    ti=tsample[i];
    for (j=0;j<nsample;j++) {
      dt=fabs(tsample[j]-ti);
      idt=(float)floor(dt);
      dtx=dt-idt;
      a0[i*nsample+j]=(1.0-dtx)*acor[idt]+dtx*acor[idt+1];
    }
    dt=fabs(tint-ti);
    idt=(float)floor(dt);
    dtx=dt-idt;
    rh0[i]=(1.0-dtx)*acor[idt]+dtx*acor[idt+1];
  }

// For the time being calculate weights for all possible masks
// Could be done implicitly if nsample becomes large
// For 4096x4096 images interpolation dominates to around nsample=20
  for (imask=0;imask<nmasks;imask++) {
    cp=coeffs+nsample*imask;
    imask1=imask;
    ngood=0;
// Find good and bad given mask
    for (i=0;i<nsample;i++) {
      if ((imask1&1)==0) { // Good data
        wgood[ngood]=i;
        ngood=ngood+1;
      }
      imask1=imask1/2;
      cp[i]=0.0f; // Zero all coefficients
    } // for i
    if (ngood>0) {
// Pick out relevant parts of a0 and rh0
      for (i=0;i<ngood;i++) {
        for (j=0;j<ngood;j++) a[ngood*i+j]=a0[nsample*wgood[i]+wgood[j]];
        coeffd[i]=rh0[wgood[i]];
      }
      dpotrf(uplo,&ngood,a,&ngood,&info); // Cholesky decomposition
      switch (nconst) {
      case 0:
        dpotrs(uplo,&ngood,&ione,a,&ngood,coeffd,&ngood,&info);
        break;
      case 1:
        for (i=0;i<ngood;i++) a1b[i]=1.;
        dpotrs(uplo,&ngood,&ione,a,&ngood,a1b,&ngood,&info);
        bta1b=0.;
        for (i=0;i<ngood;i++) bta1b=bta1b+a1b[i];
        for (i=0;i<ngood;i++) a1r[i]=rh0[wgood[i]];
        dpotrs(uplo,&ngood,&ione,a,&ngood,a1r,&ngood,&info);
        bta1r=0.0;
        for (i=0;i<ngood;i++) bta1r=bta1r+a1r[i];
        c=(bta1r-1.)/bta1b;
        for (i=0;i<ngood;i++) coeffd[i]=a1r[i]-c*a1b[i];
        break;
      default:
// Set up array of delta times
        for (i=0;i<ngood;i++) {
          xc[i]=tsample[wgood[i]]-tint;
        }
// Set up constraint matrix
// First row is 1
        for (i=0;i<ngood;i++) b[i]=1.;
        for (i=0;i<ngood;i++) a1b[i]=1.;
        rhc[0]=1.0; // Constraint is 1
// Other rows are delta times ^ polynomial order
        for (iconst=1;iconst<nconst;iconst++) {
          for (i=0;i<ngood;i++) {
            b[iconst*ngood+i]=a1b[(iconst-1)*ngood+i]*xc[i];
            a1b[iconst*ngood+i]=a1b[(iconst-1)*ngood+i]*xc[i];
          }
          rhc[iconst]=0.0;
        }
// Solve system for constraint matrix 
        dpotrs(uplo,&ngood,&nconst,a,&ngood,a1b,&ngood,&info);
// bta1b is a matrix now, so use different variable;
        for (iconst=0;iconst<nconst;iconst++) {
          for (jconst=0;jconst<nconst;jconst++) {
            sum=0.0;
            for (i=0;i<ngood;i++) {
              sum=sum+b[iconst*ngood+i]*a1b[jconst*ngood+i];
            }
            bta1b1[iconst*nconst+jconst]=sum;
          }
        }
        dpotrf(uplo,&nconst,bta1b1,&nconst,&info); // Cholesky decomposition
        for (iconst=0;iconst<nconst;iconst++) {
          sum=0.0;
          for (i=0;i<ngood;i++) sum=sum+a1b[iconst*ngood+i]*rh0[wgood[i]];
          help[iconst]=sum-rhc[iconst];
        }
        dpotrs(uplo,&nconst,&ione,bta1b1,&nconst,help,&nconst,&info);
        for (i=0;i<ngood;i++) a1r[i]=rh0[wgood[i]];
        dpotrs(uplo,&ngood,&ione,a,&ngood,a1r,&ngood,&info);
        for (i=0;i<ngood;i++) {
          sum=0.0;
          for (iconst=0;iconst<nconst;iconst++) {
            sum=sum+help[iconst]*a1b[iconst*ngood+i];
          }
          coeffd[i]=a1r[i]-sum;
        }
        break;
      } // switch (nconst)
      for (i=0;i<ngood;i++) cp[wgood[i]]=(float)coeffd[i];
/*
sum=0.0;
for (i=0;i<ngood;i++) sum=sum+coeffd[i]*tsample[wgood[i]];
printf("%d %d %f %f %f\n",imask,ngood,tint,sum,sum-tint);
*/
    } // if ngood
  } // for imask
  t1=dsecnd()-t1;

  t2=dsecnd();
// Now do actual interpolation
// Could play include games here, but for the time being not
/*
// Inner loop over time
  for (j=0;j<ny;j++) {
    for (i=0;i<nx;i++) {
      ix=i+nlead*j; // Index into array
// First encode mask in single integer
      smask=0;
      for (isample=nsample-1;isample>=0;isample--) {
        if (masks[isample][ix]!=maskgood) smask=2*smask+1; else smask=2*smask;
      }
      sum=0.0f;
      if (smask==0) { // All are there
        for (isample=0;isample<nsample;isample++) sum=sum+coeffs[isample]*images[isample][ix];
      }
      else { // Some are missing.
// Could keep wgood and index instead
        cp=coeffs+nsample*smask;
        for (isample=0;isample<nsample;isample++) {
// This if statement is costly but needed to avoid accessing NaNs
          if (masks[isample][ix]==maskgood) {
            sum=sum+cp[isample]*images[isample][ix];
          }
        }
      } // if (smask)
      image_out[ix]=sum;
    } // for i
  } // for j
*/

// Alternate version with inner loop over x
  smasks=(int *)(MKL_malloc(nx*sizeof(int),malign));
  ngoods=(int *)(MKL_malloc(nx*sizeof(int),malign));
  sums=(float *)(MKL_malloc(nx*sizeof(float),malign));
  for (j=0;j<ny;j++) {
    for (i=0;i<nx;i++) smasks[i]=0;
    for (i=0;i<nx;i++) ngoods[i]=0;
    for (isample=nsample-1;isample>=0;isample--) {
      mp=masks[isample]+nlead*j;
      for (i=0;i<nx;i++) {
        if (mp[i]!=maskgood)
          smasks[i]=2*smasks[i]+1;
        else {
          smasks[i]=2*smasks[i];
          ngoods[i]=ngoods[i]+1;
        }
      }
    }
    for (i=0;i<nx;i++) sums[i]=0.0f;
    for (isample=0;isample<nsample;isample++) {
      mp=masks[isample]+nlead*j;
      cp=coeffs+isample;
      ip=images[isample]+nlead*j;
      for (i=0;i<nx;i++) {
        if (mp[i]==maskgood) {
          sums[i]=sums[i]+cp[nsample*smasks[i]]*ip[i];
        }
      } // for i
    } // for isample
    ip=image_out+nlead*j;
    for (i=0;i<nx;i++) 
      if (ngoods[i] >= nconst) // Enough points to do interpolation
        ip[i]=sums[i];
      else // No points present
        ip[i]=fillval;
//  memcpy(ip,sums,nx*sizeof(float));
  } // for j
  MKL_free(smasks);
  MKL_free(ngoods);
  MKL_free(sums);

/*
nblock=1024;
// Alternate version with inner loop over x and blocks
  smasks=(int *)(MKL_malloc(nblock*sizeof(int),malign));
  sums=(float *)(MKL_malloc(nblock*sizeof(float),malign));
  for (j=0;j<ny;j++) {
    for (ib=0;ib<nx;ib=ib+nblock) {
      for (i=0;i<minval(nblock,nx-ib);i++) smasks[i]=0;
      for (isample=nsample-1;isample>=0;isample--) {
        mp=masks[isample]+ib;
        for (i=0;i<minval(nblock,nx-ib);i++) {
          if (mp[i]!=maskgood)
            smasks[i]=2*smasks[i]+1;
          else
            smasks[i]=2*smasks[i];
        }
      }
      for (i=0;i<minval(nblock,nx-ib);i++) sums[i]=0.0f;
      for (isample=0;isample<nsample;isample++) {
        mp=masks[isample]+ib;
        cp=coeffs+isample;
        ip=images[isample]+nlead*j+ib;
        for (i=0;i<minval(nblock,nx-ib);i++) {
          if (mp[i]==maskgood) {
            sums[i]=sums[i]+cp[nsample*smasks[i]]*ip[i];
          }
        } // for i
      } // for isample
      ip=image_out+nlead*j+ib;
      for (i=0;i<minval(nblock,nx-ib);i++) ip[i]=sums[i];
    } // for ib
  } // for j
  MKL_free(smasks);
  MKL_free(sums);
*/

  t2=dsecnd()-t2;
  printf("%d %f %f %f\n",nsample,t0,t1,t2);

  MKL_free(acor);
  MKL_free(a0);
  MKL_free(a);
  MKL_free(rh0);
  MKL_free(coeffd);
  MKL_free(coeff);
  MKL_free(a1b);
  MKL_free(a1r);
  MKL_free(coeffs);
  MKL_free(wgood);
  free (filename);
  return 0;
}


int taverage(
  int nsample, // Number of input times
  double *tsample, // Input times
  double tint, // Target time
  int nconst, // Number of polynomial terms exactly reproduced
  float **images, // Pointer array to input images
  unsigned char **masks, // Pointer array to input masks. 0=good, 1= missing
  float *image_out, // Interpolated image
  int nx, // Number of points in dimension adjacent in memory
  int ny, // Number of points in dimension not adjacent in memory
  int nlead, // Leading dimension of arrays. nlead>=nx
  int method, // Interpolation method
  char **filenamep, // Pointer to name of file to read covariance from.
                   // Set to actual file used if method > 0.
  int avmethod, // averaging method
  int order, // Interpolation order
  double tspace, // Spacing of times to interpolate to
  int hwidth, // Window width in units of tspace. Total width is 2*hwidth+1
  double par1, // In units of tspace. Meaning depends on avmethod.
  double par2, // In units of tspace. Meaning depends on avmethod.
  float fillval // Value to use if not enough points present
)
{
  const unsigned char maskgood=0; // Value of good entries in masks
  int malign=32;
  int ncor=5000;
  //double dta=1.0; // Time between samples of input autocorrelation
  int i,j,k,isample,iconst,jconst;
  FILE *fileptr;
  double ti,dt,dtx,dta;
  double *acor; // Input autocorrelation
  double *a0,*a,*rh0,*coeffd;
  double *b,*a1b,bta1b,*a1r,bta1r,c,*rhc,*xc,*bta1b1,*help;
  int idt;
  int info;
  char uplo[] = "U";
  int ione = 1;
  float *coeff;
  int nmasks,imask,imask1,ngood;
  int smask;
  int *wgood;
  int *smasks;
  float *coeffs,*cp;
  float sum,*sums,*ip;
  double t0,t1,t2;
  unsigned char *mp;
  int ix;
  //int ib, nblock;
  char *filename;
  int ibad;

  if ((order < 1) || (order > 20)) {
// Code breaks at 31 or 32 due to the use of summed masks
// Also, memory usage and efficiency drops due to all masks calculated for >20
    printf("taverage: must have 1<=order<=20\n");
    return 1;
  }
  t0=dsecnd();
  switch (method) {
  case 0:
    filename=strdup(*filenamep);
    break;
  case 1:
    filename=strdup(DEFS_MKPATH("/data/acort_hr12.txt"));
    break;
  default:
    printf("Unknown method in init_fill.\n");
    return 1;
    break;
  }
  fileptr = fopen (filename,"r");
  if (fileptr==NULL) {
    printf("File not found in taverage.\n");
    return 1;
  }
  if (method != 0) {
    if (filenamep != NULL) {
      *filenamep=strdup(filename);
    }
  }

// Read autocorrelation
  acor=(double *)(MKL_malloc(ncor*sizeof(double),malign));
//fread ((char*)(acor),sizeof(double),ncor,fileptr);
  for (i=0;i<ncor;i++) fscanf(fileptr,"%lf",acor+i);
  fclose(fileptr);
  t0=dsecnd()-t0;

int width,*ix1,*ixi,*ihelp,iwidth,iorder;
double *tw,*tx1,*txi,*thelp,*weights,wsum;
int imin,imin1,minuse,maxuse,nmiss;
double tmin;
float psum;

// Now find times to do averaging over and interpolate to.
  width=2*hwidth+1;
  tw=(double *)(MKL_malloc(width*sizeof(double),malign));
  weights=(double *)(MKL_malloc(width*sizeof(double),malign));
  wsum=0.0;
  for (i=0;i<width;i++) {
    dt=i-hwidth+0.0;
    dta=fabs(dt);
    tw[i]=tint+tspace*dt;
    switch (avmethod) {
    case tavg_boxcar:
      if (dta <= par1) {
        weights[i]=1.0; 
      }
      else {
        weights[i]=0.0;
      }
      break;
    case tavg_cosine:
      if (dta <= (par1-par2)) {
        weights[i]=1.0;
      }
      else if (dta >= (par1+par2)) {
        weights[i]=0.0;
      }
      else {
        dtx=(dta-par1)/par2;
        weights[i]=0.5*(1-sin(M_PI*dtx/2));
      }
      break;
    case tavg_fourth:
      if (dta <= (par1-par2)) {
        weights[i]=1.0;
      }
      else if (dta >= (par1+par2)) {
        weights[i]=0.0;
      }
      else {
        dtx=(dta-par1)/par2;
        weights[i]=0.5+(-0.75+0.25*dtx*dtx)*dtx;
      }
      break;
    case tavg_hathaway:
      weights[i]=(exp(dt*dt/(2*par1*par1))-exp(par2*par2/(2*par1*par1)))*(1+(par2*par2-dt*dt)/(2*par1*par1));
      break;
    default:
      printf("Unimplemented method in taverage\n");
      return 1;
    }
    printf("%d %f %f %f\n",i,dt,dta,weights[i]);
    wsum=wsum+weights[i];
  } // for i
// Normalize weigths
  for (i=0;i<width;i++) weights[i]=weights[i]/wsum;
  //for (i=0;i<width;i++) printf("%f\n",weights[i]);

// Get list of (closest) input times to use for each target time.
  ix1=(int *)(MKL_malloc(order*width*sizeof(int),malign));
  tx1=(double *)(MKL_malloc(order*width*sizeof(double),malign));
  thelp=(double *)(MKL_malloc(nsample*sizeof(double),malign));
  ihelp=(int *)(MKL_malloc(nsample*sizeof(int),malign));

  minuse=nsample; // Minimum index actually used in input
  maxuse=-1; // Maximum index actually used in input
  for (i=0;i<width;i++) { // Loop over target times and find closest order times.
    for (k=0;k<nsample;k++) {
      ihelp[k]=k;
      thelp[k]=tsample[k];
    }
    for (j=0;j<order;j++) { // Find closest order points.
      tmin=1e6;
      for (k=0;k<nsample-j;k++) { // Loop over remaining points and find closest one.
        if (fabs(thelp[k]-tw[i]) < tmin) {
          imin=k;
          tmin=fabs(thelp[k]-tw[i]);
        }
      }
      imin1=ihelp[imin]; // Index into original table
      minuse=minval(minuse,imin1);
      maxuse=maxval(maxuse,imin1);
      ix1[i*order+j]=imin1;
      tx1[i*order+j]=tsample[imin1];
      for (k=imin;k<nsample-j-1;k++) { // Remove smallest element.
        ihelp[k]=ihelp[k+1];
        thelp[k]=thelp[k+1];
      }
    }
/*
    printf("%d %f",i,tw[i]);
    for (j=0;j<order;j++) printf(" %d",ix1[i*order+j]); printf("\n");
    for (j=0;j<order;j++) printf(" %f",tx1[i*order+j]); printf("\n");
*/
  }

  t1=dsecnd();
  a0=(double *)(MKL_malloc(nsample*nsample*sizeof(double),malign));
  a=(double *)(MKL_malloc(nsample*nsample*sizeof(double),malign));
  rh0=(double *)(MKL_malloc(nsample*sizeof(double),malign));
  coeffd=(double *)(MKL_malloc(nsample*sizeof(double),malign));
  coeff=(float *)(MKL_malloc(nsample*sizeof(float),malign));
  b=(double *)(MKL_malloc(nconst*nsample*sizeof(double),malign));
  a1b=(double *)(MKL_malloc(nconst*nsample*sizeof(double),malign));
  help=(double *)(MKL_malloc(nconst*sizeof(double),malign));
  rhc=(double *)(MKL_malloc(nconst*nsample*sizeof(double),malign));
  bta1b1=(double *)(MKL_malloc(nconst*nconst*sizeof(double),malign));
  xc=(double *)(MKL_malloc(nsample*sizeof(double),malign));
  a1r=(double *)(MKL_malloc(nconst*nsample*sizeof(double),malign));

  nmasks=1;
  for (i=0;i<order;i++) nmasks=2*nmasks;

//Coeffs will now hold weights for all target and input times.
  coeffs=(float *)(MKL_malloc(width*order*nmasks*sizeof(float),malign));
  wgood=(int *)(MKL_malloc(width*sizeof(int),malign));

  for (iwidth=0;iwidth<width;iwidth++) { // Loop over target times
    txi=tx1+iwidth*order; // Times for this target
    for (i=0;i<order;i++) { // Loop over input times for that target
      ti=txi[i];
      for (j=0;j<order;j++) {
        dt=fabs(txi[j]-ti);
        idt=(float)floor(dt);
        dtx=dt-idt;
        a0[i*order+j]=(1.0-dtx)*acor[idt]+dtx*acor[idt+1];
      }
      dt=fabs(tw[iwidth]-ti);
      idt=(float)floor(dt);
      dtx=dt-idt;
      rh0[i]=(1.0-dtx)*acor[idt]+dtx*acor[idt+1];
    }

// For the time being calculate weights for all possible masks
// Could be done implicitly if order becomes large
// For 4096x4096 images interpolation dominates to around order=20
    for (imask=0;imask<nmasks;imask++) {
//printf("%d %d\n",isample,imask);
      cp=coeffs+order*imask+iwidth*order*nmasks;
      imask1=imask;
      ngood=0;
// Find good and bad given mask
      for (i=0;i<order;i++) {
        if ((imask1&1)==0) { // Good data
          wgood[ngood]=i;
          ngood=ngood+1;
        }
        imask1=imask1/2;
        cp[i]=0.0f; // Zero all coefficients
      } // for i
      if (ngood>0) {
// Pick out relevant parts of a0 and rh0
        for (i=0;i<ngood;i++) {
          for (j=0;j<ngood;j++) a[ngood*i+j]=a0[order*wgood[i]+wgood[j]];
          coeffd[i]=rh0[wgood[i]];
        }
        dpotrf(uplo,&ngood,a,&ngood,&info); // Cholesky decomposition
        switch (nconst) {
        case 0:
          dpotrs(uplo,&ngood,&ione,a,&ngood,coeffd,&ngood,&info);
          break;
        case 1:
          for (i=0;i<ngood;i++) a1b[i]=1.;
          dpotrs(uplo,&ngood,&ione,a,&ngood,a1b,&ngood,&info);
          bta1b=0.;
          for (i=0;i<ngood;i++) bta1b=bta1b+a1b[i];
          for (i=0;i<ngood;i++) a1r[i]=rh0[wgood[i]];
          dpotrs(uplo,&ngood,&ione,a,&ngood,a1r,&ngood,&info);
          bta1r=0.0;
          for (i=0;i<ngood;i++) bta1r=bta1r+a1r[i];
          c=(bta1r-1.)/bta1b;
          for (i=0;i<ngood;i++) coeffd[i]=a1r[i]-c*a1b[i];
          break;
        default:
// Set up array of delta times
          for (i=0;i<ngood;i++) {
            xc[i]=tx1[iwidth*order+wgood[i]]-tw[iwidth];
          }
// Set up constraint matrix
// First row is 1
          for (i=0;i<ngood;i++) b[i]=1.;
          for (i=0;i<ngood;i++) a1b[i]=1.;
          rhc[0]=1.0; // Constraint is 1
// Other rows are delta times ^ polynomial order
          for (iconst=1;iconst<nconst;iconst++) {
            for (i=0;i<ngood;i++) {
              b[iconst*ngood+i]=a1b[(iconst-1)*ngood+i]*xc[i];
              a1b[iconst*ngood+i]=a1b[(iconst-1)*ngood+i]*xc[i];
            }
            rhc[iconst]=0.0;
          }
// Solve system for constraint matrix 
          dpotrs(uplo,&ngood,&nconst,a,&ngood,a1b,&ngood,&info);
// bta1b is a matrix now, so use different variable;
          for (iconst=0;iconst<nconst;iconst++) {
            for (jconst=0;jconst<nconst;jconst++) {
              sum=0.0;
              for (i=0;i<ngood;i++) {
                sum=sum+b[iconst*ngood+i]*a1b[jconst*ngood+i];
              }
              bta1b1[iconst*nconst+jconst]=sum;
            }
          }
          dpotrf(uplo,&nconst,bta1b1,&nconst,&info); // Cholesky decomposition
          for (iconst=0;iconst<nconst;iconst++) {
            sum=0.0;
            for (i=0;i<ngood;i++) sum=sum+a1b[iconst*ngood+i]*rh0[wgood[i]];
            help[iconst]=sum-rhc[iconst];
          }
          dpotrs(uplo,&nconst,&ione,bta1b1,&nconst,help,&nconst,&info);
          for (i=0;i<ngood;i++) a1r[i]=rh0[wgood[i]];
          dpotrs(uplo,&ngood,&ione,a,&ngood,a1r,&ngood,&info);
          for (i=0;i<ngood;i++) {
            sum=0.0;
            for (iconst=0;iconst<nconst;iconst++) {
              sum=sum+help[iconst]*a1b[iconst*ngood+i];
            }
            coeffd[i]=a1r[i]-sum;
          }
          break;
        } // switch
        for (i=0;i<ngood;i++) cp[wgood[i]]=(float)coeffd[i];
/*
sum=0.0;
for (i=0;i<ngood;i++) sum=sum+coeffd[i]*tx1[iwidth*order+wgood[i]];
printf("%d %d %f %f\n",imask,ngood,sum,tw[iwidth]-sum);
*/
      } // if ngood
    } // for imask
/*
printf("%d",iwidth);
for (iorder=0;iorder<order;iorder++) printf(" %f",coeffs[iwidth*order*nmasks+iorder]);
printf("\n");
*/
  } // for iwidth
  t1=dsecnd()-t1;

  t2=dsecnd();

// Calculate coefficients for case of no missing data.
  float *coeff0;
  coeff0=(float *)(MKL_malloc(nsample*sizeof(float),malign));
  for (isample=0;isample<nsample;isample++) coeff0[isample]=0.0f;
  for (iwidth=0;iwidth<width;iwidth++) {
    sum=0.0;
    txi=tx1+iwidth*order; // Times for this target
    for (iorder=0;iorder<order;iorder++) {
      sum=sum+txi[iorder]*coeffs[iwidth*order*nmasks+iorder];
      isample=ix1[iwidth*order+iorder];
      coeff0[isample]=coeff0[isample]+coeffs[iwidth*order*nmasks+iorder]*weights[iwidth];
    }
    printf("%d %f %f %f\n",iwidth,tw[iwidth],sum,sum-tw[iwidth]);
  }

  sum=0.0;
  for (isample=0;isample<nsample;isample++) {
    printf("%d %f %f\n",isample,tsample[isample],coeff0[isample]);
    sum=sum+tsample[isample]*coeff0[isample];
  }
  printf("%f\n",sum);

// Now do actual interpolation
// Inner loop over time
  for (j=0;j<ny;j++) {
    for (i=0;i<nx;i++) {
      ix=i+nlead*j; // Index into array
      nmiss=0;
      for (isample=0;isample<nsample;isample++) {
        if (masks[isample][ix]!=maskgood) nmiss=nmiss+1;
      }
      sum=0.0f;
      ibad=0;
      if (nmiss==0) { // All are there
        for (isample=0;isample<nsample;isample++) sum=sum+coeff0[isample]*images[isample][ix];
      }
      else { // Some are missing.
// For the time being resort to brute force calculation.
        for (iwidth=0;iwidth<width;iwidth++) {
// First encode mask in single integer
          smask=0;
          ngood=0;
          for (iorder=order-1;iorder>=0;iorder--) { // Must go backwards because of mask definition.
            isample=ix1[iwidth*order+iorder];
            if (masks[isample][ix]!=maskgood)
              smask=2*smask+1;
            else {
              smask=2*smask;
              ngood=ngood+1;
            }
          }
          if (ngood < nconst) ibad=1; // Got a bad point!
          cp=coeffs+order*smask+iwidth*order*nmasks;
          psum=0.0f; // Partial sum
          if (smask==0) { // All are there for this target time
            for (iorder=0;iorder<order;iorder++) {
              isample=ix1[iwidth*order+iorder];
              psum=psum+cp[iorder]*images[isample][ix];
            }
          }
          else { // Some are missing.
            for (iorder=0;iorder<order;iorder++) {
              isample=ix1[iwidth*order+iorder];
              if (masks[isample][ix]==maskgood) {
                psum=psum+cp[iorder]*images[isample][ix];
              }
            }
          } // if (smask)
          sum=sum+psum*weights[iwidth];
        } // for iwidth
      } // if (nmiss)
      if (ibad == 0) {
        image_out[ix]=sum;
      }
      else {
        image_out[ix]=fillval; // Use fillval if not enough points
      }
    } // for i
  } // for j

  t2=dsecnd()-t2;
  //printf("%d %f %f %f\n",nsample,t0,t1,t2);

  MKL_free(acor);
  MKL_free(a0);
  MKL_free(a);
  MKL_free(rh0);
  MKL_free(coeffd);
  MKL_free(coeff);
  MKL_free(a1b);
  MKL_free(a1r);
  MKL_free(coeffs);
  MKL_free(wgood);
  free (filename);
  return 0;
}

char *tinterpolate_version() // Returns CVS version of tinterpolate.c
{
  return strdup("$Id: tinterpolate.c,v 1.8 2010/09/08 00:04:45 schou Exp $");
}


