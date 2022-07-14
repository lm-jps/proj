#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <mkl_blas.h>
#include <mkl_service.h>
#include <mkl_lapack.h>
#include <mkl_vml_functions.h>
#include "cholesky_down.h"
#include "gapfill.h"
#define minval(x,y) (((x) < (y)) ? (x) : (y))
#define maxval(x,y) (((x) < (y)) ? (y) : (x))

int init_fill(
int method,
double pnoise,
int order,
int targetx,
int targety,
struct fill_struct *pars
)
{
  double *acor;
  int nx=80;
  int malign=32;
  double t0,t1;
  FILE *fileptr;
  int i,j;
  int *offx,*offy,ox,oy;
  int nsample,inoise;
  double regul,pmin;
  double *a0,*rh0,*a0t,*rh0t;
  int info;
  char uplo[] = "U";

  pars->order=order;
  
  t0=dsecnd();
//fileptr = fopen ("/tmp20/schou/acor1_80x1_double.bin", "r");
  switch (method) {
  case 8:
    fileptr = fopen ("/tmp20/schou/acor1_08_80.bin", "r");
    break;
  case 9:
    fileptr = fopen ("/tmp20/schou/acor1_09_80.bin", "r");
    break;
  case 10:
    fileptr = fopen ("/tmp20/schou/acor1_10_80.bin", "r");
    break;
  case 11:
    fileptr = fopen ("/tmp20/schou/acor1_11_80.bin", "r");
    break;
  default:
    printf("Unknown method in init_fill.\n");
    return 1;
    break;
  }
  if (fileptr==NULL) {
    printf("File not found in init_fill.\n");
    return 1;
  }
  acor=(double *)(MKL_malloc(nx*nx*sizeof(double),malign));
  fread ((char*)acor,sizeof(double),nx*nx,fileptr);
  fclose(fileptr);
  t0=dsecnd()-t0;
  printf("x0 %f\n",t0);

  regul=pnoise*pnoise;
  inoise=1;
  pmin=regul*inoise;

  nsample=order*order;
  t1=dsecnd();
  offx=(int *)(MKL_malloc(nsample*sizeof(int),malign));
  offy=(int *)(MKL_malloc(nsample*sizeof(int),malign));
  for (i=0;i<nsample;i++) {
    offx[i]=(i%order)-targetx;
    offy[i]=(i/order)-targety;
  }
// Set up complete matrix
  pars->a0=(double *)(MKL_malloc(nsample*nsample*sizeof(double),malign));
  pars->a00=(double *)(MKL_malloc(nsample*nsample*sizeof(double),malign));
  pars->a0t=(double *)(MKL_malloc(nsample*nsample*sizeof(double),malign));
  a0=pars->a0;
  a0t=pars->a0t;
  for (i=0;i<nsample;i++) {
    for (j=0;j<nsample;j++) {
      ox=((offx[j]-offx[i])+nx) % nx;
      oy=((offy[j]-offy[i])+nx) % nx;
      a0[i+nsample*j]=acor[ox+nx*oy];
      a0t[i+nsample*j]=acor[ox+nx*oy];
    }
    a0[i+nsample*i]=a0[i+nsample*i]+regul;
  }
// save copy of matrix
  memcpy(pars->a00,a0,nsample*nsample*sizeof(double));
  pars->acort00=acor[0];
  pars->rh0=(double *)(MKL_malloc(nsample*sizeof(double),malign));
  pars->rh0t=(double *)(MKL_malloc(nsample*sizeof(double),malign));
  rh0=pars->rh0;
  rh0t=pars->rh0t;
  for (i=0;i<nsample;i++) {
    rh0[i]=acor[(-offx[i]+nx) % nx+nx*((-offy[i]+nx) % nx)];
    rh0t[i]=acor[(-offx[i]+nx) % nx+nx*((-offy[i]+nx) % nx)];
    if ((offx[i]==0) && (offy[i]==0)) rh0[i]=rh0[i]+pmin;
  }
// Decompose complete matrix
  dpotrf(uplo,&nsample,pars->a0,&nsample,&info); // Cholesky decomposition

  MKL_free(acor);
  MKL_free(offx);
  MKL_free(offy);

// Variables used in further calculations. Avoids multiple mallocs.
  pars->a=(double *)(MKL_malloc(nsample*nsample*sizeof(double),malign));
  pars->rh=(double *)(MKL_malloc(nsample*sizeof(double),malign));
  pars->a1b=(double *)(MKL_malloc(nsample*sizeof(double),malign));
  pars->a1r=(double *)(MKL_malloc(nsample*sizeof(double),malign));
  pars->wgood=(int *)(MKL_malloc(nsample*sizeof(int),malign));
  pars->wbad=(int *)(MKL_malloc(nsample*sizeof(int),malign));

  pars->hashmod=8191; // Hash table size
  pars->hashcount=(int *)(MKL_malloc(pars->hashmod*sizeof(int),malign));
  pars->hashtable=(struct fill_hash_struct **)(MKL_malloc(pars->hashmod*sizeof(struct fill_hash_struct *),malign));
  for (i=0;i<pars->hashmod;i++) {
    pars->hashcount[i]=0;
// Initialize to null pointers
    pars->hashtable[i]=(struct fill_hash_struct *)0;
  }
  pars->ndiff=0;
  pars->ncollision=0;

  t1=dsecnd()-t1;
  printf("x1 %f\n",t1);

  return 0;
}

int free_fill(
struct fill_struct *pars
)
{
  int i;
  struct fill_hash_struct *ptr,*oldptr;
  
  printf("Hashmod: %d\n",pars->hashmod);
  printf("Ndiff: %d\n",pars->ndiff);
  printf("Ncollision: %d\n",pars->ncollision);
/*
  for (i=0;i<pars->hashmod;i++) {
    printf("%d %d\n",i,pars->hashcount[i]);
  }
  for (i=0;i<pars->hashmod;i++) {
    ptr=(pars->hashtable[i]);
    while (ptr != NULL) {
      printf("%d %d %d\n",i,ptr->nbad,ptr->nhit);
      ptr=ptr->next;
    }
  }
*/
  MKL_free(pars->a0);
  MKL_free(pars->a00);
  MKL_free(pars->a0t);
  MKL_free(pars->rh0);
  MKL_free(pars->rh0t);
  MKL_free(pars->a);
  MKL_free(pars->rh);
  MKL_free(pars->a1b);
  MKL_free(pars->a1r);
  MKL_free(pars->wgood);
  MKL_free(pars->wbad);

// Now free hash table and contents

  MKL_free(pars->hashcount);
  for (i=0;i<pars->hashmod;i++) {
    ptr=(pars->hashtable[i]);
    while (ptr != NULL) {
      oldptr=ptr;
      MKL_free(ptr->wbad);
      MKL_free(ptr->coeff);
      ptr=ptr->next;
      MKL_free(oldptr);
    }
  }
  MKL_free(pars->hashtable);

  return 0; // Currently no error conditions

}

double fill_point(
struct fill_struct *pars,
int *mask,
double *vals,
float *cnorm,
float *ierror
)
{
  double t2,t3,t4;
  int ncopy,nsample;
  int ngood,nbad;
  int *wgood,*wbad;
  int i,info,j;
  char uplo[] = "U";
  int ione = 1;
  char transn[] = "n";
  double done=1.0;
  double dzero=0.0;
  int malign=32;
  double *a0,*a,*rh0,*rh,*a1b,bta1b,*a1r,bta1r,c,*coeff,*a0t,*a00;
  float fillval;
  int hash,ident,found;
  struct fill_hash_struct *ptr;
  double ierror2; // interpolation error squared

  nsample=(pars->order)*(pars->order);
  
// Make work copies
  a0=pars->a0;
  a00=pars->a00;
  a0t=pars->a0t;
  rh0=pars->rh0;
  a=pars->a;
  rh=pars->rh;
  a1b=pars->a1b;
  a1r=pars->a1r;
  wgood=pars->wgood;
  wbad=pars->wbad;

  ngood=0;
  nbad=0;
  hash=0;
  for (i=0;i<nsample;i++) {
    if (mask[i]==0) {
      wgood[ngood]=i;
      ngood=ngood+1;
    }
    if (mask[i]!=0) {
      hash=hash+i;
      wbad[nbad]=i;
      nbad=nbad+1;
    }
  }
  hash=hash%pars->hashmod;

  found=0;
  ptr=pars->hashtable[hash];
  while (ptr) { // Not at end of table
    if (nbad==ptr->nbad) { // Number of elements is correct, so there is hope
      ident=1;
      for (i=0;i<nbad;i++) {
        if (wbad[i] != ptr->wbad[i]) {
          ident=0;
          pars->ncollision++;
          break; // Out of for loop. Known not to match.
        }
      }
      if (ident != 0) { // Success
        coeff=ptr->coeff; // Point to existing coeficients
        ptr->nhit++;
        found=1;
        break; // Out of while loop. No need to test the rest.
      } // if
    } // if
    ptr=ptr->next;
  } // while

  if (found == 0) { // Got to calculate new coefficients
// malloc new entry
//for (i=0;i<nsample;i++) printf(" %d",mask[i]);
//printf("\n");
    ptr=(struct fill_hash_struct *)(MKL_malloc(sizeof(struct fill_hash_struct),malign));
// Stick at beginning of linked list
    ptr->next=pars->hashtable[hash];
    pars->hashtable[hash]=ptr;
// Set other variables
    ptr->nbad=nbad;
    ptr->nhit=1;
    ptr->wbad=(int *)(MKL_malloc(nbad*sizeof(int),malign));
    for (i=0;i<nbad;i++) ptr->wbad[i]=wbad[i];
    ptr->coeff=(double *)(MKL_malloc(ngood*sizeof(double),malign));
    pars->hashcount[hash]++; // Increment count for this hash
    pars->ndiff++; // Increment count for total
    coeff=ptr->coeff; // Point to new coefficients

// Calculate new coefficients

// Constants in if statement based on timings
    if (nbad < (1.50*sqrt(nsample)-6)) {
// Get new decomposition by downsampling old one
// Copy decomposed matrix.
// This dominates dcholesky_down when only removing a few points.
t2=dsecnd();
// Only need to copy relevant part of matrix
      for (i=0;i<nsample;i++) for (j=0;j<=i;j++) a[i*nsample+j]=a0[i*nsample+j];
t2=dsecnd()-t2;

// Get decomposition of reduced matrix
t3=dsecnd();
      dcholesky_down_mask(uplo, &nsample, a, &nsample, mask, &info);
t3=dsecnd()-t3;
    }
    else {
// Calculate brute force from original matrix
// Pick good points
t2=dsecnd();
      for (i=0;i<ngood;i++) for (j=0;j<=i;j++) a[i*nsample+j]=a00[wgood[i]*nsample+wgood[j]];
t2=dsecnd()-t2;

// Get decomposition of reduced matrix
t3=dsecnd();
      dpotrf(uplo,&ngood,a,&nsample,&info); // Cholesky decomposition
t3=dsecnd()-t3;
}

t4=dsecnd();
    for (i=0;i<ngood;i++) a1b[i]=1.;
    dpotrs(uplo,&ngood,&ione,a,&nsample,a1b,&nsample,&info);
    bta1b=0.;
    for (i=0;i<ngood;i++) bta1b=bta1b+a1b[i];
    for (i=0;i<ngood;i++) a1r[i]=rh0[wgood[i]];
    dpotrs(uplo,&ngood,&ione,a,&nsample,a1r,&nsample,&info);
    bta1r=0.0;
    for (i=0;i<ngood;i++) bta1r=bta1r+a1r[i];
    c=(bta1r-1.)/bta1b;
    for (i=0;i<ngood;i++) coeff[i]=a1r[i]-c*a1b[i];
    ptr->cnorm=dnrm2(&ngood,coeff,&ione);

// Now get interpolation error
    for (i=0;i<ngood;i++) rh[i]=rh0[wgood[i]];
    ierror2=pars->acort00-2*ddot(&ngood,coeff,&ione,rh,&ione);
// Reuse a to hold subset of a0t
    for (j=0;j<ngood;j++) {
      for (i=0;i<ngood;i++) {
        a[ngood*j+i]=a0t[nsample*wgood[j]+wgood[i]];
      }
    }
// Use rh to hold a0t#coeff
    dgemv(transn,&ngood,&ngood,&done,a,&ngood,coeff,&ione,&dzero,rh,&ione);
    ierror2=ierror2+ddot(&ngood,coeff,&ione,rh,&ione);
    ptr->ierror=sqrt(ierror2/pars->acort00);
t4=dsecnd()-t4;
//printf("%d %f %f %f",ngood,t2,t3,t4);
//for (i=0;i<nsample;i++) printf(" %d",mask[i]);
//if (nbad == 2) printf(" %d %d",wbad[0],wbad[1]);
//printf("\n");

// Could store coefficients in the other 7 rotations of the mask.
// Beware that the time for sholesky_down depends strongly on which
// one is calculated. See timep.pro for calculation (t2s is a good predictor).
  } // End of calculating new coefficients
  fillval=ddot(&ngood,coeff,&ione,vals,&ione);
  *cnorm=ptr->cnorm;
  *ierror=ptr->ierror;
  return fillval;
}

int fgap_fill(
  struct fill_struct *pars,
  float *im,
  int nx,
  int ny,
  int nlead,
  unsigned char *mask,
  float *cnorm,
  float *ierror
)
{
  const unsigned char maskgood=0; // Value of good entries in masks
  const unsigned char maskfill=1; // Value of entries to be filled
  const unsigned char masknofill=2; // Value of entries not to be filled
  int malign=32;
  const int ione = 1;
  int i,j,i1,j1,i0,j0;
  int order,order2;
  int ngood,nfill;
  double *im1;
  int *mask1;
  float cnorm1,ierror1;
  
  order=pars->order;
  
  im1=(double *)(MKL_malloc(order*order*sizeof(double),malign));
  mask1=(int *)(MKL_malloc(order*order*sizeof(int),malign));
  
  order2=order/2;
  nfill=0;
  for (j=0; j<ny; j++) {
    j0=j-order2;
    for (i=0; i<nx; i++) {
      if (mask[j*nlead+i]==maskgood) {
        cnorm1=1.0f; // Random noise is 1
        ierror1=0.0f; // Interpolation error is 0
      }
      if (mask[j*nlead+i]==maskfill) { // Fill if mask=1. Skip for mask>1.
        i0=i-order2;
        ngood=0;
        for (j1=j0;j1<(j0+order);j1++) {
          for (i1=i0;i1<(i0+order);i1++) {
            if ((i1>=0) && (i1<nx) && (j1>=0) && (j1<ny)) {
              mask1[i1-i0+(j1-j0)*order]=mask[i1+j1*nlead];
              if (mask[j1*nlead+i1]==maskgood) {
                im1[ngood]=im[i1+j1*nlead];
                ngood=ngood+1;
              }
            } // Inside original image
            else { // Outside original image
              mask1[i1-i0+(j1-j0)*order]=1; // Point is missing
            }
          }
        }
        if (ngood != 0) { // Have at least one good point
          im[j*nlead+i]=fill_point(pars,mask1,im1,&cnorm1,&ierror1);
        } 
        else { // Otherwise quietly leave unchanged, except set errors to 0.0
          cnorm1=0.0f;
          ierror1=0.0f;
        }
        nfill=nfill+1;
      } // end if mask==1
      if (mask[j*nlead+i]==masknofill) { // Don't fill.
        cnorm1=0.0f;
        ierror1=0.0f;
      }
      cnorm[j*nlead+i]=cnorm1;
      ierror[j*nlead+i]=ierror1;
    }
  }
  printf("%d\n",nfill);
  MKL_free(im1);
  MKL_free(mask1);
  
  return 0; // Current no error conditions
}
