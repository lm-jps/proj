#include <stdio.h>

#include <stdlib.h>

#include <time.h>

#include <string.h>

#include <math.h>

#include <mkl_blas.h>

#include <mkl_service.h>

#include <mkl_lapack.h>

#include <mkl_vml_functions.h>

#include <omp.h>

#include "fresize.h"

#define minval(x,y) (((x) < (y)) ? (x) : (y))



int main(int argc, const char **argv)

{

  const int nsubin0=100,nsubin=nsubin0+1;

  int nx,ny;

  float *image_in,*image_out;

  int i,j,ib,jb,ib1,jb1,i1,j1,k;

  double t0,t1;

  int malign=32;

  FILE *fileptr;

  int nlead;

  struct fresize_struct fresizes;

  int nt;

  int nbin,nxout,nyout;



  nlead=4128;

  nx=4096;

  ny=4096;

  image_in=(float *)(MKL_malloc(nlead*ny*sizeof(float),malign));

  image_out=(float *)(MKL_malloc(nlead*ny*sizeof(float),malign));



  t0=dsecnd();

  for (j=0; j<ny; j++) {

    for (i=0; i<nx; i++) {

//    image_in[j*nlead+i]=sin(2*M_PI*i/30.1f)+sin(2*M_PI*j/40.3f);

      image_in[j*nlead+i]=1.0f;

    }

    for (i=nx;i<nlead;i++) image_in[j*nlead+i]=1e20;

  }

  printf("Nlead: %d.\n",nlead);

  printf("Malign: %d.\n",malign);

  t0=dsecnd()-t0;

  printf("Time to set up = %f seconds.\n",t0);

 

nt=1;{

//for (nt=1;nt<=8;nt=nt+1) {

omp_set_num_threads(nt);



for (nbin=1;nbin<=1;nbin=2*nbin) {

  nxout=nx/nbin;

  nyout=ny/nbin;



//t0=dsecnd();

//init_fresize_sample(&fresizes,nbin);

//init_fresize_bin(&fresizes,nbin);

//init_fresize_boxcar(&fresizes,10,nbin);

  init_fresize_gaussian(&fresizes,20,100,nbin);

//t0=dsecnd()-t0;

//printf("Time to initialize = %f seconds.\n",t0);



  fresize(&fresizes,image_in, image_out, nx,ny,nlead,nxout,nyout,nxout,0,0,0.0f);

  t1=dsecnd();

  fresize(&fresizes,image_in, image_out, nx,ny,nlead,nxout,nyout,nxout,2,5,0.0f);



  t1=dsecnd()-t1;

  printf("%d %d %f\n",nt,nbin,t1);



  free_fresize(&fresizes);

}

}





  fileptr = fopen ("/tmp21/schou/imx.bin", "w");

  for (i=0;i<4096;i++) fwrite ((char*)(image_in+i*nlead),sizeof(float),4096,fileptr);

  for (i=0;i<nyout;i++) fwrite ((char*)(image_out+i*nxout),sizeof(float),nxout,fileptr);

  fclose(fileptr);

 

  MKL_free(image_in);

  MKL_free(image_out);

  MKL_FreeBuffers();

  

  return 0;

}





