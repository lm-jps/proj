#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <mkl_blas.h>
#include <mkl_service.h>
#include <mkl_lapack.h>
#include <mkl_vml_functions.h>
#include "polcal.h"
#define minval(x,y) (((x) < (y)) ? (x) : (y))
#define maxval(x,y) (((x) > (y)) ? (x) : (y))

static void printm(double *m)
{
  int i;

  for (i=0;i<16;i++) {
    printf(" %f",m[i]);
    if ((i%4)==3) printf("\n");
  }
}

static void mat4mm(
  double *a,
  double *b,
  double *c
)
{
  int i,j,k;
  double sum;

  for (i=0;i<4;i++) {
    for (j=0;j<4;j++) {
      sum=0.0;
      for (k=0;k<4;k++) sum=sum+a[i+4*k]*b[k+4*j];
      c[i+4*j]=sum;
    }
  }
}

static void polarizer(
  double phi,
  double *m
)
{
  double c,s;

  c=cos(2*phi);
  s=sin(2*phi);

  m[0]=0.5;
  m[1]=0.5*c;
  m[2]=0.5*s;
  m[3]=0.0;
  m[4]=0.5*c;
  m[5]=0.5*c*c;
  m[6]=0.5*s*c;
  m[7]=0.0;
  m[8]=0.5*s;
  m[9]=0.5*s*c;
  m[10]=0.5*s*s;
  m[11]=0.0;
  m[12]=0.0;
  m[13]=0.0;
  m[14]=0.0;
  m[15]=0.0;

}

static void retarder(
  double delta, // Retardance in waves
  double phi, // Rotation angle in radians
  double *m // 4x4 matrix witt result
)
{
  double cd,sd,c,s;
  int i;
  double pi=M_PI;

  for (i=1;i<4;i++) m[i]=0.0;
  for (i=1;i<4;i++) m[4*i]=0.0;

  m[0]=1.0;
  cd=cos(2*pi*delta);
  sd=sin(2*pi*delta);
  c=cos(2*phi);
  s=sin(2*phi);
  m[5]=c*c+s*s*cd;
  m[6]=c*s*(1-cd);
  m[7]=s*sd;
  m[9]=c*s*(1-cd);
  m[10]=s*s+c*c*cd;
  m[11]=-c*sd;
  m[13]=-s*sd;
  m[14]=c*sd;
  m[15]=cd;
  
}

int init_polcal(
  struct polcal_struct *pars,
  int method
)
{
  int nin=32;
  int malign=32;
  int i,j;
  double d;
  FILE *fileptr;

  if (method!=1) {
    printf("Unimplemented method in init_polcal\n");
    return 1;
  }
  pars->method=method;
  pars->nin=nin;
  pars->xin=(double *)(MKL_malloc(nin*nin*sizeof(double),malign));
  pars->yin=(double *)(MKL_malloc(nin*nin*sizeof(double),malign));
  pars->fqq_0=(double *)(MKL_malloc(nin*nin*sizeof(double),malign));
  pars->fqq_1=(double *)(MKL_malloc(nin*nin*sizeof(double),malign));
  pars->fqq_2=(double *)(MKL_malloc(nin*nin*sizeof(double),malign));
  pars->fqu_0=(double *)(MKL_malloc(nin*nin*sizeof(double),malign));
  pars->fqv_0=(double *)(MKL_malloc(nin*nin*sizeof(double),malign));
  pars->fuu_0=(double *)(MKL_malloc(nin*nin*sizeof(double),malign));
  pars->fuu_1=(double *)(MKL_malloc(nin*nin*sizeof(double),malign));
  pars->fuu_2=(double *)(MKL_malloc(nin*nin*sizeof(double),malign));
  pars->fuv_0=(double *)(MKL_malloc(nin*nin*sizeof(double),malign));
  pars->fvv_0=(double *)(MKL_malloc(nin*nin*sizeof(double),malign));
  pars->fvv_1=(double *)(MKL_malloc(nin*nin*sizeof(double),malign));
  pars->fvv_2=(double *)(MKL_malloc(nin*nin*sizeof(double),malign));
  pars->ret1_0=(double *)(MKL_malloc(nin*nin*sizeof(double),malign));
  pars->ret1_1=(double *)(MKL_malloc(nin*nin*sizeof(double),malign));
  pars->ret2_0=(double *)(MKL_malloc(nin*nin*sizeof(double),malign));
  pars->ret2_1=(double *)(MKL_malloc(nin*nin*sizeof(double),malign));
  pars->ret3_0=(double *)(MKL_malloc(nin*nin*sizeof(double),malign));
  pars->ret3_1=(double *)(MKL_malloc(nin*nin*sizeof(double),malign));
  pars->phi1_0=(double *)(MKL_malloc(nin*nin*sizeof(double),malign));
  pars->phi2_0=(double *)(MKL_malloc(nin*nin*sizeof(double),malign));
  pars->phi3_0=(double *)(MKL_malloc(nin*nin*sizeof(double),malign));

  
  fileptr = fopen ("/scr/schou/pars_090416/fit.bin", "r");
  fread (&pars->tsela,sizeof(double),1,fileptr);
  fread (&pars->tfronta,sizeof(double),1,fileptr);
  fread (pars->fqq_0,sizeof(double),nin*nin,fileptr);
  fread (pars->fqq_1,sizeof(double),nin*nin,fileptr);
  fread (pars->fqq_2,sizeof(double),nin*nin,fileptr);
  fread (pars->fqu_0,sizeof(double),nin*nin,fileptr);
  fread (pars->fqv_0,sizeof(double),nin*nin,fileptr);
  fread (pars->fuu_0,sizeof(double),nin*nin,fileptr);
  fread (pars->fuu_1,sizeof(double),nin*nin,fileptr);
  fread (pars->fuu_2,sizeof(double),nin*nin,fileptr);
  fread (pars->fuv_0,sizeof(double),nin*nin,fileptr);
  fread (pars->fvv_0,sizeof(double),nin*nin,fileptr);
  fread (pars->fvv_1,sizeof(double),nin*nin,fileptr);
  fread (pars->fvv_2,sizeof(double),nin*nin,fileptr);
  fread (pars->ret1_0,sizeof(double),nin*nin,fileptr);
  fread (pars->ret1_1,sizeof(double),nin*nin,fileptr);
  fread (pars->ret2_0,sizeof(double),nin*nin,fileptr);
  fread (pars->ret2_1,sizeof(double),nin*nin,fileptr);
  fread (pars->ret3_0,sizeof(double),nin*nin,fileptr);
  fread (pars->ret3_1,sizeof(double),nin*nin,fileptr);
  fread (pars->phi1_0,sizeof(double),nin*nin,fileptr);
  fread (pars->phi2_0,sizeof(double),nin*nin,fileptr);
  fread (pars->phi3_0,sizeof(double),nin*nin,fileptr);
  fclose(fileptr);

  for (i=0;i<nin;i++) {
    d=(i-(nin/2-0.5))/(nin/2);
    for (j=0;j<nin;j++) {
      pars->xin[j*nin+i]=d; // Major cache thrashing!!!
      pars->yin[i*nin+j]=d;
    }
  }

  return 0;
}

int free_polcal(
  struct polcal_struct *pars
)
{
  MKL_free(pars->xin);
  MKL_free(pars->yin);
  MKL_free(pars->fqq_0);
  MKL_free(pars->fqq_1);
  MKL_free(pars->fqq_2);
  MKL_free(pars->fqu_0);
  MKL_free(pars->fqv_0);
  MKL_free(pars->fuu_0);
  MKL_free(pars->fuu_1);
  MKL_free(pars->fuu_2);
  MKL_free(pars->fuv_0);
  MKL_free(pars->fvv_0);
  MKL_free(pars->fvv_1);
  MKL_free(pars->fvv_2);
  MKL_free(pars->ret1_0);
  MKL_free(pars->ret1_1);
  MKL_free(pars->ret2_0);
  MKL_free(pars->ret2_1);
  MKL_free(pars->ret3_0);
  MKL_free(pars->ret3_1);
  MKL_free(pars->phi1_0);
  MKL_free(pars->phi2_0);
  MKL_free(pars->phi3_0);

  return 0;
}

int polcal(
  struct polcal_struct *pars,
  int nframe, // Number of input polarization states
  int mode, // 1 for IQUV, 2 for LCP+RCP from full, 3 for LCP+RCP from 2 pol
  float **input, // Pointers to input filtergrams
  float **output, // Pointers to output polarization images
  int *ps1, // PS1 positions
  int *ps2, // PS2 positions
  int *ps3, // PS3 positions
  float tsel, // Polarization selector temperature
  float tfront, // Front window temperatures
  int nx, // Number of colums
  int ny, // Number of rows
  int nlead // Declared number of rows. nlead>=nx
)
{
  int malign=32;
  int nin,polsol,polout;
  int i,j,k,ix,iy,iz,iframe,ipol;
  float *demod_in;
  double ts,tf;
  double dqq,duu,dvv,dqu,dqv,duv,dret1,dret2,dret3,dphi1,dphi2,dphi3;
  double *mtel,*mhelp,*m,*m1,*m2,*m3,*dmat;
  double *pl1angle,*pl2angle,*pl3angle;
  const double dir[]={-1,-1,1};
  const double pi=M_PI;
  const double pmax=0.938508*2048./1924.*pi/180;
  const double splitmax=0.01837;
  double x0,y0,dist0,phi0,phi,delta;
  double sum;
  char uplo[] = "U";
  int ifour=4;
  int info1,info2;
  float xscale,yscale,xzero,yzero;
  float xi,yi,*cx0,*cx1,cy0,cy1;
  int *ix0,iy0;
  float c00,c01,c10,c11,demod_out;
  int iz00,iz01,iz10,iz11;

  if ((mode<1) || (mode>3)) {
    printf("Bad output polarization in polcal: %d.\n",mode);
    return 1;
  }
  if (mode == 1) {
// Do full calibration and return IQUV
    polsol=4;
    polout=4;
  }
  if (mode == 2) {
// Do full calibration and return LCP+RCP
    polsol=4;
    polout=2;
  }
  if (mode == 3) {
// Solve for I and V and return LCP+RCP
    polsol=2;
    polout=2;
  }
  nin=pars->nin;
  ts=tsel-pars->tsela;
  tf=tfront-pars->tfronta;

  demod_in=(float *)MKL_malloc(nin*nin*nframe*polsol*sizeof(double),malign);
  mtel=(double *)MKL_malloc(16*sizeof(double),malign);
  mhelp=(double *)MKL_malloc(16*sizeof(double),malign);
  m=(double *)MKL_malloc(16*sizeof(double),malign);
  m1=(double *)MKL_malloc(16*sizeof(double),malign);
  m2=(double *)MKL_malloc(16*sizeof(double),malign);
  m3=(double *)MKL_malloc(16*sizeof(double),malign);
  dmat=(double *)MKL_malloc(4*nframe*sizeof(double),malign);
  pl1angle=(double *)MKL_malloc(nframe*sizeof(double),malign);
  pl2angle=(double *)MKL_malloc(nframe*sizeof(double),malign);
  pl3angle=(double *)MKL_malloc(nframe*sizeof(double),malign);

  for (i=0;i<nframe;i++) {
    pl1angle[i]=ps1[i]*1.5*pi/180*dir[0];
    pl2angle[i]=ps2[i]*1.5*pi/180*dir[1];
    pl3angle[i]=ps3[i]*1.5*pi/180*dir[2];
  }

// Initialize constant part of telescope matrix
  mtel[0]=1;mtel[1]=0;mtel[2]=0;mtel[3]=0;
  mtel[4]=0;mtel[8]=0;mtel[12]=0;
  for (iy=0;iy<nin;iy++) {
    for (ix=0;ix<nin;ix++) {
      iz=nin*iy+ix;
      x0=(ix-(nin/2-0.5))/(nin/2);
      y0=(iy-(nin/2-0.5))/(nin/2);
      dist0=sqrt(x0*x0+y0*y0)*pmax;
      phi0=atan2(y0,x0);

      dqq=pars->fqq_0[iz]+pars->fqq_1[iz]*tf+pars->fqq_2[iz]*tf*tf;
      duu=pars->fuu_0[iz]+pars->fuu_1[iz]*tf+pars->fuu_2[iz]*tf*tf;
      dvv=pars->fvv_0[iz]+pars->fvv_1[iz]*tf+pars->fvv_2[iz]*tf*tf;
      dqu=pars->fqu_0[iz];
      dqv=pars->fqv_0[iz];
      duv=pars->fuv_0[iz];
      
      dret1=pars->ret1_0[iz]+pars->ret1_1[iz]*ts;
      dret2=pars->ret2_0[iz]+pars->ret2_1[iz]*ts;
      dret3=pars->ret3_0[iz]+pars->ret3_1[iz]*ts;
      dphi1=pars->phi1_0[iz];
      dphi2=pars->phi2_0[iz];
      dphi3=pars->phi3_0[iz];

// Set variable part of telescope matriz
      mtel[5]=dqq;mtel[6]=dqu;mtel[7]=dqv;
      mtel[9]=dqu;mtel[10]=duu;mtel[11]=duv;
      mtel[13]=-dqv;mtel[14]=-duv;mtel[15]=dvv;

      for (iframe=0;iframe<nframe;iframe++) {
        phi=pl1angle[iframe]-dphi1;
        delta=dret1*(1+0.5*dist0*dist0*cos(2*(phi0-phi)));
        retarder(delta,phi,mhelp);
        mat4mm(mhelp,mtel,m1);
      
        phi=pl2angle[iframe]-dphi2;
        delta=dret2*(1+0.5*dist0*dist0*cos(2*(phi0-phi)));
        retarder(delta,phi,mhelp);
        mat4mm(mhelp,m1,m2);
      
        phi=pl3angle[iframe]-dphi3;
        delta=dret3*(1+0.5*dist0*dist0*cos(2*(phi0-phi)));
        retarder(delta,phi,mhelp);
        mat4mm(mhelp,m2,m3);

        phi=pi/2 + y0*splitmax;
        polarizer(phi,mhelp);
        mat4mm(mhelp,m3,m);

// Set dmat to matrix used for demodulation
        for (ipol=0;ipol<4;ipol++) dmat[4*iframe+ipol]=m[4*ipol];
      
      } // for iframe

      if (polsol==4) {
// Invert for all variables
        for (i=0;i<4;i++) {
          for (j=0;j<4;j++) {
            sum=0.0;
            for (iframe=0;iframe<nframe;iframe++) sum=sum+dmat[i+4*iframe]*dmat[j+4*iframe];
            mhelp[4*i+j]=sum;
          }
        }
      } // if (polsol==4)
     
// Only invert for I and V
      if (polsol==2) {
// Get rid of Q and U in dmat
        for (i=0;i<nframe;i++) dmat[4*i+1]=dmat[4*i+3];
        for (i=0;i<2;i++) {
          for (j=0;j<2;j++) {
            sum=0.0;
            for (iframe=0;iframe<nframe;iframe++) sum=sum+dmat[i+4*iframe]*dmat[j+4*iframe];
            mhelp[4*i+j]=sum;
          }
        }
      } // if (polsol==2)
     
      dpotrf(uplo,&polsol,mhelp,&ifour,&info1);
      dpotrs(uplo,&polsol,&nframe,mhelp,&ifour,dmat,&ifour,&info2);

/*
  if ((ix==15) && (iy==20)){
    printf("\n");
    for (i=0;i<polsol;i++) {
      for (j=0;j<polsol;j++) {
        printf("%f ",mhelp[4*j+i]);
      }
      printf("\n");
    }
    for (iframe=0;iframe<nframe;iframe++) {
      for (ipol=0;ipol<polsol;ipol++) printf("%f ",dmat[ipol+4*iframe]);
      printf("\n");
    }
  }
*/

      if (mode==1) {
        for (iframe=0;iframe<nframe;iframe++) {
          for (ipol=0;ipol<4;ipol++)
            demod_in[iframe+nframe*ipol+nframe*4*iz]=dmat[ipol+4*iframe];
        }
      } // if (mode==1)

      if (mode==2) {
        for (iframe=0;iframe<nframe;iframe++) {
          demod_in[iframe+nframe*0+nframe*2*iz]=dmat[4*iframe]+dmat[3+4*iframe];
          demod_in[iframe+nframe*1+nframe*2*iz]=dmat[4*iframe]-dmat[3+4*iframe];
        }
      } // if (mode==2)

      if (mode==3) {
        for (iframe=0;iframe<nframe;iframe++) {
          demod_in[iframe+nframe*0+nframe*2*iz]=dmat[4*iframe]+dmat[1+4*iframe];
          demod_in[iframe+nframe*1+nframe*2*iz]=dmat[4*iframe]-dmat[1+4*iframe];
        }
      } // if (mode==3)

/*
  if ((ix==15) && (iy==20)) {
    printf("\n");
    for (iframe=0;iframe<nframe;iframe++) {
      for (ipol=0;ipol<polout;ipol++) printf("%f ",demod_in[iframe+nframe*ipol+nframe*polout*iz]);
      printf("\n");
    }
  }
*/
      
    } // ix
  } // iy

  xscale=(nin+0.0)/nx;
  xzero=-xscale*(nx-1)/2+(nin-1.)/2;
  yscale=(nin+0.0)/ny;
  yzero=-yscale*(ny-1)/2+(nin-1.)/2;

  i=floor(xzero);
  j=floor(xscale);

  ix0=(int *)MKL_malloc(nx*sizeof(int),malign);
  cx0=(float *)MKL_malloc(nx*sizeof(float),malign);
  cx1=(float *)MKL_malloc(nx*sizeof(float),malign);
  for (ix=0;ix<nx;ix++) {
// Nor that numbers are extrapolated by up to 0.5 pixels on the edges
    xi=xzero+ix*xscale;
    ix0[ix]=minval(maxval(0,floor(xi)),(nin-2));
    cx1[ix]=xi-ix0[ix];
    cx0[ix]=1.0-cx1[ix];
  }

  for (iy=0;iy<ny;iy++) {
    yi=yzero+iy*yscale;
    iy0=minval(maxval(0,floor(yi)),(nin-2));
    cy1=yi-iy0;
    cy0=1.0-cy1;
    for (ix=0;ix<nx;ix++) {
      iz=nlead*iy+ix;
      iz00=nin*iy0+ix0[ix];
      iz01=nin*iy0+ix0[ix]+1;
      iz10=nin*(iy0+1)+ix0[ix];
      iz11=nin*(iy0+1)+ix0[ix]+1;
      c00=cy0*cx0[ix];
      c01=cy0*cx1[ix];
      c10=cy1*cx0[ix];
      c11=cy1*cx1[ix];
      for (ipol=0;ipol<polout;ipol++) {
        sum=0.0;
        for (iframe=0;iframe<nframe;iframe++) {
// First interpolate demodulation matrix
          demod_out=c00*demod_in[iframe+nframe*ipol+nframe*polout*iz00]+
                    c01*demod_in[iframe+nframe*ipol+nframe*polout*iz01]+
                    c10*demod_in[iframe+nframe*ipol+nframe*polout*iz10]+
                    c11*demod_in[iframe+nframe*ipol+nframe*polout*iz11];
// Then apply it
          sum=sum+demod_out*input[iframe][iz];
        }
        output[ipol][iz]=sum;
      }
    } // ix
  } // iy

  MKL_free(ix0);
  MKL_free(cx0);
  MKL_free(cx1);
  MKL_free(demod_in);
  MKL_free(mtel);
  MKL_free(mhelp);
  MKL_free(m);
  MKL_free(m1);
  MKL_free(m2);
  MKL_free(m3);
  MKL_free(dmat);
  MKL_free(pl1angle);
  MKL_free(pl2angle);
  MKL_free(pl3angle);

  return 0;
}

