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
#include "polcal.h"
#include "fresize.h"
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

int init_polcal(struct polcal_struct *pars, int method, const char *paramFile)
{
  int nin=32;
  int malign=32;
  int i,j;
  double d;
  FILE *fileptr = NULL;

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

  
  /* Remove hard-coding: fileptr = fopen ("/home/schou/hmi/anapol/pars_131222/fit.bin", "r"); */
  if (paramFile)
  {
    fileptr = fopen(paramFile, "r");
  }
  
  if (fileptr==NULL) {
    printf("File not found in init_polcal.\n");
    return 1;
  }

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
  double pscale;
  const float tscale=0.0f; // Scale factor to multiply on window temp diff.
  float fiq,fiu,fiv; // Leaks from I to Q, U and V
  int minps2,maxps2;
  int leakcor; // 0 to do nothing. 1 for constant correction. 2 for 4th order in r^2
  int psfcor; // 0 to do nothing. 1 do spatially constant correction.
  struct fresize_struct psf_q;
  struct fresize_struct psf_u;
  float *helpq,*helpu;
  float xx,yy,yy2,rr2x,*xx2s;
  float fiq0,fiq1,fiq2,fiq3,fiq4;
  float fiu0,fiu1,fiu2,fiu3,fiu4;
  float fiv0,fiv1,fiv2,fiv3,fiv4;
  float int0,int1,int2,int3,int4,int5;

// Arrays to get rid of I->Q,U leakage.
// New from averages of 0503, 0603, 0703, 0803, 0818, 8019 and 0903, 2010
// Values of cqa and cua psf_cor_100922.pro and psf_comb_100922.pro

/* Older version
  const float kerq_new[5][5] = 
{{0.0025447795,0.0014029023,-0.0068429443,-0.0036107973,0.0030788508},
{-0.0021226902,0.0024857402,-0.00038104460,0.0033163148,-0.00041170042},
{-0.0072690692,-0.0025919788,0.026803747,0.0042591058,-0.0095582193},
{-0.0018942680,0.0010363222,-0.0091301401,0.0011089003,-0.0036943740},
{0.0015510311,0.00070910644,-0.0033853462,0.0011713266,0.0013558073}};

  const float keru_new[5][5] =
{{-0.00050576922,0.0027637050,0.0021685222,-0.0048109307,-0.0017113943},
{-0.0037961382,-0.0036678872,-0.00043845913,-0.0025575728,0.0055228121},
{0.00068532645,-0.0033599591,0.018978875,-0.0026890615,0.0019440106},
{0.0031915392,-0.0020511113,-0.00015412126,-0.0033419297,-0.0048845723},
{-8.3558650e-05,-0.0034679275,0.0010754580,0.0022082227,-0.0010739051}};

*/
// Same data but now with version forcing total(kernel)=0
// Using psf_cor_100923 and psf_comb_100923. cqa and cua.

  const float kerq_new[5][5] =
{{0.0025540328,0.0014054127,-0.0068338746,-0.0036082337,0.0030879086},
{-0.0021203885,0.0024757898,-0.00038367994,0.0033062588,-0.00040921161},
{-0.0072601518,-0.0025927805,0.026828298,0.0042583267,-0.0095492689},
{-0.0018917886,0.0010262981,-0.0091327691,0.0010989052,-0.0036920641},
{0.0015600444,0.00071167337,-0.0033762701,0.0011738308,0.0013651198}};

  const float keru_new[5][5] =
{{-0.00049870912,0.0027657090,0.0021755355,-0.0048088858,-0.0017044651},
{-0.0037942851,-0.0036754557,-0.00044049794,-0.0025652278,0.0055248008},
{0.00069222959,-0.0033606194,0.018997891,-0.0026897083,0.0019509337},
{0.0031935224,-0.0020587455,-0.00015615637,-0.0033495235,-0.0048827137},
{-7.6658680e-05,-0.0034658817,0.0010824736,0.0022102219,-0.0010668119}};


// Old from averages of 100419 and 100416
  const float kerq_old[5][5] =
{{0.0037206091,-0.0019643126,-0.016193939,-0.0076218857,0.0043843131},
{0.00049905806,0.0040439727,-0.017494203,0.0024333910,0.0030044689},
{-0.0067546715,0.0080826366,0.062275348,0.010767046,-0.0089655220},
{-0.0024383606,0.00057499661,-0.023057758,0.00071254324,-0.0036152480},
{0.0026420443,-0.0017713332,-0.013929130,-0.0018103380,0.0024885890}};

  const float keru_old[5][5] =
{{0.0018666587,0.0075391244,0.00049939188,-0.0097649409,-0.0032506137},
{-0.0033866379,0.0057003782,-0.0046985594,-0.012012604,0.0052348375},
{-0.00060864512,-0.0028413328,0.023487552,-0.0011891184,0.0022289085},
{0.0029176405,-0.013011020,-0.0013825372,0.0083993347,-0.0044858224},
{0.00094548694,-0.0086396909,-0.00012516038,0.0078249956,-0.0012269566}};

// New averages not made.

// Set up instrumental polarization correction.
// Probably is actually something else.
// Only correct for mode==1 (IQUV).
  leakcor=0;
  psfcor=0;
  fiq=0.0f;
  fiu=0.0f;
  fiv=0.0f;
  minps2=300;
  maxps2=-10;
  for (i=0;i<nframe;i++) {
    minps2=minval(minps2,ps2[i]);
    maxps2=maxval(maxps2,ps2[i]);
  }
  if ((minps2 == 91) && (maxps2 == 91)) { // Old settings
// Stick to simple correction here.
    leakcor=1;
    fiq= 2.15e-4;
    fiu= 1.10e-4;
    fiv= 0.00e-4;
    psfcor=1;
    init_fresize_user(&psf_q,2,1,(float *)kerq_old); // Cast to get rid of idiotic compiler warning
    init_fresize_user(&psf_u,2,1,(float *)keru_old);
  }
  if ((minps2 == 53) && (maxps2 == 53)) { // New settings
/* Constant correction. No longer used.
    leakcor=1;
    fiq=-0.10e-4;
    fiu=-1.30e-4;
    fiv= 0.00e-4;
*/
    leakcor=2;
// Coefficients from fit_rdep_100922a
    fiq0=   2.0411683e-06;
    fiq1=   0.00012040193;
    fiq2=   0.00026537064;
    fiq3=  -0.00059574417;
    fiq4=  -0.0026085394;

    fiu0=  -0.00011803840;
    fiu1=   0.00011883301;
    fiu2=   0.00023335908;
    fiu3=  -0.0011413839;
    fiu4=  -0.0034897879;

// Note that V is currently not corrected!
    fiv0=   4.2255156e-06;
    fiv1=   1.7374540e-05;
    fiv2=   1.3426406e-05;
    fiv3=   0.00012586214;
    fiv4=   0.00045829690;
    fiv0=sqrt(-1.0); // Just in case someone tries doing V by accident

    psfcor=1;
    init_fresize_user(&psf_q,2,1,(float *)kerq_new);
    init_fresize_user(&psf_u,2,1,(float *)keru_new);
  }
  if (minps2 != maxps2) {
    printf("Variable PS2 positions. Not correcting for instrumental polarization.\n");
  }
  if ((minps2 != 53) && (minps2 !=91)) {
    printf("Non-standard PS2 positions. Not correcting for instrumental polarization.\n");
  }

  pscale=0.5; // Scaling of polarization. 0.5 to make I equal to the input intensity.
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
  tf=tscale*(tfront-pars->tfronta);

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

#pragma omp parallel for default(none) \
private(ix,iy,iz,yi,iy0,cy0,cy1,iz00,iz01,iz10,iz11,c00,c01,c10,c11,ipol,sum,iframe,demod_out) \
shared(nlead,output,nx,ny,yzero,yscale,nin,ix0,cx0,cx1,nframe,polout,demod_in,pscale,fiq,fiu,fiv,mode,input)
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
        output[ipol][iz]=pscale*sum;
      } // ipol
/* Moved to loop below
      if (mode == 1) {
        output[1][iz]=output[1][iz]-fiq*output[0][iz];
        output[2][iz]=output[2][iz]-fiu*output[0][iz];
        output[3][iz]=output[3][iz]-fiv*output[0][iz];
      }
*/
    } // ix
  } // iy

  if ((mode == 1) && (leakcor == 1)) { // Do simple correction for leak.
#pragma omp parallel for default(none) \
private(ix,iy,iz) \
shared(nlead,output,helpq,helpu,nx,ny,fiq,fiu,fiv)
    for (iy=0;iy<ny;iy++) {
      for (ix=0;ix<nx;ix++) {
        iz=nlead*iy+ix;
        output[1][iz]=output[1][iz]-fiq*output[0][iz];
        output[2][iz]=output[2][iz]-fiu*output[0][iz];
        output[3][iz]=output[3][iz]-fiv*output[0][iz];
      } // ix
    } // iy
  }

  if ((mode == 1) && (leakcor == 2)) { // Leak is given by 4th order in r^2.
    xx2s=(float *)MKL_malloc(nx*sizeof(float),malign); // Array to hold x values.
    for (ix=0;ix<nx;ix++) {
      xx=(ix-nx/2+0.5f)/(nx/2);
      xx2s[ix]=xx*xx;
    }
#pragma omp parallel for default(none) \
private(ix,iy,iz,yy,yy2,rr2x) \
shared(int0,int1,int2,int3,int4) \
shared(fiq0,fiq1,fiq2,fiq3,fiq4) \
shared(fiu0,fiu1,fiu2,fiu3,fiu4) \
shared(fiv0,fiv1,fiv2,fiv3,fiv4) \
shared(nlead,output,helpq,helpu,nx,ny,fiq,fiu,fiv,xx2s)
    for (iy=0;iy<ny;iy++) {
      yy=(iy-ny/2+0.5f)/(ny/2);
      yy2=yy*yy;
      for (ix=0;ix<nx;ix++) {
        iz=nlead*iy+ix;
        rr2x=xx2s[ix]+yy2-0.5f; // Subtract 0.5 to limit coefficients.
// Make intensities times power of radius.
        int0=output[0][iz];
        int1=rr2x*int0;
        int2=rr2x*int1;
        int3=rr2x*int2;
        int4=rr2x*int3;
        output[1][iz]=output[1][iz]-fiq0*int0-fiq1*int1-fiq2*int2-fiq3*int3-fiq4*int4;
        output[2][iz]=output[2][iz]-fiu0*int0-fiu1*int1-fiu2*int2-fiu3*int3-fiu4*int4;
// For now don't do V
//      output[3][iz]=output[3][iz]-fiv0*int0-fiv1*int1-fiv2*int2-fiv3*int3-fiv4*int4;
      } // ix
    } // iy
    MKL_free(xx2s);
  }

  if ((mode == 1) && (psfcor == 1)) {
    helpq=(float *)MKL_malloc(ny*nlead*sizeof(float),malign);
    helpu=(float *)MKL_malloc(ny*nlead*sizeof(float),malign);
    fresize(&psf_q,output[0],helpq,nx,ny,nlead,nx,ny,nlead,0,0,0.0f);
    fresize(&psf_u,output[0],helpu,nx,ny,nlead,nx,ny,nlead,0,0,0.0f);
#pragma omp parallel for default(none) \
private(ix,iy,iz) \
shared(nlead,output,helpq,helpu,nx,ny)
    for (iy=0;iy<ny;iy++) {
      for (ix=0;ix<nx;ix++) {
        iz=nlead*iy+ix;
        output[1][iz]=output[1][iz]-helpq[iz];
        output[2][iz]=output[2][iz]-helpu[iz];
      } // ix
    } // iy
    MKL_free(helpq);
    MKL_free(helpu);
  }

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

char *polcal_version() // Returns CVS version of polcal.c
{
  return strdup("$Id: polcal.c,v 1.6 2016/01/21 20:30:48 arta Exp $");
}

