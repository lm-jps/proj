
//all functions that are needed to calculate the rotational flatfield



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <jsoc_main.h>
#include <string.h>
#include <time.h>
#include </home/jsoc/include/fftw3.h>
#include <omp.h>
#include "module_flatfield.h"






long sign(double);
double mat_rh(long[], double[], double[], int);
void nine_diag(long[], double[], int, int, double[], long[], double);
void blockiter(double[], double[], double[], double*, double[], int, int, double);
void tridag(double veca[nx], double vecb[nx], double vecc[nx], double vecr[nx], double vecu[nx]);
void mat_mult(double[nx], double[nx], double[nx], double[nx], double[nx], int);
void printtime();
void highpass(int, int, double, double[nx][ny]);
void highpass_2d(int M, int N, double fwhm1, double fwhm2, double phi, double a[nx][ny]);
void derotation(double, double, double, double, double, double, double *, int, double *);
void derotation_test(double time, double radius, double cent_x, double cent_y, double p0, double b0, double *rotf);




int flatfield(double *rhsp, double *rhsm, short *badpix, int pairs, double *flati, double *param, struct code_param cpa,
               double deltat)
{

  
 
  double convergence=cpa.convergence;
  int maxiter=cpa.maxiter;
  double omega=cpa.omega;
  double croprad=cpa.croprad;
  double norm=cpa.norm;

  long count;
  int status;
  long i,j,k, l;
  int datum;
  double rdatum;
  double ddatum;
  double gout[nx*ny];
  

  for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) gout[j*nx+i]=flati[j*nx+i];
 
  int crop[ny][2];
  int xmin, xmax, ymin, ymax;
  long dirx, diry;

  double rhs_a[nx*ny+1];

  double aa_val[4*(nx-1)*ny+1];
  long aa_rowa[4*(nx-1)*ny+1], aa_rowb[4*(nx-1)*ny+1]; 
  long colarra[4*nx*ny], colarrb[4*nx*ny];

  double ata_vala[5*nx*ny], ata_valb[5*nx*ny], temp_ata[5*nx*ny];

  long aaro[4];
  double aava[4];

  double rh_a[nx*ny], rh_b[nx*ny];
 
  double gam0, gam1, gam2, gam3;

 
  double rsq;
       double gout_t[nx*ny];

       int iter;

 
	double res;
	double resdu[nx*ny];



	double R_SUN=*(param+0); //R_SUN in pixel
	double XX0=*(param+1);   //center x in pixel
	double YY0=*(param+2);   //cener y in pixel
	double P_ANG=*(param+3); //P-angle in rad
	double B_ANG=*(param+4); //B-angle in rad
	double dist=*(param+5);  //distance in AU


 // calculate rotation vector
	double *rot_coef, *rotf;
	rot_coef=(double *)(malloc(3*sizeof(double)));

	rot_coef[0]=cpa.rotcoef0;
	rot_coef[1]=cpa.rotcoef1;
	rot_coef[2]=cpa.rotcoef2;

	rotf=(double *)(malloc(nx*ny*2*sizeof(double)));

	//derotation(R_SUN, XX0, YY0, dist,  P_ANG, B_ANG, rot_coef, 2, rotf); 
	printf("rotpat param %lf %lf %lf %lf %lf %lf\n", deltat, R_SUN, XX0, YY0, P_ANG, B_ANG);
	derotation_test(deltat, R_SUN, XX0, YY0, P_ANG, B_ANG, rotf);
 
	//void derotation(double radius, double cent_x, double cent_y, double dist, double p0, double b0, double *rot_coef, int order2_rot_coef, double *shift, int nx, int ny){
    	double rota[nx][ny][2];
	for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) rota[i][j][0]=rotf[j*nx+i];
	for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) rota[i][j][1]=rotf[nx*ny+j*nx+i]; //rotf contains deltat




	


      //**************************calculate circular crop table (R_SUN*croprad)

	double ax1, ax2, ay1, ay2;


      ay1=YY0-(int)R_SUN*croprad;
      ay2=YY0+(int)R_SUN*croprad;

      ymin=(ay1 > 0 ? ay1:0);
      ymax=(ay2 < ny ? ay2:ny-1);
    
     
      for (j=0; j<ny; ++j){
	int xq=pow(R_SUN*croprad,2)-pow(j-YY0,2);
	if (xq > 0){ax1=XX0-sqrt(xq); ax2=XX0+sqrt(xq);  crop[j][0]=(ax1 > 0 ? ax1:0); crop[j][1]=(ax2 < nx ? ax2:nx-1);} else {crop[j][0]=XX0; crop[j][1]=XX0-1;}
      }

    
 
        const long dcount=4*(nx-1)*ny;
        const long drow=nx*ny;

	for (i=0; i<5*nx*ny; ++i){ata_vala[i]=0.0; ata_valb[i]=0.0;}

	// l-loop: forward (1) and backward (0) interpolation

	printtime();
    
	for (l=1; l>=0; --l){
	  count=0;

	for (i=0; i<(4*nx*ny); ++i){colarra[i]=dcount; colarrb[i]=dcount;}  //initialize colarr with dcount (=missing value)
        for (i=0; i<(4*(nx-1)*ny); ++i){aa_val[i]=0.0; aa_rowa[i]=0; aa_rowb[i]=0;}
        for (i=0; i<=nx*ny; ++i) rhs_a[i]=0.0;


	aa_val[dcount]=0.0;
	aa_rowa[dcount]=drow;
	aa_rowb[dcount]=drow;

	for (j=ymin; j<=ymax; ++j){
	  xmin=crop[j][0];
	  xmax=crop[j][1];
	 
     	  for (i=xmin; i<=xmax; ++i){
	  
            diry=(l*2-1)*sign(rota[i][j][1]);
	    dirx=(l*2-1)*sign(rota[i][j][0]);

	    if (!((i == xmax && dirx == 1) || (i == xmin && dirx == -1) || (j == ymax && diry == 1) || (j == ymin && diry == -1))){
	      //	      if (badpix[j*nx+i] && badpix[j*nx+i+dirx] && badpix[(j+diry)*nx+i] && badpix[(j+diry)*nx+i+dirx]){
	       
	
		aa_val[count]=-1.0+(1.0-fabs(rota[i][j][0]))*(1.0-fabs(rota[i][j][1]));
	   
	 
		aa_rowa[count]=j*nx+i;

	 
		aa_rowb[count]=i*ny+j;
             

		colarra[4*(j*nx+i)]=count;
		colarrb[4*(i*ny+j)]=count;
		count=count+1;
	      //
	    
		aa_val[count]=fabs(rota[i][j][0])*(1.0-fabs(rota[i][j][1]));
      
	 
		aa_rowa[count]=j*nx+i;

	 
		aa_rowb[count]=i*ny+j;

		colarra[4*(j*nx+i+dirx)+1]=count;
		colarrb[4*((i+dirx)*ny+j)+1]=count;
		count=count+1;
	      //
	      
		aa_val[count]=(1.0-fabs(rota[i][j][0]))*fabs(rota[i][j][1]);
             
	 
		aa_rowa[count]=j*nx+i;

	 
		aa_rowb[count]=i*ny+j;

		colarra[4*((j+diry)*nx+i)+2]=count;
		colarrb[4*(i*ny+j+diry)+2]=count;
		count=count+1;
	      //
	      
		aa_val[count]=fabs(rota[i][j][0])*fabs(rota[i][j][1]);
	     
	 
		aa_rowa[count]=j*nx+i;

	 
		aa_rowb[count]=i*ny+j;


		colarra[4*((j+diry)*nx+i+dirx)+3]=count;
		colarrb[4*((i+dirx)*ny+j+diry)+3]=count;
		count=count+1;
		// }
	     }
	  }
	}

 

	
	nine_diag(aa_rowa, aa_val, nx, ny, temp_ata, colarra, norm);
	for (i=0; i<5*nx*ny; ++i) ata_vala[i]=ata_vala[i]+temp_ata[i];   



        nine_diag(aa_rowb, aa_val, ny, nx, temp_ata, colarrb, norm);
	for (i=0; i<5*nx*ny; ++i) ata_valb[i]=ata_valb[i]+temp_ata[i];

	



	//right hand side


	  for (j=ymin; j<=ymax; ++j){
	     xmin=crop[j][0];
	     xmax=crop[j][1];
	      for (i=xmin; i<=xmax; ++i){
	        diry=(2*l-1)*sign(rota[i][j][1]);
	        dirx=(2*l-1)*sign(rota[i][j][0]);
		if (!((i == xmax && dirx == 1) || (i == xmin && dirx == -1) || (j == ymax && diry == 1) || (j == ymin && diry == -1))){
		  //	  if (badpix[j*nx+i] && badpix[j*nx+i+dirx] && badpix[(j+diry)*nx+i] && badpix[(j+diry)*nx+i+dirx]){

		    gam0=(1.0-fabs(rota[i][j][0]))*(1.0-fabs(rota[i][j][1]));
		    gam1=fabs(rota[i][j][0])*(1.0-fabs(rota[i][j][1]));
		    gam2=(1.0-fabs(rota[i][j][0]))*fabs(rota[i][j][1]);
		    gam3=fabs(rota[i][j][0])*fabs(rota[i][j][1]);

		    switch(l){

		    case 1: rhs_a[j*nx+i]=gam0*rhsp[j*nx+i]
		        	+gam1*rhsp[j*nx+(i+dirx)]
			        +gam2*rhsp[(j+diry)*nx+i]
			        +gam3*rhsp[(j+diry)*nx+(i+dirx)]
		                -rhsm[j*nx+i];
		      break;

                    case 0: rhs_a[j*nx+i]=gam0*rhsm[j*nx+i]
			                 +gam1*rhsm[j*nx+i+dirx]
			                 +gam2*rhsm[(j+diry)*nx+i]
			                 +gam3*rhsm[(j+diry)*nx+i+dirx]
			                 -rhsp[j*nx+i];
		      break;
		    }
		    // }
		}
	      }
	  }
	


	 
	      for (i=0; i<nx*ny; ++i){
		for (j=0; j<4; ++j){
	          aaro[j]=aa_rowa[colarra[4*i+j]];
	          aava[j]=aa_val[colarra[4*i+j]];
	        }
		switch(l){

		case 1: rh_a[i]=mat_rh(aaro, aava, rhs_a, 4); break;
		case 0: rh_a[i]=rh_a[i]+mat_rh(aaro, aava, rhs_a, 4); break;
	      
		}
	      }

	     
	      
	} /// end l-loop


 


	for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) rh_b[i*ny+j]=rh_a[j*nx+i];



   



       for (j=0; j<ny; ++j) for (i=0; i<nx; ++i){gout_t[i*ny+j]=gout[j*nx+i];}
 
    
	int icount=0; // reset counter

	if (fabs(rota[nx/2][ny/2][0]) >= fabs(rota[nx/2][ny/2][1])){ 
	  do {
	
	    blockiter(ata_vala, rh_a, gout, &res, resdu, nx, ny, omega);
	   	    
	      printf("%d\t%g a\n", icount, res);
	      ++icount;
	  }
	  while (res > convergence && icount <maxiter);
	}
	else
	  {
	    do{
	      blockiter(ata_valb, rh_b, gout_t, &res, resdu, ny, nx, omega);
	 
	    
	    printf("%d\t%g b\n", icount, res);
	    ++icount;
	    }
	    while (res > convergence && icount <maxiter);
	    
	    for (i=0; i<nx; ++i) for (j=0; j<ny; ++j) gout[j*nx+i]=gout_t[i*ny+j];
	  
	  }

	if (icount < maxiter) status=0; else status=1; 

	for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) flati[j*nx+i]=gout[j*nx+i];


	free(rot_coef);
	free(rotf);
	return status;
  
}





double mat_rh(long row_1[], double vec1[], double rhs[], int n)
//vector multiplication
{
  double res=0.0;
  int i;
  for (i=0; i<n; ++i) res=res+vec1[i]*rhs[row_1[i]];

  return res;
}


double mat_ata(long row_1[], double vec1[], long row_2[], double vec2[], int n)
{
  //vector dot product for sparse vectors
  double res=0.0;
  int i, c, kk;

  for (i=0; i<n; ++i){
    kk=-1;
    for (c=0; c<n; ++c) if (row_1[i] == row_2[c]) kk=c;
    if (kk != -1) res+=vec1[i]*vec2[kk];
  }
  return res;
}


void nine_diag(long aa_row[], double aa[], int nnx, int nny, double ata[], long colarr[], double norm)
{
  //performs A^T A matrix multiplication for for block-tridiag matrix
  long i, j, qi, qj, loc;
  long aaro1[4], aaro2[4];
  double aa1[4], aa2[4];
  long nn=nnx*nny;


  for (i=0; i<5*nx*ny; ++i) ata[i]=0.0;
  
  for (i=0; i<(nn-1); ++i){
    
    qi=i;
    qj=i;
   
    loc=i;
    for (j=0; j<4; ++j){aaro1[j]=aa_row[colarr[4*qi+j]]; aaro2[j]=aa_row[colarr[4*qj+j]]; aa1[j]=aa[colarr[4*qi+j]];  aa2[j]=aa[colarr[4*qj+j]];}
    ata[loc]=mat_ata(aaro1, aa1, aaro2, aa2, 4)+norm;

    //
    qi=i+1;
    for (j=0; j<4; ++j){aaro1[j]=aa_row[colarr[4*qi+j]]; aa1[j]=aa[colarr[4*qi+j]];}
 
    loc=1*nn+i;
    ata[loc]=mat_ata(aaro1, aa1, aaro2, aa2, 4);
    
    //
    if (i <= nnx*(nny-1)){
      qi=i+nnx-1;

      for (j=0; j<4; ++j){aaro1[j]=aa_row[colarr[4*qi+j]]; aa1[j]=aa[colarr[4*qi+j]];}

      loc=2*nn+i;
      ata[loc]=mat_ata(aaro1, aa1, aaro2, aa2, 4);
      
    }
    if (i <= nnx*(nny-1)-1){
      qi=i+nnx;
      for (j=0; j<4; ++j){aaro1[j]=aa_row[colarr[4*qi+j]]; aa1[j]=aa[colarr[4*qi+j]];}

      loc=3*nn+i;
      ata[loc]=mat_ata(aaro1, aa1, aaro2, aa2, 4);
     
    }
      //
      if (i <= nnx*(nny-1)-2){
	qi=i+nnx+1;
	for (j=0; j<4; ++j){aaro1[j]=aa_row[colarr[4*qi+j]]; aa1[j]=aa[colarr[4*qi+j]];}

	loc=4*nn+i;
	ata[loc]=mat_ata(aaro1, aa1, aaro2, aa2, 4);
	
      }
  }

  i=nnx*nny-1;
  qi=i;
  qj=i;

  for (j=0; j<4; ++j){aaro1[j]=aa_row[colarr[4*qi+j]]; aaro2[j]=aa_row[colarr[4*qj+j]]; aa1[j]=aa[colarr[4*qi+j]];  aa2[j]=aa[colarr[4*qj+j]];}

  loc=0*nn+i;
  ata[loc]=mat_ata(aaro1, aa1, aaro2, aa2, 4)+norm;

 
 	
}


long sign(double num)
//sign of num, zero is considered positive
{
  int rnum=floor((num*1000.0)+0.5);
  long out=1;
  if (rnum < 0) out=-1;
  if (rnum > 0) out=1;
  return out;

}



//long sign(double num)
//{
//  long out=1;
//  if (num > 0.0){out=1;}
//  if (num < 0.0){out=-1;}
//
//  return out;

//}


void blockiter(double ata[], double rh[], double x[], double  *res, double resdu[], int nnx, int nny, double omega)

//Gauss Seidel method for Block-Tridiagonal matrices
{
  long nn=nnx*nny;
  long j,i;
  double fa[nnx], fb[nnx], fc[nnx], fx[nnx], ex[nnx], ax[nnx];
  double aa[nnx], ab[nnx], ac[nnx];
  double xa[nnx], rhs[nnx];
  double xi[nnx];
  double xh[nn];

  double rhx[nn], rdatum;
 
  




  for (j=0; j<nny; ++j){

  for (i=0; i<nnx; ++i){fx[i]=0.0; ex[i]=0.0;}   

    if (j < (nny-1)){

      for (i=0; i<nnx; ++i){
	fa[i]=-ata[2*nn+j*nnx+i];
	fb[i]=-ata[3*nn+j*nnx+i];
	fc[i]=-ata[4*nn+j*nnx+i];
	xa[i]=x[(j+1)*nnx+i];
      }
      mat_mult(fa, fb, fc, xa, fx, nnx);
    }



    if (j > 0){
      fa[0]=0.0;
      fb[0]=-ata[3*nn+(j-1)*nnx];
      fc[0]=-ata[2*nn+(j-1)*nnx+1];
      xa[0]=x[(j-1)*nnx];

     
      for (i=1; i<(nnx-1); ++i){
	fa[i]=-ata[4*nn+(j-1)*nnx+i-1];
	fb[i]=-ata[3*nn+(j-1)*nnx+i];
	fc[i]=-ata[2*nn+(j-1)*nnx+i+1];
	
        xa[i]=x[(j-1)*nnx+i];
      }
	fa[nnx-1]=-ata[4*nn+j*nnx-2];
	fb[nnx-1]=-ata[3*nn+j*nnx-1];
	fc[nnx-1]=0.0;
	xa[nnx-1]=x[j*nnx-1];

      
    
     
      mat_mult(fa, fb, fc, xa, ex, nnx);
    }

 
   
 
  
   
    aa[0]=0.0;
    ab[0]=ata[j*nnx];
    ac[0]=ata[nn+j*nnx];
    xa[0]=x[j*nnx];

    for (i=1; i<nnx; ++i){
      aa[i]=ata[nn+j*nnx-1+i];
      ab[i]=ata[j*nnx+i];
      ac[i]=ata[nn+j*nnx+i];
      xa[i]=x[j*nnx+i];
    }

    mat_mult(aa, ab, ac, xa, ax, nnx);

    for (i=0; i<nnx; ++i) rhs[i]=omega*(fx[i]+ex[i]+rh[j*nnx+i])+(1-omega)*ax[i];

    tridag(aa, ab, ac, rhs, xi);
  
    

    
    for (i=0; i<nnx; ++i) x[j*nnx+i]=xi[i];
  }





    // calculate residuum


  for (j=0; j<nny; ++j){
  
    for (i=0; i<nnx; ++i){fx[i]=0.0; ex[i]=0.0; ax[i]=0.0;}  
 
    if (j < (nny-1)){
      for (i=0; i<nnx; ++i){
	fa[i]=ata[2*nn+j*nnx+i];
	fb[i]=ata[3*nn+j*nnx+i];
	fc[i]=ata[4*nn+j*nnx+i];
	xa[i]=x[(j+1)*nnx+i];
      }
      mat_mult(fa, fb, fc, xa, fx, nnx);
    }

    if (j > 0){
      fa[0]=0.0;
      fb[0]=ata[3*nn+(j-1)*nnx];
      fc[0]=ata[2*nn+(j-1)*nnx+1];
      xa[0]=x[(j-1)*nnx];

     
      for (i=1; i<(nnx-1); ++i){
	fa[i]=ata[4*nn+(j-1)*nnx+i-1];
	fb[i]=ata[3*nn+(j-1)*nnx+i];
	fc[i]=ata[2*nn+(j-1)*nnx+i+1];
	
        xa[i]=x[(j-1)*nnx+i];

      }

	fa[nnx-1]=ata[4*nn+j*nnx-2];
	fb[nnx-1]=ata[3*nn+j*nnx-1];
	fc[nnx-1]=0.0;
	xa[nnx-1]=x[j*nnx-1];
     
      mat_mult(fa, fb, fc, xa, ex, nnx);
    }

  
    aa[0]=0.0;
    ab[0]=ata[j*nnx];
    ac[0]=ata[nn+j*nnx];
    xa[0]=x[j*nnx];

    for (i=1; i<nnx; ++i){
      aa[i]=ata[nn+j*nnx-1+i];
      ab[i]=ata[j*nnx+i];
      ac[i]=ata[nn+j*nnx+i];

      xa[i]=x[j*nnx+i];

    }
    mat_mult(aa, ab, ac, xa, ax, nnx);

    for (i=0; i<nnx; ++i) xh[j*nnx+i]=fx[i]+ex[i]+ax[i]; 
  }

  *res=0.0;
  for (j=0; j<nny; ++j){
    for (i=0; i<nnx; ++i){
      // *res=*res+(double(xh[j*nnx+i])-double(rh[j*nnx+i]))*(double(xh[j*nnx+i])-double(rh[j*nnx+i]));
      *res=*res+fabs((double)xh[j*nnx+i]-(double)rh[j*nnx+i]);
         resdu[j*nnx+i]=xh[j*nnx+i]-rh[j*nnx+i];
    }
  }

}




void tridag(double veca[nx], double vecb[nx], double vecc[nx], double vecr[nx], double vecu[nx])
//Thomas algorithm to solve tridiagonal linear system
  {
    double bet, gam[nx];
    int j;
    
    if (vecb[0] == 0.0){printf("Error 1 in tridag\n"); return;};
  
    bet=vecb[0];
    vecu[0]=vecr[0]/bet;
    for (j=1;j<nx;++j){
      gam[j]=vecc[j-1]/bet;
      bet=vecb[j]-veca[j]*gam[j];
      if (bet == 0.0){printf("Error 2 in tridag\n");  return;};
      vecu[j]=(vecr[j]-veca[j]*vecu[j-1])/bet;
    }

    for (j=nx-2; j>=0; --j){
      vecu[j]=vecu[j]-gam[j+1]*vecu[j+1];
    }
  }


void mat_mult(double a[nx], double b[nx], double c[nx], double x[nx], double fx[nx], int n)
//multiplies a vector with a tridiagonal matrix
{
  //fx=Ax (A tridiagonal with vectors a, b, c)
  int i;
  double x1, x2, x3;
  for (i=0; i<n; ++i){
  
    x2=0.0;
    x3=0.0;

    x1=b[i]*x[i];
    if (i > 0) x2=a[i]*x[i-1];
    if (i < (n-1)) x3=c[i]*x[i+1];

    fx[i]=x1+x2+x3;
  }
}






void printtime()	
// print time
{
  time_t timer, timerd;
  char *timestring;
  int i;

	timerd=time(&timer);
	timestring=ctime(&timer);
        for (i=0; i<24; ++i) printf("%c", *(timestring+i));
	printf("\n");
}


 void derotation_test(double time, double radius, double cent_x, double cent_y, double p0, double b0, double *rotf)
{
    

  int i,j; //loop variables


   
    double svecp[2];
    int datum_int;

    double xy[2], inr;
    double xyp[3], xym[3], xys[2], norm;

    double slat, slat2, omeg;

    for (i=0; i<nx; ++i){
       for (j=0; j<ny; ++j){

	
	 xy[0]=((double)i+0.5-cent_x)/radius;
	 xy[1]=((double)j+0.5-cent_y)/radius;

         inr=sqrt(xy[0]*xy[0]+xy[1]*xy[1]);

	  if (inr <= 1.0)
	    {
	      // p-angle rotation
	      xyp[0]=cos(p0)*xy[0]+sin(p0)*xy[1];
	      xyp[1]=-sin(p0)*xy[0]+cos(p0)*xy[1];
	      // project onto sphere 
	      xym[0]=sqrt(1.0-(xyp[0]*xyp[0]+xyp[1]*xyp[1]));
              xym[1]=xyp[0];
	      xym[2]=xyp[1];
	      
              slat=sin(b0)*xym[0]+cos(b0)*xym[2];
              slat2=slat*slat;
	      omeg=2*M_PI*1e-9*(452. - 49.0*slat2 -84.0*slat2*slat2 - 31.7);

              svecp[0]=omeg*(xym[0]*cos(b0)-xym[2]*sin(b0))*time*radius;
              svecp[1]=omeg*xym[1]*sin(b0)*time*radius;

	    }

  
													     
	   else
	    {
              xyp[0]=cos(p0)*xy[0]+sin(p0)*xy[1];
	      xyp[1]=-sin(p0)*xy[0]+cos(p0)*xy[1];

	      norm=sqrt(xyp[0]*xyp[0]+xyp[1]*xyp[1]);
	      xys[0]=xys[0]/norm;
              xys[1]= xys[1]/norm;
	      
              slat=cos(b0)*xys[1];
              slat2=slat*slat;
              omeg=2*M_PI*1e-9*(452. - 49.0*slat2 -84.0*slat2*slat2 - 31.7);
	      
              svecp[0]=-omeg*xys[1]*time*radius*sin(b0);
              svecp[1]=omeg*xys[0]*time*radius*sin(b0);

	    }


	  rotf[j*nx+i]=cos(p0)*svecp[0]-sin(p0)*svecp[1];
	  rotf[nx*ny+j*nx+i]=sin(p0)*svecp[1]+cos(p0)*svecp[1];
      }
    }


  
}




void derotation(double radius, double cent_x, double cent_y, double dist, double p0, double b0, double *rot_coef, int order2_rot_coef, double *shift){
  //vector field of solar rotation in the image plain (linearized)

  int i, j, k;
  double xy[2], xyp[2];
  double xyz[3], sxyz[3];
  double sxyp[2], sxy[2];
  double inr, ind, ikf, phie, ing, beta;

   
  double slat, omeg;
  
      
    
   const double au=215.020*dist; // (1 astonomical unit in solar radii) au=149597870.691 km, R_SUN=695740 km (Kuhn 2004)
   const double maxphi=asin(1.0/au);
 
   
   double sinp0=sin(p0);
   double cosp0=cos(p0);
   double sinb0=sin(b0);
   double cosb0=cos(b0);

   int nlead=nx;


#pragma omp parallel for private(i,j,xyp,inr,phie,ind,ikf,beta,xyz,slat,omeg,sxyz,sxyp,sxy)
   for (j=0; j<ny; ++j){
      for (i=0; i<nx; ++i){
   
	  xyp[0]=((double)i-cent_x)/radius;
	  xyp[1]=((double)j-cent_y)/radius;


	  inr=sqrt(xyp[0]*xyp[0]+xyp[1]*xyp[1]); // radial coordinate in the sky
         
 
             
	  if (inr >= 1.0)
	    {
	      xyp[0]=xyp[0]/inr;
	      xyp[1]=xyp[1]/inr;  // force off-limb point to radial distance 1.0
	      inr=1.0;

	      phie=maxphi;
	      ind=sqrt(1.0-1.0/au/au);
            }
	  else
            {
	      phie=inr*maxphi;
	      beta=tan(phie);
              ind=(beta*au-beta*sqrt(1.0+beta*beta-au*au*beta*beta))/(1+beta*beta);
            }

	  if (inr != 0.0) ikf=ind/inr; else ikf=(au-1.0)/au; //limit for exact center


       	  xyz[0]=sqrt(1.0-ind*ind);
	  xyz[1]=xyp[0]*ikf;
	  xyz[2]=xyp[1]*ikf;

	  	      
	  slat=sinb0*xyz[0]+cosb0*xyz[2];  //sine of heliographic latitude (finite distance approximation)

		
	  omeg=1e-6*(rot_coef[1]*pow(slat,2)+rot_coef[2]*pow(slat,4)+rot_coef[0]);

	  sxyz[0]=-omeg*xyz[1]*cosb0;
	  sxyz[1]=omeg*(xyz[0]*cosb0-xyz[2]*sinb0); // linearized rotation
	  sxyz[2]=omeg*xyz[1]*sinb0; //
	 
	 
          sxyp[0]=sxyz[1]/ikf;
	  sxyp[1]=sxyz[2]/ikf;

	  sxy[0]=cosp0*sxyp[0]-sinp0*sxyp[1];
	  sxy[1]=sinp0*sxyp[0]+cosp0*sxyp[1];  // inverse p-angle rotation 

	  


	  shift[j*nlead+i]=sxy[0]*radius;
	  shift[ny*nlead+j*nlead+i]=sxy[1]*radius;
	   
     }
   }

}






void highpass(int M, int N, double fwhm, double a[nx][ny])
//highpass filter for 2d-array a
  {
 
    double b[M][N];
    fftw_complex ac[M][N/2+1], bc[M][N/2+1];
    fftw_plan p;
    long i,j;

    double sigma=fwhm/2.0/sqrt(2.0*log(2.0));
    double scale=1.0/((double)(M*N));
   

     //fft(a)
     p = fftw_plan_dft_r2c_2d(M, N, &a[0][0], &ac[0][0], FFTW_ESTIMATE);
     fftw_execute(p);
     fftw_destroy_plan(p);

     // Gaussian
     for (j = 0; j < N; ++j){
          for (i = 0; i < M; ++i) {
	    //b[i][j] = exp(-(((double)i-(double)M/2.0)*((double)i-(double)M/2.0)+((double)j-(double)N/2.0)*((double)j-(double)N/2.0))/2.0/sigma/sigma)/2.0/M_PI/sigma*sigma;
	    b[i][j] = exp(-(pow((double)i-(double)M/2.0,2)+pow((double)j-(double)N/2.0,2))/2.0/pow(sigma,2))/2.0/M_PI/pow(sigma,2);
            }
     }
     
     //fft(gaussian)
     p = fftw_plan_dft_r2c_2d(M, N, &b[0][0], &bc[0][0], FFTW_ESTIMATE);
     fftw_execute(p);
     fftw_destroy_plan(p);  


     //fft(a)*fft(b)
     for (i = 0; i < M; ++i){
          for (j = 0; j < N/2+1; ++j) {
	    ac[i][j]=ac[i][j]*bc[i][j]*scale;
          }
     }

     //  fft(fft(a)*fft(b),1)
     p = fftw_plan_dft_c2r_2d(M, N, &ac[0][0], &b[0][0], FFTW_ESTIMATE);
     fftw_execute(p);
     fftw_destroy_plan(p);

     //  highpass=data-lowpass
     for (i = 0; i < M; ++i){
          for (j = 0; j < N; ++j){
	    a[i][j]=a[i][j]-b[(i+M/2) % M][(j+N/2) % N];
	  }
     }


    

  }


void highpass_2d(int M, int N, double fwhm1, double fwhm2, double phi, double a[nx][ny])
//highpass filter for 2d-array a
  {
 
    double b[M][N];
    fftw_complex ac[M][N/2+1], bc[M][N/2+1];
    fftw_plan p;
    long i,j;
    double x,y;

    double sigma1=fwhm1/2.0/sqrt(2.0*log(2.0));
    double sigma2=fwhm2/2.0/sqrt(2.0*log(2.0));

    double scale=1.0/((double)(M*N));
    double sumgauss;

     //fft(a)
     p = fftw_plan_dft_r2c_2d(M, N, &a[0][0], &ac[0][0], FFTW_ESTIMATE);
     fftw_execute(p);
     fftw_destroy_plan(p);

     // Gaussian
     sumgauss=0.0;
     for (i = 0; i < N; ++i){
          for (j = 0; j < M; ++j) {
	    x=cos(phi/180.0*M_PI)*((double)i-(double)M/2.0)-sin(phi/180.0*M_PI)*((double)j-(double)N/2.0);
	    y=sin(phi/180.0*M_PI)*((double)i-(double)M/2.0)+cos(phi/180.0*M_PI)*((double)j-(double)N/2.0);

	    b[i][j] = exp(-(pow(x,2)/2.0/pow(sigma1,2)+pow(y,2)/2.0/pow(sigma2,2)));
	    sumgauss=sumgauss+b[i][j];
            }
     }

for (i = 0; i < N; ++i){
  for (j = 0; j < M; ++j){
    b[i][j]=b[i][j]/sumgauss;
  }
 }




     
     //fft(gaussian)
     p = fftw_plan_dft_r2c_2d(M, N, &b[0][0], &bc[0][0], FFTW_ESTIMATE);
     fftw_execute(p);
     fftw_destroy_plan(p);  


     //fft(a)*fft(b)
     for (i = 0; i < M; ++i){
          for (j = 0; j < N/2+1; ++j) {
	    ac[i][j]=ac[i][j]*bc[i][j]*scale;
          }
     }

     //  fft(fft(a)*fft(b),1)
     p = fftw_plan_dft_c2r_2d(M, N, &ac[0][0], &b[0][0], FFTW_ESTIMATE);
     fftw_execute(p);
     fftw_destroy_plan(p);

     //  highpass=data-lowpass
     for (i = 0; i < M; ++i){
          for (j = 0; j < N; ++j){
	    a[i][j]=a[i][j]-b[(i+M/2) % M][(j+N/2) % N];
	  }
     }


    

  }





void limb_darkening(double radius, double cent_x, double cent_y, double *b, int order, double *limb_dark)
{
//limb darkening function for solar intensity
  int i, j, k;
  double *mu;
  double rad;
  // const int order=5;

  mu=(double *)(malloc(nx*ny*sizeof(double)));

  for (i=0; i<nx; ++i){
    for (j=0; j<ny; ++j){
      rad=sqrt(pow((double)i-cent_x,2)+pow((double)j-cent_y,2))/radius;
      if (rad >=1.0)
	{
	  mu[j*nx+i]=1.0;
	}
      else
	{
	  if (rad <= 0.998) mu[j*nx+i]=sqrt(1.0-rad*rad);
	  if (rad > 0.998) mu[j*nx+i]=sqrt(0.002);
	}

      limb_dark[j*nx+i]=0.0;
    }
  }



  for (i=0; i<nx; ++i){
    for (j=0; j<ny; ++j){
      for (k=0; k<=order; ++k)  limb_dark[j*nx+i]=limb_dark[j*nx+i]+b[k]*pow(mu[j*nx+i],k);
    }
  }

  free(mu);
}





void apod_circ(float rad, float nb, float offx, float offy, float *vd)               
{
  float *rarr;
  rarr=(float *)(malloc(nx*ny*sizeof(float)));
  int i, j;

  for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) rarr[j*nx+i]=sqrt(((float)i-((float)nx/2+offx))*((float)i-((float)ny/2+offx))+((float)j-((float)nx/2+offy))*((float)j-((float)ny/2+offy)));
	 
 
  for (j=0; j<ny; ++j){
    for (i=0; i<nx; ++i){
      if (rarr[j*nx+i] < rad) 
	vd[j*nx+i]=1.0;

      if (rarr[j*nx+i] >= rad && rarr[j*nx+i] < (rad+nb))
	vd[j*nx+i]=0.5*cos(M_PI/nb*(rarr[j*nx+i]-rad))+0.5;

      if (rarr[j*nx+i] >= (rad+nb))
	vd[j*nx+i]=0.0;

      

    }
  }
	  
}


