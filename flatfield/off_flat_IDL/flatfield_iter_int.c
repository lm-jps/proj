

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>







const int niter=10;
const float omega=1.2;
const int residuum=0;
const int n_foc=1;



void apod_circ(int, float, float, int, float, float, float *, int , int);
void printtime();

int main(int argc, char *argv[]){

  int camera;
  int nx, ny;
  int nl;
  float fac,facr;

  if (argc != 6){fputs("Wrong number of arguments", stderr); exit(1);}

  nx=atoi(argv[1]);
  ny=atoi(argv[1]);
  nl=atoi(argv[2]);
  camera=atoi(argv[3]);
  fac=(float)atoi(argv[4])/100.0;
  facr=(float)atoi(argv[5])/100.0;

  printf("nx ny nl cam %d %d %d %d\n", nx, ny, nl, camera);

  float samp=4096.0/(float)nx;
 

  float **imr; // image variable
  float *imp;
  float *imm;

  float *imr_im, *imr_jm;
  
  //float g[nx][ny];        // flatfield iteration
  //float kap[nx][ny];       // initial value

  float *vdo, *vdl;  // 
  float *vdp; 

  vdo=(float *)(malloc(nx*ny*sizeof(float)));
  vdl=(float *)(malloc(nx*ny*sizeof(float)));
  vdp=(float *)(malloc(nx*ny*nl*sizeof(float)));
  
  float xx[nl], yy[nl];

  int nbad;
  //side
  // int nbad=32;
  //  int ix_bad[32];

  //  int ix_bad_front[32]={738657 ,    1100979 ,    1100980 ,    1420552 ,    1424648   ,  1428744   ,  1502539 ,    1502540  ,   1506635 ,    1506636  ,   5961670   ,  6051854, 6330150 ,   10171860,    10175955 ,   10175956};

  //  int ix_bad_side[32]={1800075     ,1800076     ,1804173     ,9958325     ,9958326     ,9958327     ,9962423     ,9962424     ,9962425     ,9966520    ,12176519    ,12176520 ,12176521    ,12180613    ,12180614    ,12180615    ,12180616    ,12180617    ,12184709   , 12184710    ,12184711    ,12184712    ,12184713    ,12188806 ,12188807   , 12188808   , 12188809    ,12426166    ,12426167    ,12430263   , 12430264    ,12434360};



 



  //float n[nx][ny];
  //float gn[nx][ny];
  //float nn[nx][ny];


  float *g, *kap, *n, *gn;
  g=(float *)(malloc(nx*ny*sizeof(float)));
  kap=(float *)(malloc(nx*ny*sizeof(float)));
 
  n=(float *)(malloc(nx*ny*sizeof(float)));
  gn=(float *)(malloc(nx*ny*sizeof(float)));
 

  float res[niter];
 
  // reading images
  float datum;
  int im_number;

  float delta;

  int i, j, k, im, jm, iter;  // loop variables

  int valfoc[7]={0,2,3,4,5,9,13}; // 0 replacing 1 



 //front
  //  if (camera == 2){nbad=16;
  //    for (k=0; k<nbad; ++k) ix_bad[k]=ix_bad_front[k];
  //  }

  //side camera
  //  if (camera == 1){nbad=32;
  //    for (k=0; k<nbad; ++k) ix_bad[k]=ix_bad_side[k];
  //  }


  printf("reading data\n");
  int fi;
  for (fi=0; fi < n_foc; ++fi){  //0 replacing 6
 
    //reading data

  
    ///get input name
    char *fname="imr_hmi0.bin";
 
     //////////////////////////
   // get output flatfield name
 


    
    char  *ffname="flatfield_out0.bin";

   
    

    printf("output name %s \n", ffname);
    ////////////////////////////// 

    imm=(float *)(malloc(nx*sizeof(float)));
   imr=(float **)(malloc(nl*sizeof(float*)));
    
    FILE *imfile;
    size_t bytes_read;
    imfile=fopen(fname, "rb"); // read binary data

    if (imfile==NULL){fputs("File error", stderr); exit(1);};
    {
      for (k=0; k<nl; ++k){
	printf("%u \n", k);
	imp=(float *)(malloc(nx*ny*sizeof(float)));
        *(imr+k)=imp;
	bytes_read=fread(imp,sizeof(float), nx*ny,imfile);
      }
    }
    fclose(imfile);

 
    //get log of input data


    // reading leg positions
    FILE *legpos;
    legpos=fopen("legpos2", "r");  

    if (legpos != NULL){
      for (i=0; i<nl; i++){
	fscanf(legpos, "%f", &datum);
	xx[i]=datum/samp;
	printf("%f\n", xx[i]);
      }
      printf("\n");
      for (i=0; i<nl; i++){
	fscanf(legpos, "%f", &datum);
	yy[i]=datum/samp;
	printf("%f\n", yy[i]);
      }
    }
    fclose(legpos);



    // constructing apodization

    int low9=0; // int(192.0/samp);
    float nr=nx/2*fac;
    float nb=nx/2*fac*0.0; // nb=0 !!

    float nrs=nx/2*facr*0.95;
    float nbs=nx/2*facr*0.05;

    apod_circ(nx, nr, nb, 0, 0.0, 0.0, vdl,nx,ny);   


    //   for (k=0; k<nbad; ++k) vdl[ix_bad[k]]=0.0;
      
    for (k=0; k<nl; ++k){
      printf("%d\n", k);
      apod_circ(nx, nrs, nbs, 0, xx[k], yy[k], vdo,nx,ny); 
          
      for (j=0; j<nx; ++j) for (i=0; i<ny; ++i) vdp[k*nx*ny+j*nx+i]=vdl[j*nx+i]*vdo[j*nx+i];
    }

       
  
    // zero iteration

    int xf, yf, xg, yg;
    float rx, ry;
    float xm, ym;
    float xn, yn;
  



  int nthreads;
  //nthreads=omp_get_num_procs();                                      //number of threads supported by the machine where the code is running
  nthreads=1;
  omp_set_num_threads(nthreads);                                     //set the number of threads to t
  printf("Number of threads run in parallel by the subroutine= %d \n",nthreads);


    printtime();

    printf("zero iteration\n");

    for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) kap[j*nx+i]=0.0;

#pragma omp parallel for private(i,j,im,jm,xm,xf,ym,yf,yn,yg,xn,xg)
 for (j=0; j<ny; ++j){
   for (i=0; i<nx; ++i){
     n[j*nx+i]=0.0;
     for (jm=0; jm<nl; ++jm){ //exchanged loops
      for (im=jm+1; im<nl; ++im){

	imr_im=*(imr+im);         
	imr_jm=*(imr+jm);
  
	    xm=i-xx[im]+xx[jm];
	    xf=(int)(xm+0.5);
	    ym=j-yy[im]+yy[jm];
	    yf=(int)(ym+0.5);

	    if (xf >= 0 && xf <= nx-2 && yf >= 0 && yf <= ny-2){
	      if (vdp[im*nx*ny+j*nx+i] != 0.0 && vdp[jm*nx*ny+yf*nx+xf] != 0.0 && imr_im[j*nx+i] > 0.0 && imr_jm[yf*nx+xf] > 0.0){
		  
		kap[j*nx+i]=kap[j*nx+i]+(vdp[im*nx*ny+j*nx+i])*(vdp[jm*nx*ny+yf*nx+xf])*(log(imr_im[j*nx+i])-log(imr_jm[yf*nx+xf]));
 		n[j*nx+i]=n[j*nx+i]+(vdp[im*nx*ny+j*nx+i])*(vdp[jm*nx*ny+yf*nx+xf]);
		  
	      }
	    }
 
	      yn=j+yy[im]-yy[jm];
	      yg=(int)(yn+0.5);
	      xn=i+xx[im]-xx[jm];
	      xg=(int)(xn+0.5);

	      if (xg >= 0 && xg <= nx-2 && yg >=0 && yg <= ny-2){
		
		if (vdp[jm*nx*ny+j*nx+i] != 0.0 && vdp[im*nx*ny+yg*nx+xg] !=0.0 && imr_jm[j*nx+i] > 0.0 && imr_im[yg*nx+xg] > 0.0){
                
		  		 
		  kap[j*nx+i]=kap[j*nx+i]+(vdp[jm*nx*ny+j*nx+i])*(vdp[im*nx*ny+yg*nx+xg])*(log(imr_jm[j*nx+i])-log(imr_im[yg*nx+xg]));

		  n[j*nx+i]=n[j*nx+i]+(vdp[jm*nx*ny+j*nx+i])*(vdp[im*nx*ny+yg*nx+xg]);
		}
	      }
      }
    } 							    
	
    if (n[j*nx+i] > 0.0) kap[j*nx+i]=kap[j*nx+i]/n[j*nx+i]; else kap[j*nx+i]=0.0;

   }
 }


    
 

    for (j=0; j<ny; ++j) for (i=0; i<nx; ++i) g[j*nx+i]=kap[j*nx+i];
    for (i=0; i<niter; ++i) res[i]=0.0;
	 

    printtime();

    // iterations
   
    for (iter=0; iter<niter; ++iter){
	    
      printf("iteration %u \n", iter+1);
     
#pragma omp parallel for private(i,j,im,jm,xm,xf,ym,yf,yn,yg,xn,xg)
	  for (j=0; j<ny; ++j){
	      for (i=0; i<nx; ++i){
		n[j*nx+i]=0.0;
		gn[j*nx+i]=0.0;
	

		for (jm=0; jm<nl; ++jm){
		  for (im=jm+1; im<nl; ++im){

		    xm=i-xx[im]+xx[jm];
		    xf=(int)(xm+0.5);
		    ym=j-yy[im]+yy[jm];
		    yf=(int)(ym+0.5);

		    if (yf >= 0 && yf <= ny-2 && xf >= 0 && xf <= nx-2){
		      if (vdp[im*nx*ny+j*nx+i] != 0.0 && vdp[jm*nx*ny+yf*nx+xf] !=0.0){

			n[j*nx+i]=n[j*nx+i]+vdp[im*nx*ny+j*nx+i]*vdp[jm*nx*ny+yf*nx+xf];
			gn[j*nx+i]=gn[j*nx+i]+vdp[im*nx*ny+j*nx+i]*vdp[jm*nx*ny+yf*nx+xf]*g[yf*nx+xf];
		  

		  }
		}
	      
	    
	  
		    xn=i+xx[im]-xx[jm];
		    xg=(int)(xn+0.5);
		    yn=j+yy[im]-yy[jm];
		    yg=(int)(yn+0.5);

		    if (yg >= 0 && yg <= ny-2 && xg >= 0 && xg <= nx-2){
		      if (vdp[jm*nx*ny+j*nx+i] !=0.0 && vdp[im*nx*ny+yg*nx+xg] !=0.0){
		
			  n[j*nx+i]=n[j*nx+i]+vdp[jm*nx*ny+j*nx+i]*vdp[im*nx*ny+yg*nx+xg];
			  gn[j*nx+i]=gn[j*nx+i]+vdp[jm*nx*ny+j*nx+i]*vdp[im*nx*ny+yg*nx+xg]*g[yg*nx+xg];
			 

		      }
		    }
	      
		  }
		}
		  
		if (n[j*nx+i] > 0.0) gn[j*nx+i]=gn[j*nx+i]/n[j*nx+i]; else gn[j*nx+i]=0.0;
		delta=kap[j*nx+i]-g[j*nx+i]+gn[j*nx+i];
		g[j*nx+i]=g[j*nx+i]+omega*delta;

		
	      }
	  }


	
      if (residuum){
	for (j=0; j<ny; ++j){
	  for (i=0; i<nx; ++i){

	    for (jm=0; jm<(nl-1); ++jm){
	      for (im=jm+1; im<nl; ++im){

		imr_im=*(imr+im);        
		imr_jm=*(imr+jm);

		xf=(int)(i+xx[im]);
		yf=(int)(j+yy[im]);
		xg=(int)(i+xx[jm]);
		yg=(int)(j+yy[jm]);

		if (xf >=0 && xf <nx && xg >=0 && xg <nx && yf>=0 && yf<ny && yg >=0 && yg < ny){
		  if (vdp[im*nx*ny+yf*nx+xf]*vdp[jm*nx*ny+yg*nx+xg] !=0.0 && imr_im[yf*nx+xf] > 0.0 && imr_jm[yg*nx+xg] > 0.0){
		    res[iter]=res[iter]+vdp[im*nx*ny+yf*nx+xf]*vdp[jm*nx*ny+yg*nx+xg]*pow(-g[yf*nx+xf]+g[yg*nx+xg]+log(imr_im[yf*nx+xf])-log(imr_jm[yg*nx+xg]),2);
		  }
		}

	      }
	    }

	  }
	}
	  
	printf("residuum: %f\n", res[iter]);
      }
    }

 
    printtime();


 FILE *outname;
     outname = fopen (ffname, "w");
     for (j=0;j<ny;j++){
       for (i=0; i<nx; ++i) imm[i]=exp(g[j*nx+i]);
       fwrite ((char*)(imm),sizeof(float),nx,outname);
     }
     fclose(outname);


	
    for (k=0; k<nl; ++k) free(*(imr+k));
    free(imr);       
  }
  
  free(imm);

  return 0;
}     


void apod_circ(int nn, float rad, float nb, int low9, float offx, float offy, float *vd, int nx, int ny)
{
  float *rarr;
  rarr=(float *)(malloc(nx*ny*sizeof(float)));
  int i, j;

  for (i=0; i<nn; ++i) for (j=0; j<nn; ++j) rarr[j*nx+i]=sqrt(((float)i-((float)nn/2+offx))*((float)i-((float)nn/2+offx))+((float)j-((float)nn/2+offy))*((float)j-((float)nn/2+offy)));
	 
  for (i=0; i<nn; ++i){
    for (j=0; j<nn; ++j){

      if (rarr[j*nx+i] < rad) 
	vd[j*nx+i]=1.0;

      if (rarr[j*nx+i] >= rad && rarr[j*nx+i] < (rad+nb))
	vd[j*nx+i]=0.5*cos(M_PI/nb*(rarr[j*nx+i]-rad))+0.5;

      if (rarr[j*nx+i] >= (rad+nb))
	vd[j*nx+i]=0.0;

      if (low9 != 0){
	     
	if (j < low9) vd[j*nx+i]=0.0;
	if (j >= low9 && j <= (low9+(int)nb)) vd[j*nx+i]=(-0.5*cos(M_PI/nb*(j-low9))+0.5)*vd[j*nx+i];
      }

 


    }
  }
	  
}


			      

void printtime()	// print time
{
  time_t timer, timerd;
  char *timestring;
  int i;

	timerd=time(&timer);
	timestring=ctime(&timer);
        for (i=0; i<24; ++i) printf("%c", *(timestring+i));
	printf("\n");
}

			      
