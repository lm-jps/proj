/* 
	I.Scholl 

		from Solar Astrometry Program (solarSDO.c)
				Marcelo Emilio  Jan 2010 


	#define CODE_NAME 		"limbfit"
	#define CODE_VERSION 	"V0.1r8" 
	#define CODE_DATE 		"Fri Apr  9 20:20:36 HST 2010" 
*/

#include "limbfit.h"

int limbfit(LIMBFIT_INPUT *info,LIMBFIT_OUTPUT *results,int skipc,int debug)
{
printf(" in limfit\n");

/************************************************************************
					INIT VARS coming from build_lev1.c
************************************************************************/	
#ifdef HMI
	float *data=info->data;
//	printf("hmi selected\n");
//	float MIN_VALUE=BAD_PIXEL_VALUE;
//	float MAX_VALUE=(float)fabs(BAD_PIXEL_VALUE);
#elif  MDI
	short *data=info->data;
//	printf("mdi selected\n");
//	short MIN_VALUE=BAD_PIXEL_VALUE;
//	short MAX_VALUE=(short)abs(BAD_PIXEL_VALUE);
#endif

int 	ret_code = 0;

int 	w		 = info->annw;
long 	S		 = info->mxsz_annv;
int 	nang	 = info->nb_ldf;
int 	nprf	 = info->nb_rdb;
int 	nreg	 = info->nb_abb;
float 	rsi		 = info->lol;
float 	rso		 = info->upl;
int 	r_size	 = info->nb_fit_pnts;
int		grange	 = info->guess_rng;
float 	dx		 = info->incx;
float 	dy		 = info->incy;
int		naxis_row= info->img_szx;
int		naxis_col= info->img_szy;
int		iter	 = info->iter;
//float 	flag	 = BAD_PIXEL_VALUE;

if (debug) 
{
	fprintf(info->fd,"in limbfit annw: %d \n",info->annw);
	fprintf(info->fd,"in limbfit mxsz_annv: %ld \n",info->mxsz_annv);
	fprintf(info->fd,"in limbfit nb_ldf: %d \n",info->nb_ldf);
	fprintf(info->fd,"in limbfit nb_rdb: %d \n",info->nb_rdb);
	fprintf(info->fd,"in limbfit nb_abb: %d \n",info->nb_abb);
	fprintf(info->fd,"in limbfit nb_fit_pnts: %d \n",info->nb_fit_pnts);
	fprintf(info->fd,"in limbfit guess range: %d \n",info->guess_rng);
	fprintf(info->fd,"in limbfit lol: %f\n",info->lol);
	fprintf(info->fd,"in limbfit upl: %f\n",info->upl);
	fprintf(info->fd,"in limbfit incx: %f\n",info->incx);
	fprintf(info->fd,"in limbfit incy: %f\n",info->incy);
	fprintf(info->fd,"in limbfit szx: %d\n",info->img_szx);
	fprintf(info->fd,"in limbfit szy: %d\n",info->img_szy);
	fprintf(info->fd,"in limbfit bad pix: %f\n",BAD_PIXEL_VALUE);
	fprintf(info->fd,"end init............\n");
} 
	
/************************************************************************/
/*                        set parameters                               */
/************************************************************************/

long 	ii, jj, jk, i,j; 
float cmx, cmy,r, guess_cx, guess_cy, guess_r;
int nitr, ncut, ifail;
int jreg, jang, jprf;

int cx = naxis_row/2;
int cy = naxis_col/2;
r  = (float)naxis_row/2;
int ng = naxis_row/16;//4; 
int nb_lines=3;
float *Dx, *xx;
#ifdef HMI
	float *line;
	line = (float *) malloc(sizeof(float)*(nb_lines*ng));
#elif  MDI
	short *line;
	line = (short *) malloc(sizeof(short)*(nb_lines*ng));
#endif
	
if (skipc == 1) 
{
	/* Initial guess estimate of center position from Richard's code*/
		cmx=(float)info->ix;
		cmy=(float)info->iy;
		r=(float)info->ir;
		guess_cx=cmx;
		guess_cy=cmy;
		guess_r=r;
}
else
{
		/***********************************************************************
		   Find the guess for the center position and radius                 
		***********************************************************************/
		int nb_iters=3;
		
		
			if(!line) return ERR_MALLOC_FAILED;
		
			for (i=0;i<naxis_row*naxis_col;i++)
				if (data[i] <0.0)
					data[i] =0.0;
		
			Dx = (float *) malloc(sizeof(float)*(nb_lines*ng));
				if(!Dx) return ERR_MALLOC_FAILED;
			xx = (float *) malloc(sizeof(float)*(nb_lines*ng));
				if(!xx) return ERR_MALLOC_FAILED;
			
			int k, kk, dh, vv; 
			float x1, x2, x3, y1, y2, y3, ma, mb, sav_max=0., mx[nb_lines], maximo;
			
			for(jj = 0; jj < nb_iters ; jj++) 
			{		
				for(ii = 0; ii < ng ; ii++)
				{
					xx[ii] = (float)ii;
					line[ii] = 0;
					if (data[ii*naxis_col+cy] > 0) line[ii] = data[ii*naxis_col+cy];
				}
		
				for(ii = naxis_row-1; ii > naxis_row-ng-1 ; ii--) 
				{
					vv=abs(ii-naxis_row+1);
					xx[ng+vv] = (float)ii;
					line[ng+vv] = 0;
					if (data[ii*naxis_col+cy] > 0) line[ng+vv] = data[ii*naxis_col+cy];
				 }
				
				for(ii = 0; ii < ng ; ii++) 
				{
					xx[(2*ng)+ii] =  (float)ii;
					line[(2*ng)+ii] = 0;
					if (data[cx*naxis_col+ii] > 0) line[(2*ng)+ii] = data[cx*naxis_col+ii];
				}
				dh = 1;
				for(k = 0; k < nb_lines ; k++) 
				{
					kk=k*ng;
					 /* Calculate the Derivative */ 
					Dx[kk]		=(float)(-3*line[kk+4]+16*line[kk+3]-36*line[kk+2]+48*line[kk+1]-25*line[kk])/(12*dh);
					Dx[kk+1]	=(float)(line[kk+4]-6*line[kk+3]+18*line[kk+2]-10*line[kk+1]-3*line[kk])/(12*dh);
					for (i=2; i< ng-2; i++) 
						Dx[kk+i]=(float)(line[kk+i-2]-8*line[kk+i-1]+8*line[kk+i+1]-line[kk+i+2])/(12*dh);
					Dx[kk+ng-2] =(float)(-line[kk+ng-5] +6*line[kk+ng-4]-18*line[kk+ng-3]+10*line[kk+ng-2]+3*line[kk+ng-1])/(12*dh);
					Dx[kk+ng-1] =(float)(3*line[kk+ng-5]-16*line[kk+ng-4]+36*line[kk+ng-3]-48*line[kk+ng-2]+25*line[kk+ng-1])/(12*dh);
					
					/* square the result */
					for(ii = 0; ii < ng ; ii++) 
						Dx[kk+ii]=Dx[kk+ii]*Dx[kk+ii];
		
					for(ii = 0; ii < ng ; ii++) 
						if (line[kk+ii] < 0.0001)   // IS: because a float cannot be compared to 0
							Dx[kk+ii]=0; 
		
					 /* smooth */
					for(ii = 1; ii < ng-1 ; ii++) 
						Dx[kk+ii]=(Dx[kk+ii-1]+Dx[kk+ii]+Dx[kk+ii+1])/3;
					
					 /* find the maximum */
					mx[k]=0.;					
					maximo=-1.;
					for(ii = 1; ii < ng-1 ; ii++) {
						if (Dx[kk+ii] > maximo) {
							maximo=Dx[kk+ii];
							mx[k]=xx[kk+ii];
						}
					}
				//fprintf(info->fd,"\n -jj %ld k %d ii %ld kk %d mx[k] %d\n",jj,k,ii,kk,mx[k]);
				} // endfor k
		
				/* find the geometric center */
				
				y1 = (float)mx[0];
				x1 = (float)cy;
				y3 = (float)mx[1];
				x3 = (float)cy;
				y2 = (float)cx;
				x2 = (float)mx[2];
				//fprintf(info->fd,"ici : %f %f %f %f %f %f\n",y1,x1,y3,x3,y2,x2);
				// for tuning purpose only
				if (mx[0] > sav_max) sav_max=mx[0];
				if (mx[1] > sav_max) sav_max=mx[1];
				if (mx[2] > sav_max) sav_max=mx[2];
				
				ma=(y2-y1)/((x2-x1)*1);
				mb=(y3-y2)/((x3-x2)*1);
				//fprintf(info->fd,"ici2 : %f %f \n",ma,mb);
				
				cx=(int) ((ma*mb*(y1-y3)+mb*(x1+x2)-ma*(x2+x3))/(2*(mb-ma)));
				cy=(int) (-(cx-(x1+x2)/2.)/ma+(y1+y2)/2.);
				r=(float)sqrt((x1-cx)*(x1-cx)+(y1-cy)*(y1-cy));		
				//fprintf(info->fd,"ici3 : %i %i %f \n",cx,cy,r);
			} // endfor jj
		
		guess_cx=(float)cx;
		guess_cy=(float)cy;
		guess_r=r;
		
		if (debug)
		{
			fprintf(info->fd,"Guess : cx = %6i, cy = %6i, r = %6.2f\n", cx, cy, r);
		}
	
	/************************************************************************/
	/*                        set parameters                               */
	/************************************************************************/		
	/* Initial guess estimate of center position */
		cmx=(float)cx;
		cmy=(float)cy;
} //skipc

jreg = nreg ;
jang = nang + 1;
jprf = nprf;

/************************************************************************/
/*                        Make annulus data                             */
/************************************************************************/
float d, *x, *y, *v;

x = (float *) malloc(sizeof(float) * S);
	if(!x) return ERR_MALLOC_FAILED;
y = (float *) malloc(sizeof(float) * S);
	if(!y)return ERR_MALLOC_FAILED;
v = (float *) malloc(sizeof(float) * S);
	if(!v) return ERR_MALLOC_FAILED;

	/* Select Points */
	jk=-1;
	for(ii = 0; ii < naxis_row; ii++)
	{
		for(jj = 0; jj < naxis_col; jj++) 
		{ 
			d=(float)sqrt((double)(ii-cmy)*(ii-cmy)+(jj-cmx)*(jj-cmx));
			if (d<=(r+w/2.))
			{
				if (d>=(r-w/2.)) 
				{
					 jk=jk+1;
					 x[jk]=(float)jj;
					 y[jk]=(float)ii;
					 v[jk]=data[ii*naxis_col+jj];
				}
			}
		}
	}
	
	if (debug)
	{
		fprintf(info->fd,"jk = %ld\n", jk);
	}
		printf("jk = %ld\n", jk);
	if (jk >= S) return ERR_SIZE_ANN_TOO_BIG;

/************************************************************************/
/*                  Call Fortran Code Subrotine limb.f                  */
/************************************************************************/
float *anls, *rprf, *lprf, *alph, *beta, *b0, *beta1, *beta2, *beta3;
int nbc=3;
anls = (float *) malloc(sizeof(float)*(nbc*jk));
	if(!anls) return ERR_MALLOC_FAILED;
rprf = (float *) malloc(sizeof(float)*(nprf));
	if(!rprf) return ERR_MALLOC_FAILED;
lprf = (float *) malloc(sizeof(float)*((nang+1)*nprf));
	if(!lprf) return ERR_MALLOC_FAILED;
alph = (float *) malloc(sizeof(float)*(nreg));
	if(!alph) return ERR_MALLOC_FAILED;
beta = (float *) malloc(sizeof(float)*(nreg));
	if(!beta) return ERR_MALLOC_FAILED;
beta1 = (float *) malloc(sizeof(float)*(nreg));
	if(!beta1) return ERR_MALLOC_FAILED;
beta2 = (float *) malloc(sizeof(float)*(nreg));
	if(!beta2) return ERR_MALLOC_FAILED;
beta3 = (float *) malloc(sizeof(float)*(nreg));
	if(!beta3) return ERR_MALLOC_FAILED;
b0 = (float *) malloc(sizeof(float)*(nreg));
	if(!b0) return ERR_MALLOC_FAILED;

	for (ii = 0; ii < jk; ii++) 
	{
		anls[nbc*ii]=x[ii];
		anls[nbc*ii+1]=y[ii];
		anls[nbc*ii+2]=v[ii];
	}
	for (ii=0; ii <nreg; ii++) beta[ii]=0.;

	// fortran call:
	if (debug) fprintf(info->fd,"entering limb_-----------\n");
	int centyp;
	//#1
	if (skipc == 1) centyp=1; else centyp=0;
	limb_(&anls[0],&jk, &cmx, &cmy, &r, &nitr, &ncut, &nang, &nprf, &rprf[0], &lprf[0], &nreg, &rsi, &rso, 
		&dx, &dy, &jreg, &jang, &jprf, &alph[0], &b0[0], &ifail, &beta[0], &centyp); 
		//&dx, &dy, &jreg, &jang, &jprf, &alph[0], &beta[0], &ifail, &b0[0], &centyp); 
	if (debug) fprintf(info->fd,"exiting limb--------%d---\n",ifail);
	centyp=1;
	if (iter > 1 && ifail == 0)
	{	
		if (debug) fprintf(info->fd,"re-entering-----------\n");
		//#2
		limb_(&anls[0],&jk, &cmx, &cmy, &r, &nitr, &ncut, &nang, &nprf, &rprf[0], &lprf[0], &nreg, &rsi, &rso, 
			&dx, &dy, &jreg, &jang, &jprf, &alph[0], &beta1[0], &ifail, &b0[0], &centyp); 
		if (debug) fprintf(info->fd,"exiting limb_2nd time--------%d---\n",ifail);
		for (ii=0; ii<nreg; ii++) b0[ii]=b0[ii]+beta1[ii];// <<< TO CONVERT THIS WITH POINTERS
		//#3
		if(iter > 2 && ifail ==0)
		{	
			limb_(&anls[0],&jk, &cmx, &cmy, &r, &nitr, &ncut, &nang, &nprf, &rprf[0], &lprf[0], &nreg, &rsi, &rso, 
				&dx, &dy, &jreg, &jang, &jprf, &alph[0], &beta2[0], &ifail, &b0[0], &centyp); 
			if (debug) fprintf(info->fd,"exiting limb_3rd time--------%d---\n",ifail);
			for (ii=0; ii<nreg; ii++) b0[ii]=b0[ii]+beta2[ii];// <<< TO CONVERT THIS WITH POINTERS
			if ((iter > 3) && (ifail ==0))
			{	
				//#4
				for (ii=0; ii<nreg; ii++) b0[ii]=beta[ii]+beta1[ii];// <<< TO CONVERT THIS WITH POINTERS
				limb_(&anls[0],&jk, &cmx, &cmy, &r, &nitr, &ncut, &nang, &nprf, &rprf[0], &lprf[0], &nreg, &rsi, &rso, 
					&dx, &dy, &jreg, &jang, &jprf, &alph[0], &beta3[0], &ifail, &b0[0], &centyp); 
				if (debug) fprintf(info->fd,"exiting limb_4th time--------%d---\n",ifail);
				for (ii=0; ii <nreg; ii++) b0[ii]=b0[ii]+beta3[ii];// <<< TO CONVERT THIS WITH POINTERS
			}
		} //ici faire qqch si ifail >0 sauver ldf,a,b de la version d'avant...
	}

	if (debug)
	{
		fprintf(info->fd,"ifail = %6d\n", ifail);
		fprintf(info->fd,"nitr = %6d\n", nitr);
		fprintf(info->fd,"cx = %6.2f, cy = %6.2f\n", cmx, cmy);
	}
	
	if(ifail == 0) // || ifail == 7) 
	{		
		/************************************************************************/
		/*          Compute the mean of the center of the image                 */
		/************************************************************************/
		double cmean,ctot=0.0;
		long nbp=0;
		float limx_m=cmx-250;
		float limx_p=cmx+250;
		float limy_m=cmy-250;
		float limy_p=cmy+250;
		for (i=limx_m;i<limx_p;i++)
		{
			for (j=limy_m;j<limy_p;j++)
			{
				ctot=ctot+data[i*naxis_col+j]; 
				nbp++;
			}
		}
		cmean=ctot/nbp;
		if (debug)
		{
			fprintf(info->fd,"cmean = %6.4f (ctot= %6.4f , nbp=%ld)\n", cmean,ctot,nbp); //comparer avec idl
		}

		/************************************************************************/
		/*                  Compute the Inflection Point                      */
		/************************************************************************/
		float *LDF, *D;
		LDF = (float *) malloc(sizeof(float)*(nprf));
			if(!LDF) return ERR_MALLOC_FAILED;
		D = (float *) malloc(sizeof(float)*(nprf));
			if(!D) return ERR_MALLOC_FAILED;

		double ip, ip1, maxim; 
		double radius = {0.0};
		double A[6] = { 0., 0., 0., 0., 0., 0. };
		double erro[6] = { 0., 0., 0., 0., 0., 0. };
		int cont, c, ret;
		float h;
			
		// for saving them in FITS file
		float *save_ldf, *save_alpha, *save_beta, *save_fulldf; //,*save_beta1,*save_beta2;
		double *save_as, *save_es, *save_radius, *save_ip;
		long aeris_nrow=nang+1, aes_ncol=6;
		long ldf_nrow=nang+2, ldf_ncol=nprf;  //ii -nbr of ldfs, jj -nbr of points for each ldf
		long ab_nrow=nreg, ab_ncol=2;
		long fulldf_nrow=nang+2;
		long fulldf_ncol=nprf ;
	    save_as  	= (double *) malloc((aeris_nrow*aes_ncol)*sizeof(double));  
			if(!save_as) return ERR_MALLOC_FAILED;
	    save_es  	= (double *) malloc((aeris_nrow*aes_ncol)*sizeof(double));
			if(!save_es) return ERR_MALLOC_FAILED;
	    save_radius = (double *) malloc(aeris_nrow*sizeof(double));
			if(!save_radius) return ERR_MALLOC_FAILED;
	    save_ip		= (double *) malloc(aeris_nrow*sizeof(double));
			if(!save_ip) return ERR_MALLOC_FAILED;
	    save_ldf 	= (float *) malloc((ldf_nrow*ldf_ncol)*sizeof(float));
			if(!save_ldf) return ERR_MALLOC_FAILED;
	    save_alpha  = (float *) malloc((ab_nrow)*sizeof(float));
			if(!save_alpha) return ERR_MALLOC_FAILED;
	    save_beta   = (float *) malloc((ab_nrow)*sizeof(float));
			if(!save_beta) return ERR_MALLOC_FAILED;
		//for later use
	    save_fulldf 	= (float *) malloc((fulldf_nrow*fulldf_ncol)*sizeof(float));
			if(!save_fulldf) return ERR_MALLOC_FAILED;
			
		for (cont=0; cont<=nang; cont++)
		{	 
			 ip1=0.;
			 /* dx */
			 h=(float)rprf[1]-rprf[0];
			
			 /* Take the last (mean) LDF from lprf vector */
			 for(ii = 0; ii < nprf ; ii++) 
			 	//LDF[ii]=lprf[(cont*(nang+1))+ii];
			 	LDF[ii]=lprf[(nprf*cont)+ii];
			/*
			if (cont==180)
			{
				FILE *opf1;
		      	opf1 = fopen("fout1.txt", "w");      
				for(ii = 0; ii < nprf ; ii++) 
				{
				 	LDF[ii]=lprf[(nprf*cont)+ii];
					fprintf(opf1,"%f\n",LDF[ii]);
				}
				fclose(opf1);
			}
			*/
			 /* Calculate the Derivative */
			 D[0]=(-3*LDF[4]+16*LDF[3]-36*LDF[2]+48*LDF[1]-25*LDF[0])/(12*h);
			 D[1]=(LDF[4]-6*LDF[3]+18*LDF[2]-10*LDF[1]-3*LDF[0])/(12*h);
			 for (i=2; i< nprf-2; i++) 
				D[i]=(LDF[i-2]-8*LDF[i-1]+8*LDF[i+1]-LDF[i+2])/(12*h);
			 D[nprf-2]=(-LDF[nprf-5]+6*LDF[nprf-4]-18*LDF[nprf-3]+10*LDF[nprf-2]+3*LDF[nprf-1])/(12*h);
			 D[nprf-1]=(3*LDF[nprf-5]-16*LDF[nprf-4]+36*LDF[nprf-3]-48*LDF[nprf-2]+25*LDF[nprf-1])/(12*h);
			
			 /* square the result */
			 for(ii = 0; ii < nprf ; ii++) 
				 D[ii]=D[ii]*D[ii];
			
			/* smooth */
			 for(ii = 1; ii < nprf-1 ; ii++) 
				 D[ii]=(D[ii-1]+D[ii]+D[ii+1])/3;
			
			 /* find the maximum */
			 jj=0;
			 maxim=-1;
			 for(ii = 1; ii < nprf-1 ; ii++) 
			 {
			 	if (D[ii] > maxim) 
			 	{
			  		maxim=(double)D[ii];
			   		jj=ii;
			  	}
			 }
			
			/* improve the maximum estimation looking at the 2 neibors points */
			ip=rprf[jj]-h/2.*(D[jj-1]-D[jj+1])/(2*D[jj]-D[jj-1]-D[jj+1]);
			
			if (debug)
			{
				if (cont==nang)
					fprintf(info->fd,"Inflection Point 1: %ld %d %8.5f %8.5f\n", jj, nprf, rprf[jj], ip);
			}
			ip1=ip;
			/* FIT A GAUSSIAN PLUS QUADRATIC FUNCTION */
			
			/* Select Region */
			int m_in, m_ex, N;
			 /* Gaussian + quadratic function */
			 /* After Launch verify if this number is OK */
			 /* Should contain the pixels around the 3*FWHM level */
			m_in=jj-r_size;
			m_ex=jj+r_size;
			
			if (m_in < 0) m_in = 0;
			if (m_ex > nprf-1) m_ex =nprf-1;
			N=m_ex-m_in+1;
			/* parameters:- initial guess */
			A[0]= maxim;
			A[1]= ip;
			A[2]= 0.5;
			A[3]= 12*maxim;
			A[4]= 0;
			A[5]= 0.;
			for (c=0;c<6;c++) erro[c]=0.;
		
			if (N >= 6) 	// degree of freedom 6 parameters, 6 values minimum
			{
				double px[N], der[N] , sigma[N];
				
				jj=-1;
				for(ii = m_in; ii <= m_ex ; ii++)
				{
					   jj++;
					   px[jj]=rprf[ii];
					   der[jj]=D[ii];
					   sigma[jj] = 1.;
			 	}
				A[4]= der[N-1]-der[0];
				if (debug==2) fprintf(info->fd,"#: %d\n",cont);
				/*
				int is;
				fprintf(info->fd,"\n--> %d: %d %d %d %ld\n A: ",cont,N,m_in,m_ex,jj);
				for (is=0;is<6;is++) fprintf(info->fd,"%6.3f ",A[is]);
				fprintf(info->fd,"\n px: ");
				for (is=0;is<N;is++) fprintf(info->fd,"%6.3f ",px[is]);
				fprintf(info->fd,"\n erro: ");
				for (is=0;is<6;is++) fprintf(info->fd,"%6.3f ",erro[is]);
				fprintf(info->fd,"\n sigma: ");
				for (is=0;is<N;is++) fprintf(info->fd,"%6.3f ",sigma[is]);
				fprintf(info->fd,"\n der: ");
				for (is=0;is<N;is++) fprintf(info->fd,"%6.3f ",der[is]);
				fprintf(info->fd,"\n");
				if (cont >= 179) debug=3;
				*/
				ret=gaussfit(der, px, sigma, A, erro, N, debug,info->fd);
				if (ret < 0)
				{
					if (debug) fprintf(info->fd,"in error #1 %d err:%d\n",cont,ret);
					for (c=0;c<6;c++) 
					{	
						A[c]=0.;
						erro[c]=0.;
					}
					radius=0.;
				}
				else
				{			
				 	/* FIND THE MAXIMUM OF THE GAUSSIAN PLUS QUADRATIC FUNCTION */
				 	radius = A[1]; /* Initial Guess */
				 	radius = fin_min(A, radius, grange, debug,info->fd);
					if (radius < 0)
					{
						if (debug) fprintf(info->fd,"in error #2 %d\n",cont);
						for (c=0;c<6;c++) 
						{	
							A[c]=0.;
							erro[c]=0.;
						}
						radius=0.;
					}
				}
			} 
			else 
			{ 
				if (debug)
				{	
					fprintf(info->fd,"Inflection point too close to annulus data border\n");
				}
				for (c=0;c<6;c++) erro[c]=0.;
				radius=ip;
				save_ip[cont]=0.;
			}
			// save them
			for (c=0;c<aes_ncol;c++)
			{
				save_as[c*aeris_nrow+cont]=A[c];
				save_es[c*aeris_nrow+cont]=erro[c];
			}
			save_radius[cont]=radius;
			save_ip[cont]=ip1;

		} // endfor-cont
		if (debug)
		{	
			fprintf(info->fd,"Inflection Point 2: %8.5f %8.5f %8.5f\n", A[1], erro[1], radius);
			fprintf(info->fd,"-----: %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n ",A[0],A[1],A[2],A[3],A[4],A[5]);
			fprintf(info->fd,"-----: %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n ",erro[0],erro[1],erro[2],erro[3],erro[4],erro[5]);
		}

	//**********************************************************************
	//						Save LDFs & AB data in a FITS file                      
	//**********************************************************************
	
	// full LDF
		float *p_sfldf=&save_fulldf[0];
		float *pl_sfldf=&save_fulldf[(fulldf_nrow*fulldf_ncol)-1];
		//right now just initialize it, later move its values
		while (p_sfldf <= pl_sfldf)
			*(p_sfldf++)=0;
	
	// LDF 				// lprf & rprf come from limbfit.f			
		float *p_sldf=&save_ldf[0];
		float *pl_sldf=&save_ldf[(ldf_nrow*ldf_ncol)-1];
		float *p_rprf=&rprf[0];
		float *pl_rprf=&rprf[nprf-1];
		while (p_rprf <= pl_rprf)
			*(p_sldf++)=*(p_rprf++);

		float *p_lprf=&lprf[0];
		p_sldf=&save_ldf[ldf_ncol];
		while (p_sldf <= pl_sldf)
			*(p_sldf++)=*(p_lprf++);

	// AB 
		float *p_alph=&alph[0];
		float *pl_alph=&alph[ab_nrow-1];
		float *p_beta=&b0[0]; // sum of beta, beta1, beta2
		float *p_save_alpha=&save_alpha[0];
		float *p_save_beta=&save_beta[0];
		while (p_alph <= pl_alph)
		{
			*(p_save_alpha++)=*(p_alph++);
			*(p_save_beta++)=*(p_beta++);
		}

		// Update Returned Structure when process succeeded                    
		results->cenx=cmx;
		results->ceny=cmy;
		results->radius=radius;
		results->cmean=cmean;
		results->error=0;		
		results->numext=5;		
		results->fits_ext0_naxis=2;
		results->fits_ext0_naxis1=ldf_ncol;	
		results->fits_ext0_naxis2=ldf_nrow;
		results->fits_ext1_nrows=aeris_nrow;
		results->fits_ext1_tfields=aes_ncol;
		results->fits_ext2_nrows=aeris_nrow;
		results->fits_ext2_tfields=aes_ncol;
		results->fits_ext3_nrows=aeris_nrow;
		results->fits_ext3_tfields=1;
		results->fits_ext4_nrows=ab_nrow;
		results->fits_ext4_tfields=ab_ncol;
		results->fits_ext5_nrows=aeris_nrow;
		results->fits_ext5_tfields=1;
		results->fits_ext6_nrows=fulldf_nrow;
		results->fits_ext6_tfields=fulldf_ncol;
		if (ifail == 0) 
			results->quality=9; 
		else //IS to be tuned but in function of what?
			results->quality=5; 
		results->fits_ldfs_data=save_ldf; 
		results->fits_fulldfs=save_fulldf; 
		results->fits_as=save_as; 
		results->fits_es=save_es; 
		results->fits_r=save_radius; 
		results->fits_ip=save_ip; 
		results->fits_ab_dataa=save_alpha; 
		results->fits_ab_datab=save_beta; 
		results->gxc=guess_cx;
		results->gyc=guess_cy;
		results->gr=guess_r;
		//results->max_limb=sav_max;

		if (debug == 3)
		{
			results->annulus=anls;
			results->annulus_jk=jk;
		}
		else 
		{
			results->annulus=0;
			results->annulus_jk=0;
		}
		free(D);
		free(LDF);
	} 
	else 
	{
		if (debug)
		{	
			fprintf(info->fd,"limb.f routine returned ifail = %2d\n", ifail);
 		}     	
		// Update Returned Structure when process failed                    
		results->cenx=guess_cx;
		results->ceny=guess_cy;
		results->radius=guess_r;
		results->fits_ext0_naxis=0;
		results->numext=0;		
		results->error=ifail;
		results->quality=0;
		results->gxc=guess_cx;
		results->gyc=guess_cy;
		results->gr=guess_r;
		results->max_limb=0.;
      	ret_code=ERR_LIMBFIT_FAILED;

		if (debug == 3)
		{
			results->annulus=anls;
			results->annulus_jk=jk;
		}
		else 
		{
			results->annulus=0;
			results->annulus_jk=0;
		}
	}

	// Write common keywords
	time_t t = time(NULL);
	char tdate[50];
    struct tm *timeptr;
	//  Fri Feb 12 16:13:59 2010 -> '2009.10.26_04:00:30_TAI'   
	timeptr=gmtime(&t);
	sprintf(tdate,"%4d.%d.%d_%d:%d:%2d_TAI\0",
		1900+timeptr->tm_year,timeptr->tm_mon,	timeptr->tm_mday,
		timeptr->tm_hour,timeptr->tm_min,timeptr->tm_sec);
	strcpy(results->proc_date,tdate);
	strcpy(results->code_name,CODE_NAME);
	strcpy(results->code_version,CODE_VERSION);		
printf("ici ok\n");
	// IS: do not free 'data' et others passed from or to the structure !
	free(x);
	free(y);
	free(v); 
	if (skipc == 0) 
	{
		free(line);
		free(Dx);
		free(xx);
	}
	free(rprf);
	free(lprf);
	free(alph);
	free(beta);
	free(beta1);
	free(beta2);
	free(beta3);
	free(b0);
	if (debug != 3) free(anls);

printf("ici ok\n");
	
	if (debug)
	{	
		fprintf(info->fd,">>>>end of limbfit with: %d\n",ret_code);
	}
	
return(ret_code);

} /* end of main */


/*--------------------------------------------------------------------------*/
int gaussfit(double y[], double t[],double sigma[], double A[6], double erro[6], long N, int debug, FILE *fd)
/* Calculate a Least SqrareGaussian + Quadratic fit              */
/*      Uses the GNU Scientific Library                          */
/* Marcelo Emilio (c) v 1.0 Jan 2009                             */
/* fits A[0]exp(-z^2/2)+A[3]+A[4]*x+A[5]*x^2                     */
/* z=(t-A[1])/A[2]                                               */
/* Need to add :                                                 */
/* #include <gsl/gsl_rng.h>                                      */
/* #include <gsl/gsl_randist.h>                                  */
/* #include <gsl/gsl_vector.h>                                   */
/* #include <gsl/gsl_blas.h>                                     */
/* #include <gsl/gsl_multifit_nlin.h>                            */
/* #include "expfit.c"                                           */
/* compiles as  gcc program.c --lm -lgsl -lgslcblas              */
/* y --> f(x) - ordinate values                                  */
/* t --> x(i) - abscissa values                                  */
/* sigma --> independent gaussian erros                          */
/* N --> Number of points in the vector                          */
/* A --> Input a initial guess and output the  result            */
/* erro --> Chi-square of A                                      */

{
	int ret_code=0;
	
	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;
	int status;
	unsigned int i, iter = 0;
	const size_t n = N;
	const size_t p = 6;
		 
	gsl_matrix *covar = gsl_matrix_alloc (p, p);
	/* double t[N], y[N], sigma[N]; */
	struct dataF d = { n, t, y, sigma};
	gsl_multifit_function_fdf f;
	
	/* Initial Guess */
	double x_init[6];
	for (i=0; i <=5; i++) 
	   x_init[i]=A[i];
	
	gsl_vector_view x = gsl_vector_view_array (x_init, p);
	const gsl_rng_type * type;
	gsl_rng * r;
		 
	gsl_rng_env_setup();
		 
	type = gsl_rng_default;
	r = gsl_rng_alloc (type);
		 
	f.f = &expb_f;
	f.df = &expb_df;
	f.fdf = &expb_fdf;
	f.n = n;
	f.p = p;
	f.params = &d;
		 
	//double d1,d2,d3,d4;

	T = gsl_multifit_fdfsolver_lmsder;
	s = gsl_multifit_fdfsolver_alloc (T, n, p);
	gsl_multifit_fdfsolver_set (s, &f, &x.vector);  

	do
	{
	
		iter++;
		status = gsl_multifit_fdfsolver_iterate (s);     
		if (!status)	
			status = gsl_multifit_test_delta (s->dx, s->x, 1e-6, 1e-2);  		/* IS: change epsrel (last arg of gsl_multifit_test_delta) after the SDO Launch */ 
		else 
			fprintf (fd,"iter: %u, status: %d (%s)\n",iter,status,gsl_strerror (status));
	}
	while (status == GSL_CONTINUE && iter < 500);
			
	if (status >=0)  
	{
		gsl_multifit_covar (s->J, 0.0, covar);
			 
		#define FIT(i) gsl_vector_get(s->x, i)
		#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
			 
		double chi = gsl_blas_dnrm2(s->f);
		double dof = (double)(n - p);
		double c = GSL_MAX_DBL(1, chi / sqrt(dof)); 
		  
		
		for (i=0; i <=5; i++) {
			A[i]=FIT(i);
			erro[i]=c*ERR(i);
		}  
		/*
		if (debug==3)
		{
			int is;
		  fprintf(fd,"\n chisq/dof = %6.3f %6.3f %6.3f\n",  chi,dof,c);     
			for (is=0;is<6;is++) fprintf(fd,"%6.3f ",A[is]);
			fprintf(fd,"\n erro: ");
			for (is=0;is<6;is++) fprintf(fd,"%6.3f ",erro[is]);
			fprintf(fd,"\n");
		}
		*/
	}
	else
	{
		fprintf (fd,"(gsl_multifit_fdfsolver_iterate) failed, (iter=%u) gsl_errno=%d\n",iter,status); ///to change to a debug printf
		ret_code=ERR_GSL_GAUSSFIT_FDFSOLVER_FAILED;			
	}

	
	if (debug==2)
	{																		
		fprintf (fd,"Nonlinear LSQ fitting status = %s (%d)\n", gsl_strerror (status),status); 
	}
	gsl_multifit_fdfsolver_free (s);
	gsl_matrix_free (covar);
	gsl_rng_free (r);
	return(ret_code);
}


/*--------------------------------------------------------------------------*/

double fin_min(double A[6], double m, int range, int debug, FILE *fd)
/* Calculate the maximum of the quadratic + gaussian function    */
/*                 using Brent algorithm                         */
/* Marcelo Emilio (c) v 1.0 Jan 2009                             */
/* A are the parameters of the gaussian                          */
/* m is the guess where the function is maximum                  */
/* Calculate the maximum around                                  */
/*                 (-5 + m) < m > (5 - m)  ( in pixel units )    */
/* Returns the m vaule where the function is maximum             */
/* requires: #include "expmax.c" --> definition of the function  */
/*           #include <gsl/gsl_errno.h>                          */
/*           #include <gsl/gsl_math.h>                           */
/*           #include <gsl/gsl_min.h>                            */
{
	int status;
	gsl_set_error_handler_off();
	
	 int iter = 0, max_iter = 1000;
	 const gsl_min_fminimizer_type *T;
	 gsl_min_fminimizer *s;
	 
	 gsl_function F;
	
	//double m_expected = m+0.01; //1e-8;
	//double a = m-5000, b = m+5000;
	double a = m-range, b = m+range;      // initially = 2
	struct exp_plus_quadratic_params params= { A[0], A[1], A[2], A[3], A[4], A[5] }; 
	
	  F.function = &exp_plus_quadratic_function;
	  F.params = &params;
		 
	  T = gsl_min_fminimizer_brent;
	  s = gsl_min_fminimizer_alloc (T);
	
	  status=gsl_min_fminimizer_set (s, &F, m, a, b);     
	  if (status)
	  {
			if (status == GSL_EINVAL) {
				fprintf (fd,"(gsl_min_fminimizer_set) doesn't find a minimum, m=%f\n", m );
			} else {
				fprintf (fd,"(gsl_min_fminimizer_set) failed, gsl_errno=%d\n",status);
			}
			return ERR_GSL_FIN_MIN_SET_FAILED;
	  }
	  
	  do
	  {
	   iter++;
	   status = gsl_min_fminimizer_iterate (s);
		 
	   m = gsl_min_fminimizer_x_minimum (s);
	   a = gsl_min_fminimizer_x_lower (s);
	   b = gsl_min_fminimizer_x_upper (s);
		 
	   status = gsl_min_test_interval (a, b, 1E-4, 0.0); 
	  }
	  while (status == GSL_CONTINUE && iter < max_iter);
	  
	  if (status)		// regarder lequel plantait et verifier si on a autre chose que GSL_EINVAL comme erreur
		{				// ajouter debug pour les messages apres
			if (status == GSL_EINVAL) {
				fprintf (fd,"invalid argument, n=%8.5f\n", a);
			} else {
				fprintf (fd,"failed, gsl_errno=%d\n",status);
			}
			return ERR_GSL_FIN_MIN_PRO_FAILED; //IS a arranger! et a recuperer dans le main
		}
	  
	  if (debug==2)
	  {
			fprintf (fd," a:%8.5f b:%8.5f m:%8.5f iter=%d, b-a:%f\n",a,b,m,iter,b-a)     ;
	  }


	  gsl_min_fminimizer_free (s);
		 
	  return m;

}
