/* 
	I.Scholl 

		limbfit(): adapted from Solar Astrometry Program (solarSDO.c)
									Marcelo Emilio - Jan 2010 


	#define CODE_NAME 		"limbfit_tas"
	#define CODE_VERSION 	"V5.03" 
	#define CODE_DATE 		"Thu Apr  4 10:50:46 HST 2013" 
*/

#include "limbfit_tas.h"

#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

void sav_b0(float *pf_sb0, float *pl_sb0, float *pf_b0)
{
	float *p_sb0=pf_sb0;
	float *p_b0=pf_b0;
	while(p_sb0<=pl_sb0) *(p_sb0++)=*(p_b0++);
}

void sum_b0(float *beta, float *pf_b0, float *pl_b0)
{
	float *p_b0=pf_b0;
	float *p_beta=beta;
	while(p_b0 <= pl_b0) *(p_b0++)=*(p_b0)+*(p_beta++);				
}

int limbfit(LIMBFIT_INPUT *input, LIMBFIT_OUTPUT *results, LIMBFIT_IO_PUT *ios)
{
static char *log_msg_code="limbfit";
char log_msg[200];

if (results->debug) lf_logmsg("DEBUG", "APP", 0, 0, "", log_msg_code, results->opf);

/************************************************************************
					INIT VARS coming from build_lev1.c
************************************************************************/	
	float *data=input->data;
//	float MIN_VALUE=BAD_PIXEL_VALUE;
//	float MAX_VALUE=(float)fabs(BAD_PIXEL_VALUE);


int 	ret_code = 0;

int 	w		 = ANNULUS_WIDTH;
long 	S		 = MAX_SIZE_ANN_VARS;
int 	nang	 = NUM_LDF; //180
int 	nprf	 = NUM_RADIAL_BINS; //64
int 	nreg	 = NUM_AB_BINS;
float 	rsi		 = LO_LIMIT;
float 	rso		 = UP_LIMIT;
int 	r_size	 = NUM_FITPNTS;
int		grange	 = GUESS_RANGE;
float 	dx		 = INC_X;
float 	dy		 = INC_Y;
int		naxis_row= input->img_sz0;
int		naxis_col= input->img_sz1;
float	lahi	 = AHI;
int		fldf;
//int		skip	 = SKIPGC;
//float 	flag	 = BAD_PIXEL_VALUE; 
int		iter	 = NB_ITER;
if (input->iter!=NB_ITER) iter=input->iter; 
fldf=input->fldf; 

if (input->spe==1)
{
	iter=NB_ITER2;
	lahi=AHI2;
	rsi=LO_LIMIT2;
	rso=UP_LIMIT2;
}
	
/************************************************************************/
/*                        set parameters                               */
/************************************************************************/
long npixels=naxis_row*naxis_col;
long ii, jj, jk, i,j; 
float cmx, cmy,r;//, guess_cx, guess_cy, guess_r;
int nitr=0, ncut=0;

r  = (float)naxis_row/2;
	
/* Initial guess estimate of center position from Richard's code*/
	cmx=(float)input->ix;
	cmy=(float)input->iy;
	r=(float)input->ir;


/************************************************************************/
/*                        Make annulus data                             */
/************************************************************************/
	int nbc=3;
	float *p_anls=ios->pf_anls;		

	if (ios->is_firstobs == 0) 
	{
		ios->is_firstobs=1;
		float d;
		float w2p=(double)r+w/2.;
		float w2m=(double)r-w/2.;
		double iimcmy2,jjmcmx;
		/* Select Points */
		jk=-1;
		
		for(ii = 0; ii < naxis_row; ii++)
		{
			iimcmy2=(ii-cmy)*(ii-cmy);
			for(jj = 0; jj < naxis_col; jj++) 
			{ 
				jjmcmx=jj-cmx;
				d=(float)sqrt(iimcmy2+(jjmcmx)*(jjmcmx));
				if (d<=w2p && d>=w2m) 
				{
					jk++;
					*(p_anls++)=(float)jj;
					*(p_anls++)=(float)ii;
					*(p_anls++)=data[ii*naxis_col+jj];					 
				}
			}
		}
		ios->anls_nbpix=jk;
		ios->pl_anls=p_anls-1;
	
		if (results->debug)
		{
			sprintf(log_msg," jk = %ld", jk);
			lf_logmsg("DEBUG", "APP", 0, 0, log_msg, log_msg_code, results->opf);
		}
		if ((jk*3) >= S) 
		{
			lf_logmsg("ERROR", "APP", ERR_SIZE_ANN_TOO_BIG, 0,"nbc>S", log_msg_code, results->opf);
			return ERR_SIZE_ANN_TOO_BIG;
		}
	}
	else 
	{
		jk=ios->anls_nbpix;
		ii=0;
		p_anls=p_anls+2;
		while(p_anls<=ios->pl_anls) 
		{
			*(p_anls)=data[(int)ios->anls[nbc*ii+1]*naxis_col+(int)ios->anls[nbc*ii]];
			ii++;
			p_anls=p_anls+3;
		}		
	}

/************************************************************************/
/*                  Call Fortran Code Subrotine limb.f                  */
/************************************************************************/
	float *rprf, *lprf, *alph, *beta, *b0, *sb0; //, *beta1, *beta2, *beta3;
	long ab_nrow=nreg, ab_ncol=2;

	// contains the LDFs AXIS
	rprf = (float *) malloc(sizeof(float)*(nprf));
		if(!rprf) 
		{
			lf_logmsg("ERROR", "APP", ERR_MALLOC_FAILED, 0,"malloc failed (rprf)", log_msg_code, results->opf);
			return ERR_MALLOC_FAILED;
		}
	float *p_rprf, *pl_rprf;
	// contains LDFs inc. average ldf
	lprf = (float *) malloc(sizeof(float)*((nang+1)*nprf));
		if(!lprf) 
		{
			lf_logmsg("ERROR", "APP", ERR_MALLOC_FAILED, 0,"malloc failed (lprf)", log_msg_code, results->opf);
			return ERR_MALLOC_FAILED;
		}
	float *p_lprf, *pl_lprf;
	
	alph = (float *) malloc(sizeof(float)*(nreg));
		if(!alph) 
		{
			lf_logmsg("ERROR", "APP", ERR_MALLOC_FAILED, 0,"malloc failed (alph)", log_msg_code, results->opf);
			return ERR_MALLOC_FAILED;
		}
	float *pf_alph=&alph[0];
	float *p_alph;
	float *pl_alph=&alph[ab_nrow-1];
	
	beta = (float *) malloc(sizeof(float)*(nreg));
		if(!beta) 
		{
			lf_logmsg("ERROR", "APP", ERR_MALLOC_FAILED, 0,"malloc failed (xbeta)", log_msg_code, results->opf);
			return ERR_MALLOC_FAILED;
		}
	float *pf_beta=&beta[0];
	float *p_beta;
	
	b0 = (float *) malloc(sizeof(float)*(nreg));
		if(!b0) 
		{
			lf_logmsg("ERROR", "APP", ERR_MALLOC_FAILED, 0,"malloc failed (b0)", log_msg_code, results->opf);
			return ERR_MALLOC_FAILED;
		}
	float *pf_b0=&b0[0];
	float *p_b0;
	float *pl_b0=&b0[ab_nrow-1]; 
	
	/*
	sb0 = (float *) malloc(sizeof(float)*(nreg));
		if(!sb0) 
		{
			lf_logmsg("ERROR", "APP", ERR_MALLOC_FAILED, 0,"malloc failed (sb0)", log_msg_code, results->opf);
			return ERR_MALLOC_FAILED;
		}
	float *pf_sb0=&sb0[0];
	float *p_sb0;
	float *pl_sb0=&sb0[ab_nrow-1]; 
	*/
	
	p_b0=pf_b0;	
	while(p_b0<=pl_b0) *(p_b0++)=0.;

	// fortran call:
	if (results->debug) lf_logmsg("DEBUG", "APP", 0, 0, "entering limb", log_msg_code, results->opf);

	int centyp=0; //=0 do center calculation, =1 skip center calculation
	//#1
	if (input->cc == 0) centyp=1;
	
	//in limb_ call: beta contains the output beta, b0 is initialized with 0.0
	int it=1, ifail=0;	
	
	do
	{
		// init everything that is passed to fortran
		//
		p_alph=pf_alph;
		p_beta=pf_beta;
		while(p_alph<=pl_alph) 
		{
			*(p_beta++)=0.;
			*(p_alph++)=0.;
		}
		p_rprf=&rprf[0];
		pl_rprf=&rprf[nprf-1];
		while(p_rprf<=pl_rprf) *(p_rprf++)=0.;	
		p_lprf=&lprf[0];
		pl_lprf=&lprf[((nang+1)*nprf)-1];
		while(p_lprf<=pl_lprf) *(p_lprf++)=0.;

			if (results->debug)
			{
				sprintf(log_msg,"entering limb# %d", it);
				lf_logmsg("DEBUG", "APP", 0, 0, log_msg, log_msg_code, results->opf);
			}
		ifail=0;
		limb_(&ios->anls[0],&jk, &cmx, &cmy, &r, &nitr, &ncut, &rprf[0], &lprf[0],  &rsi, &rso, 
			&dx, &dy, &alph[0], &beta[0], &ifail, &b0[0], &centyp, &lahi); 
			if (results->debug)
			{
				sprintf(log_msg,"exiting limb ifail= %d", ifail);
				lf_logmsg("DEBUG", "APP", 0, 0, log_msg, log_msg_code, results->opf);
			}
		if(ifail==0)	
		{
			centyp=1;
			//if(it==iter && it>1) sav_b0(pf_sb0,pl_sb0,pf_b0);
			sum_b0(&beta[0],pf_b0,pl_b0);
		}
		it++;
	}
	while(it<=iter && ifail==0);



	if (results->debug)
	{
		sprintf(log_msg," nitr = %6d", nitr);
		lf_logmsg("DEBUG", "APP", 0, 0, log_msg, log_msg_code, results->opf);
		sprintf(log_msg," cmx = %6.2f, cmy = %6.2f", cmx, cmy);
		lf_logmsg("DEBUG", "APP", 0, 0, log_msg, log_msg_code, results->opf);
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
		if (results->debug)
		{
			sprintf(log_msg," cmean = %6.4f (ctot= %6.4f , nbp=%ld)", cmean,ctot,nbp);
			lf_logmsg("DEBUG", "APP", 0, 0, log_msg, log_msg_code, results->opf);
		}

		/************************************************************************/
		/*                  Compute the Inflection Point                      */
		/************************************************************************/
		float *LDF, *D, *t_ip, *p_t_ip;
		LDF = (float *) malloc(sizeof(float)*(nprf));
			if(!LDF) 
			{
				lf_logmsg("ERROR", "APP", ERR_MALLOC_FAILED, 0,"malloc failed (LDF)", log_msg_code, results->opf);
				return ERR_MALLOC_FAILED;
			}

		D = (float *) malloc(sizeof(float)*(nprf));
			if(!D)	
			{
				lf_logmsg("ERROR", "APP", ERR_MALLOC_FAILED, 0,"malloc failed (D)", log_msg_code, results->opf);
				return ERR_MALLOC_FAILED;
			}
	    t_ip 	= (float *) malloc((nang)*sizeof(float));
			if(!t_ip) 
			{
				lf_logmsg("ERROR", "APP", ERR_MALLOC_FAILED, 0,"malloc failed (save_ip)", log_msg_code, results->opf);
				return ERR_MALLOC_FAILED;
			}
		p_t_ip=&t_ip[0];

		double ip, maxim; 
		double radius = 0.;
		int cont, c, ret_gsl;
		float h;
			
		//if (iter > 1) ab_ncol=3; else ab_ncol=2;
		int degf=6;
		//if (iter > 2) degf=4; else degf=6; // problem with degf = 4...
		double A[degf];
		double erro[degf];


		 /* copy the last pixel of the mean ldf in the prev one to correct a fortran issue */
		lprf[(nprf*nang)+63]=lprf[(nprf*nang)+62];

		/* -----------------------------------------------------------------------------
		
			3 radius:
				* A[1]	 		= maximum of the Gaussian
				* radius 		= maximum of the Gaussian + quadratic fit -> R_LFS
				* ip, t_ip[*]	= inflection point 
		
		------------------------------------------------------------------------------- */			
			
		for (cont=0; cont<=nang; cont++)
		{ 
		 //cont=nang;
			ret_gsl=0;
			 /* dx */
			 h=(float)rprf[1]-rprf[0];
			
			 /* Take the last (mean) LDF from lprf vector */
			 for(ii = 0; ii < nprf ; ii++) LDF[ii]=lprf[(nprf*cont)+ii];

			 /* Calculate the Derivative */
			 D[0]=(-3*LDF[4]+16*LDF[3]-36*LDF[2]+48*LDF[1]-25*LDF[0])/(12*h);
			 D[1]=(LDF[4]-6*LDF[3]+18*LDF[2]-10*LDF[1]-3*LDF[0])/(12*h);
			 for (i=2; i< nprf-2; i++) D[i]=(LDF[i-2]-8*LDF[i-1]+8*LDF[i+1]-LDF[i+2])/(12*h);
			 D[nprf-2]=(-LDF[nprf-5]+6*LDF[nprf-4]-18*LDF[nprf-3]+10*LDF[nprf-2]+3*LDF[nprf-1])/(12*h);
			 D[nprf-1]=(3*LDF[nprf-5]-16*LDF[nprf-4]+36*LDF[nprf-3]-48*LDF[nprf-2]+25*LDF[nprf-1])/(12*h);

			 /* square the result */
			// for(ii = 0; ii < nprf ; ii++) 
			//	 D[ii]=D[ii]*D[ii];	//pointers here!
			float* p_D=&D[0];
			float* pl_D=&D[nprf-1];
			while(p_D<=pl_D) *(p_D++)=(*(p_D) * *(p_D));

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
			t_ip[cont]=rprf[jj]-h/2.*(D[jj-1]-D[jj+1])/(2*D[jj]-D[jj-1]-D[jj+1]);
			if (results->debug)
			{
				if (cont==nang)
				{
					sprintf(log_msg," Inflection Point 1: %ld %d %8.5f %8.5f", jj, nprf, rprf[jj], t_ip[cont]);
					lf_logmsg("DEBUG", "APP", 0, 0, log_msg, log_msg_code, results->opf);
				}
			}
		} // endfor-cont
		cont--;
		// now only the mean ldf
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
		A[1]= t_ip[cont];
		A[2]= 0.5;
		A[3]= 12*maxim;
		for (c=0;c<degf;c++) erro[c]=0.;
	
		if (N >= degf) 	// degree of freedom 6 parameters, 6 values minimum
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
			//A[4]= der[N-1]-der[0];
			//if (debug==2) fprintf(opf,"%s DEBUG_INFO in limbfit: #: %d\n",LOGMSG1,cont);

			ret_gsl=gaussfit(der, px, sigma, A, erro, N, degf,results->debug,results->opf);

			if (ret_gsl < 0)
			{
				if (cont==nang)
				{
					sprintf(log_msg," gaussfit failed for the averaged LDF %d err:%d", cont,ret_gsl);
					lf_logmsg("ERROR", "APP", 0, ret_gsl, log_msg, log_msg_code, results->opf);
				}
				for (c=0;c<degf;c++) 
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
				radius = fin_min(A, radius, grange, degf, results->debug,results->opf);
				if (radius < 0)
				{
					ret_gsl=(int)radius;
					if (cont==nang)
					{
						sprintf(log_msg," fin_min failed for the averaged LDF %d err:%d", cont, ret_gsl);
						lf_logmsg("ERROR", "APP", 0, ret_gsl, log_msg, log_msg_code, results->opf);
					}
					for (c=0;c<degf;c++) 
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
			sprintf(log_msg," Inflection point too close to annulus data border %d", cont);
			lf_logmsg("WARNING", "APP", 0, 0, log_msg, log_msg_code, results->opf);
			for (c=0;c<degf;c++) erro[c]=0.;
			// in this case: A[1]=radius=t_ip[last] 
			radius=t_ip[cont];
		}
		// save them
		for (c=0;c<degf;c++)
		{
			results->fits_as[c]=(float)A[c];
			results->fits_es[c]=(float)erro[c];
		}
		if (results->debug)
		{	
			sprintf(log_msg," Inflection Point 2: %8.5f %8.5f %8.5f", A[1], erro[1], radius);
			lf_logmsg("DEBUG", "APP", 0, 0, log_msg, log_msg_code, results->opf);
			if (degf==4)
			{
				sprintf(log_msg," -----: %8.5f %8.5f %8.5f %8.5f ", A[0],A[1],A[2],A[3]);
				lf_logmsg("DEBUG", "APP", 0, 0, log_msg, log_msg_code, results->opf);
				sprintf(log_msg," -----: %8.5f %8.5f %8.5f %8.5f ", erro[0],erro[1],erro[2],erro[3]);
				lf_logmsg("DEBUG", "APP", 0, 0, log_msg, log_msg_code, results->opf);
			}
			else 
			{
				sprintf(log_msg," -----: %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f", A[0],A[1],A[2],A[3],A[4],A[5]);
				lf_logmsg("DEBUG", "APP", 0, 0, log_msg, log_msg_code, results->opf);
				sprintf(log_msg," -----: %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f", erro[0],erro[1],erro[2],erro[3],erro[4],erro[5]);
				lf_logmsg("DEBUG", "APP", 0, 0, log_msg, log_msg_code, results->opf);
			}
		}
	
		//	Full ldfs                     
		int fldfr=0;	//proc result: 0=ok; 1=failed; 2=cannot be processed; 3=malloc pb; 4=not processed; 
		int fulldf_nrows, fulldf_ncols=2;
		int bins1=0, bins2=0;
		int retcode=0;
		float *save_full_ldf;
		if (fldf==1)
		{
			// full LDF

			//if (ret_gsl<0 || cmx < 0.01 || cmy < 0.01 || radius < 0.01) // that I don't remember what for	
			if (cmx < 0.01 || cmy < 0.01 || radius < 0.01)
			{
				fulldf_nrows=1;
				bins1=0;
				save_full_ldf  = (float *) malloc((2)*sizeof(float));
					if(!save_full_ldf) 
					{
						lf_logmsg("ERROR", "APP", ERR_MALLOC_FAILED, 0,"malloc failed (save_full_ldf)", log_msg_code, results->opf);
						return ERR_MALLOC_FAILED;
					}				
				save_full_ldf[0]=0;
				save_full_ldf[1]=0;
				fldfr=3;
			}
			else
			{
				retcode=mk_fldfs(cmx, cmy, radius, naxis_row, naxis_col, npixels, data, &save_full_ldf, &bins1, &bins2, results->opf, results->debug);
				if (retcode == ERR_MALLOC_FAILED)
				{
					fldfr=2;
					return ERR_MALLOC_FAILED;
				}
				else
					if (retcode == ERR_NR_STACK_TOO_SMALL)
					{
						ret_code=ERR_LIMBFIT_FLDF_FAILED;
						fldfr=1;
						//shouldnt exist?
					}
				fulldf_nrows=bins2;
			}
		}			


		/* ---------------------Save LDFs & AB data in a FITS file--------------
		
		 LDF = lprf(ldfs) & rprf (axis) come from limbfit.f			
		
					 to be read in IDL as: FLOAT     = Array[65, 182]

		  			NAXIS1          FLOAT               65
   					NAXIS2          FLOAT              182
		
			-----------------
			| axis	 |  fsn |
			-----------------
			| lfds	 |  ip  |
			| ...	 | 	... |
			| ....	 |  ... |
			-----------------
			| avg ldf|  ip  |
			-----------------
		
		--------------------------------------------------------------------- */
		// for saving them in FITS file
		float *save_ldf, *save_alpha_beta; 
		long ldf_nrow=nang+2, ldf_ncol=nprf+1;  //ii -nbr of ldfs, jj -nbr of points for each ldf
	    save_ldf 	= (float *) malloc((ldf_nrow*ldf_ncol)*sizeof(float));
			if(!save_ldf) 
			{
				lf_logmsg("ERROR", "APP", ERR_MALLOC_FAILED, 0,"malloc failed (save_ldf)", log_msg_code, results->opf);
				return ERR_MALLOC_FAILED;
			}
	    save_alpha_beta = (float *) malloc((ab_nrow*ab_ncol)*sizeof(float));
			if(!save_alpha_beta) 
			{
				lf_logmsg("ERROR", "APP", ERR_MALLOC_FAILED, 0,"malloc failed (save_alpha_beta)", log_msg_code, results->opf);
				return ERR_MALLOC_FAILED;
			}		

		float *p_sldf=&save_ldf[0];
		float *pl_sldf=&save_ldf[(ldf_nrow*ldf_ncol)-1];
		float *plc_sldf;

		// axis
		p_rprf=&rprf[0];
		pl_rprf=&rprf[nprf-1];
		while (p_rprf <= pl_rprf)
			*(p_sldf++)=*(p_rprf++);

//		*(p_sldf++)=(float)input->fsn;
		*(p_sldf++)=0.0;

/*old		
			
		//ldfs
		p_lprf=&lprf[0];
		p_sldf=&save_ldf[ldf_ncol]; //		shouldn't be? pl_lprf=&lprf[ldf_ncol]; ???
		while (p_sldf <= pl_sldf)
			*(p_sldf++)=*(p_lprf++);
*/			
		int cpt_row=0;
		p_lprf=&lprf[0];
		pl_lprf=&lprf[((nang+1)*nprf)-1];

		while (p_lprf <= pl_lprf)
		{
			plc_sldf=&lprf[((cpt_row+1)*nprf)-1];
			while (p_lprf <= plc_sldf)
				*(p_sldf++)=*(p_lprf++);
			*(p_sldf++)=*(p_t_ip++);
			cpt_row++;
		}

		/* ---------------------------------------------------------------------
		
		 AB 		write first alpha then beta
					to be read in IDL as:  Array[2, 256]
		
		  			NAXIS1          FLOAT                2
   					NAXIS2          FLOAT              256
 
		--------------------------------------------------------------------- */
		p_alph=pf_alph;
		p_b0=pf_b0;
		float *p_save_alpha_beta=&save_alpha_beta[0];
		while (p_alph <= pl_alph) 
		{
			*(p_save_alpha_beta++)=*(p_alph++);
			*(p_save_alpha_beta++)=*(p_b0++);
		}
		/*
		if (iter>1) 
		{
			p_sb0=pf_sb0;
			while (p_sb0 <= pl_sb0) *(p_save_alpha_beta++)=*(p_sb0++);
		}
		*/ 

		// Update Returned Structure when process succeeded                    
		results->cenx=cmx;
		results->ceny=cmy;
		results->radius=radius;
		results->cmean=cmean;
		results->numext=3;		
		results->fits_ldfs_naxis1=ldf_ncol;	
		results->fits_ldfs_naxis2=ldf_nrow;
		results->fits_fldfs_tfields=fulldf_ncols;
		results->fits_fldfs_nrows=fulldf_nrows;
		results->fits_ab_naxis1=ab_ncol;
		results->fits_ab_naxis2=ab_nrow;
		results->nb_fbins=bins1;
		results->ann_wd=w;
		results->mxszannv=S;
		results->nb_ldf=nang;
		results->nb_rdb=nprf;
		results->nb_abb=nreg;
		results->up_limit=rsi;
		results->lo_limit=rso;
		results->inc_x=dx;
		results->inc_y=dy;
		results->nfitpnts=r_size;
		results->nb_iter=iter;
		results->fldfr=fldfr;
		results->ahi=lahi;
		results->error1=ifail;		
		results->error2=ret_gsl;		
		if (ifail == 0) 
			results->quality=9; 
		//? if (cont>=nang && ret_gsl<0)
		if (ret_gsl<0)
		{
			results->quality=1; 
			ret_code=ERR_LIMBFIT_FIT_FAILED;
			if (results->debug)
			{	
				sprintf(log_msg,"  ret_gsl<0 = %2d", ret_gsl);
				lf_logmsg("DEBUG", "APP", 0, 0, log_msg, log_msg_code, results->opf);
	 		}     	
		}
		results->fits_ldfs_data=save_ldf; 
		results->fits_alpha_beta=save_alpha_beta; 
		if (fldf == 1) results->fits_fulldfs=save_full_ldf; else results->fits_fulldfs=0;
		free(D);
		free(LDF);
		free(t_ip);
	} // end limb OK
	else 
	{
		if (results->debug)
		{	
			sprintf(log_msg," limb.f routine returned ifail = %2d", ifail);
			lf_logmsg("DEBUG", "APP", 0, 0, log_msg, log_msg_code, results->opf);
 		}     	
		// Update Returned Structure when process failed                    
		results->numext=0;		
		results->error1=ifail;
		results->error2=-1;
		results->quality=0;
      	ret_code=ERR_LIMBFIT_FAILED;
		results->ann_wd=w;
		results->mxszannv=S;
		results->nb_ldf=nang;
		results->nb_rdb=nprf;
		results->nb_abb=nreg;
		results->up_limit=rsi;
		results->lo_limit=rso;
		results->inc_x=dx;
		results->inc_y=dy;
		results->nfitpnts=r_size;
		results->nb_iter=iter;
		results->ahi=lahi;
		results->cmean=0;
		results->nb_fbins=0;
		results->fldfr=4;
	} // end limb failed
	// IS: do not free those (save_ldf,save_params,save_alpha_beta,save_full_ldf) passed from or to the structure !

	free(rprf);
	free(alph);
	free(beta);
	free(b0);
	//free(sb0);
//	free(lprf);

	if (results->debug)
	{	
		sprintf(log_msg," >>>>end of limbfit with: %d", ret_code);
		lf_logmsg("DEBUG", "APP", 0, 0, log_msg, log_msg_code, results->opf);
	}
	sprintf(log_msg," end: RC: %d - limb:%d - fit:%d - quality: %d", ret_code, results->error1, results->error2, results->quality);
	lf_logmsg("INFO", "APP", 0, 0, log_msg, log_msg_code, results->opf);

return(ret_code);

} /* end of main */


/*--------------------------------------------------------------------------*/
int gaussfit(double y[], double t[],double sigma[], double A[], double erro[], long N, int degf, int debug, FILE *opf)
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
	static char *log_msg_code="gaussfit";
	int ret_code=0;
	char log_msg[200];

	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;
	int status;
	unsigned int i, iter = 0;
	const size_t n = N;
	const size_t p = degf;
		 
	gsl_matrix *covar = gsl_matrix_alloc (p, p);
	/* double t[N], y[N], sigma[N]; */
	struct dataF d = { n, t, y, sigma};
	gsl_multifit_function_fdf f;
	
	/* Initial Guess */
	double x_init[degf];
	for (i=0; i <degf; i++) 
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

	T = gsl_multifit_fdfsolver_lmsder;
	s = gsl_multifit_fdfsolver_alloc (T, n, p);
	gsl_multifit_fdfsolver_set (s, &f, &x.vector);  

	do
	{
		iter++;
		status = gsl_multifit_fdfsolver_iterate (s);     

		if (status==0 || status == GSL_ENOPROG)	
			status = gsl_multifit_test_delta (s->dx, s->x,1e-4, 1e-4);  		/* IS: change epsrel (last arg of gsl_multifit_test_delta) after the SDO Launch */ 
		else 
		{
			ret_code=ERR_GSL_GAUSSFIT_FDFSOLVER_FAILED;	
			if (debug) 
			{
				sprintf(log_msg," iter: %u, status: %d (%s)\n", iter,status,gsl_strerror (status));
				lf_logmsg("WARNING", "APP", ret_code, status, log_msg, log_msg_code, opf);			
			}
		}

	}
	while (status == GSL_CONTINUE && iter < 500);

			
	if (status == 0)  
	{

		gsl_multifit_covar (s->J, 0.0, covar);
			 
		#define FIT(i) gsl_vector_get(s->x, i)
		#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
			 
		double chi = gsl_blas_dnrm2(s->f);
		double dof = (double)(n - p);
		double c = GSL_MAX_DBL(1, chi / sqrt(dof)); 
		  
		
		for (i=0; i <degf; i++) {
			A[i]=FIT(i);
			erro[i]=c*ERR(i);
		}  
	}
	else
	{
		ret_code=ERR_GSL_GAUSSFIT_FDFSOLVER_FAILED;			
		if (debug)
		{
			sprintf(log_msg," (gsl_multifit_fdfsolver_iterate) failed, (iter=%u) gsl_errno=%d", iter,status);
			lf_logmsg("ERROR", "APP", ret_code, status, log_msg, log_msg_code, opf);
		}
	}

	
	if (debug==2)
	{																		
		sprintf(log_msg," Nonlinear LSQ fitting status = %s (%d)", gsl_strerror (status),status);
		lf_logmsg("DEBUG", "APP", 0, 0, log_msg, log_msg_code, opf);
	}
	if (debug) 
	{
		sprintf(log_msg," finished at iter# %u\n", iter);
		lf_logmsg("DEBUG", "APP", 0, 0, log_msg, log_msg_code, opf);			
	}
	gsl_multifit_fdfsolver_free (s);
	gsl_matrix_free (covar);
	gsl_rng_free (r);
	return(ret_code);
}


/*--------------------------------------------------------------------------*/

double fin_min(double A[], double m, int range, int degf, int debug, FILE *opf)
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
	static char *log_msg_code="fin_min";
	int status;
	gsl_set_error_handler_off();
	int ret_code=0;
	char log_msg[200];
	
	 int iter = 0, max_iter = 1000;
	 const gsl_min_fminimizer_type *T;
	 gsl_min_fminimizer *s;
	 
	 gsl_function F;
	
	//double m_expected = m+0.01; //1e-8;
	//double a = m-5000, b = m+5000;
	double a = m-range, b = m+range;      // initially = 2
//	if (degf == 4)
	//struct exp_plus_quadratic_params params= { A[0], A[1], A[2], A[3] }; 
//	else
	struct exp_plus_quadratic_params params= { A[0], A[1], A[2], A[3], A[4], A[5] }; 
	  F.function = &exp_plus_quadratic_function;
	  F.params = &params;
		 
	  T = gsl_min_fminimizer_brent;
	  s = gsl_min_fminimizer_alloc (T);
	
	  status=gsl_min_fminimizer_set (s, &F, m, a, b);     
	  if (status)
	  {
			if (status == GSL_EINVAL) 
			{
				ret_code=ERR_GSL_FINMIN_NOMIN_FAILED;
				if (debug)
				{
					sprintf(log_msg," (gsl_min_fminimizer_set) doesn't find a minimum, m=%f", m);
					lf_logmsg("DEBUG", "APP", 0, status, log_msg, log_msg_code, opf);
				}
			} 
			else 
			{
				ret_code=ERR_GSL_FINMIN_SET_FAILED;
				if (debug)
				{
					sprintf(log_msg," (gsl_min_fminimizer_set) failed, gsl_errno=%d", status);
					lf_logmsg("ERROR", "APP", ret_code, status, log_msg, log_msg_code, opf);			
				}
			}
			gsl_min_fminimizer_free(s);
			return ret_code;
	  }
	  
	  do
	  {
	   iter++;
	   status = gsl_min_fminimizer_iterate (s);
		 
	   m = gsl_min_fminimizer_x_minimum (s);
	   a = gsl_min_fminimizer_x_lower (s);
	   b = gsl_min_fminimizer_x_upper (s);
		 
	   status = gsl_min_test_interval (a, b,1E-3, 0.0);
	  }
	  while (status == GSL_CONTINUE && iter < max_iter);
	  
	  if (status)		// regarder lequel plantait et verifier si on a autre chose que GSL_EINVAL comme erreur
		{				// ajouter debug pour les messages apres
			if (status == GSL_EINVAL) 
			{
				ret_code=ERR_GSL_FINMIN_PRO_FAILED;
				if (debug)
				{
					sprintf(log_msg," invalid argument, n=%8.5f", a);
					lf_logmsg("ERROR", "APP", ret_code, status, log_msg, log_msg_code, opf);			
				}
			} 
			else 
			{			
				ret_code=ERR_GSL_FINMIN_PRO_FAILED;
				if (debug)
				{
					sprintf(log_msg," failed, gsl_errno=%d", status);
					lf_logmsg("ERROR", "APP", ret_code, status, log_msg, log_msg_code, opf);			
				}
			}
			gsl_min_fminimizer_free(s);
			return ret_code; 
		}
	  
	  if (debug==2)
	  {
			sprintf(log_msg," a:%8.5f b:%8.5f m:%8.5f iter=%d, b-a:%f", a,b,m,iter,b-a);
			lf_logmsg("DEBUG", "APP", 0, 0, log_msg, log_msg_code, opf);
	  }


	  gsl_min_fminimizer_free (s);
		 
	  return m;

}
//---------------------------------------------------------------------------------------------------------------------
// Make full ldfs
//---------------------------------------------------------------------------------------------------------------------
float median(float * tmed, int siz)
{
	int m=siz%2;
	int s=siz/2;
	return m==0?(tmed[s-1]+tmed[s])/2:(tmed[s]);
}

int mk_fldfs(float cmx, float cmy, double radius, int naxis_row, int naxis_col, long npixels, 
					float *data, float **save_full_ldf, int *bins1, int *bins2, FILE *opf, int debug)
{
static char *log_msg_code="mk_fldfs";

	int status=0;
	int retcode=0;
	int fulldf_nrows, fulldf_ncols=2;
	// compute  radius array
	float *data2 = (float *) malloc(sizeof(float) * npixels);
		if(!data2) 
		{
			lf_logmsg("ERROR", "APP", ERR_MALLOC_FAILED, 0,"malloc failed (data2)", log_msg_code, opf);
			return ERR_MALLOC_FAILED;
		}

	unsigned long ti,tx,ty;
	float tdx,tdy;
	unsigned long cnt=0;
	for (tx=0;tx<naxis_col;tx++)
	{
		for (ty=0;ty<naxis_row;ty++)
		{
			ti=naxis_col*ty+tx;
			if (data[ti]>-2147483648.)
			{	
				tdx=(float)fabs(tx-cmx);
				tdy=(float)fabs(ty-cmy);
				data2[ti]=(float)sqrt(tdx*tdx + tdy*tdy);
				cnt++;
			} else data2[ti]=900000.;	
		}
	}
		if (debug) 
			{
				lf_logmsg("DEBUG", "APP", DEBUG_MSG, status, "building array", log_msg_code, opf);			
			}
	// index them
	unsigned long *indx;
	indx=lvector(1,npixels,&status); 
		if (status<0) 
		{
			lf_logmsg("ERROR", "APP", ERR_MALLOC_FAILED, 0,"malloc failed", "lvector(indx)", opf);
			return ERR_MALLOC_FAILED;
		}
		if (debug) 
			{
				lf_logmsg("DEBUG", "APP", DEBUG_MSG, status, "vector", log_msg_code, opf);			
			}	retcode=indexx(npixels,data2,indx);	
		if (retcode<0) 
		{
			lf_logmsg("ERROR", "APP", ERR_NR_STACK_TOO_SMALL, 0,"stack too small", "indexx", opf);
			return ERR_NR_STACK_TOO_SMALL;
		}
		if (debug) 
			{
				lf_logmsg("DEBUG", "APP", DEBUG_MSG, status, "indexx", log_msg_code, opf);			
			}
	// make bins        	
	int rc=4; // size in pixels of DeltaR of the last bin before the limb: 2*the distorsion
	float v1=(float)(1.-(rc/radius));
	float v2=1/(1-(v1*v1));
	int bins=(int)v2-1;		
	unsigned int st=(unsigned int)(M_PI*radius*radius);
	int sb=st/bins;
	int ns=round(sb);
	int sr=st%bins;
	if (sr==0) *bins1=bins; else *bins1=bins+1;
	int mbins=bins+5; // to get what is over the limb: +1 = residual, then +4 outside
	fulldf_nrows=mbins;
	*bins2=mbins;
	// compute them	
	float *ttmed1, *ttmed2;
	// printf("pixels %ld cnt: %lu BINS: %d, MBINS: %d rs: %f st=%u, sb=%d, sr=%d, ns=%d \n",npixels,cnt,bins,mbins,radius,st,sb,sr,ns);

	ttmed1=vector(1,ns,&status); 
		if (status<0) 
		{
			lf_logmsg("ERROR", "APP", ERR_MALLOC_FAILED, 0,"malloc failed", "vector(ttmed1)", opf);
			return ERR_MALLOC_FAILED;
		}
	ttmed2=vector(1,ns,&status); 
		if (status<0) 
		{
			lf_logmsg("ERROR", "APP", ERR_MALLOC_FAILED, 0,"malloc failed", "vector(ttmed2)", opf);
			return ERR_MALLOC_FAILED;
		}
	float *t_med_int 	= (float *) malloc((mbins)*sizeof(float));
		if(!t_med_int) 
		{
			lf_logmsg("ERROR", "APP", ERR_MALLOC_FAILED, 0,"malloc failed (t_med_int)", log_msg_code, opf);
			return ERR_MALLOC_FAILED;
		}
	float *t_med  		= (float *) malloc((mbins)*sizeof(float));
		if(!t_med) 
		{
			lf_logmsg("ERROR", "APP", ERR_MALLOC_FAILED, 0,"malloc failed (t_med)", log_msg_code, opf);
			return ERR_MALLOC_FAILED;
		}
	int tk,next_i=0;
	for(tk=0;tk<mbins;tk++)
	{
		if (tk<bins || tk>bins)
		{
			for(ti=0;ti<ns;ti++) 
			{        		
				ttmed1[ti]=data[indx[next_i+ti]];
				ttmed2[ti]=data2[indx[next_i+ti]];
			}
			retcode=sort(ns,ttmed1);
				if (retcode<0) 
				{
					lf_logmsg("ERROR", "APP", ERR_NR_STACK_TOO_SMALL, 0,"stack too small", "sort(1)", opf);
					return ERR_NR_STACK_TOO_SMALL;
				}
/*				if (debug) 
				{
					lf_logmsg("DEBUG", "APP", NULL, status, "sort 1a", log_msg_code, opf);			
				}			retcode=sort(ns,ttmed2);
*/
			retcode=sort(ns,ttmed2);
				if (retcode<0) 
				{
					lf_logmsg("ERROR", "APP", ERR_NR_STACK_TOO_SMALL, 0,"stack too small", "sort(2)", opf);
					return ERR_NR_STACK_TOO_SMALL;
				}
/*				if (debug) 
				{
					lf_logmsg("DEBUG", "APP", NULL, status, "sort 2a", log_msg_code, opf);			
				}
*/
			t_med_int[tk]=median(ttmed1,ns);
			t_med[tk]=median(ttmed2,ns);
			next_i=(tk+1)*ns;
		} 
		else if (tk == bins && sr != 0)
		{
			for(ti=0;ti<sr;ti++) 
			{        		
				ttmed1[ti]=data[indx[next_i+ti]];
				ttmed2[ti]=data2[indx[next_i+ti]];
			}
			retcode=sort(sr,ttmed1);
				if (retcode<0) 
				{
					lf_logmsg("ERROR", "APP", ERR_NR_STACK_TOO_SMALL, 0,"stack too small", "sort(3)", opf);
					return ERR_NR_STACK_TOO_SMALL;
				}
/*				if (results->debug) 
				{
					lf_logmsg("DEBUG", "APP", NULL, status, "sort 1b", log_msg_code, results->opf);			
				}
*/
			retcode=sort(sr,ttmed2);
				if (retcode<0) 
				{
					lf_logmsg("ERROR", "APP", ERR_NR_STACK_TOO_SMALL, 0,"stack too small", "sort(4)", opf);
					return ERR_NR_STACK_TOO_SMALL;
				}
/*				if (results->debug) 
				{
					lf_logmsg("DEBUG", "APP", NULL, status, "sort 2b", log_msg_code, results->opf);			
				}
*/
			t_med_int[tk]=median(ttmed1,sr);
			t_med[tk]=median(ttmed2,sr);
			next_i=next_i+sr;
		}
	}
	if (debug) 
		{
			lf_logmsg("DEBUG", "APP", DEBUG_MSG, status, "after all sorts", log_msg_code, opf);			
		}
	// combine the 2 arrays = > pas necessaire pour l'integration dans le vrai fits	
	*save_full_ldf  = (float *) malloc((fulldf_nrows*fulldf_ncols)*sizeof(float));
		if(!save_full_ldf) 
		{
			lf_logmsg("ERROR", "APP", ERR_MALLOC_FAILED, 0,"malloc failed (save_full_ldf)", log_msg_code, opf);
			return ERR_MALLOC_FAILED;
		}


	float *p_full1=&t_med_int[0];
	float *pl_full1=&t_med_int[fulldf_nrows-1];
	float *p_full2=&t_med[0];
	float *pl_full2=&t_med[fulldf_nrows-1]; 
	float *p_save_full_ldf=save_full_ldf[0]; //&save_full_ldf2[0];
	while (p_full1 <= pl_full1) *(p_save_full_ldf++)=*(p_full1++);
	while (p_full2 <= pl_full2) *(p_save_full_ldf++)=*(p_full2++);	
	
	// free
	free(data2);
	free(t_med_int);
	free(t_med);
	// plus all others...  vectors...
	free_lvector(indx,1,npixels);
	free_vector(ttmed1,1,ns);
	free_vector(ttmed2,1,ns);
				if (debug) 
				{
					lf_logmsg("DEBUG", "APP", DEBUG_MSG, status, "end flds", log_msg_code, opf);			
				}
return 0;
}
