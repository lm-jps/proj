/* I.Scholl Fri Apr  9 20:20:36 HST 2010 */

/* limbfit.h */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_rng.h>                                      
#include <gsl/gsl_randist.h>                                  
#include <gsl/gsl_vector.h>                                   
#include <gsl/gsl_blas.h>                                     
#include <gsl/gsl_multifit_nlin.h>                            
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_min.h>    
#include <time.h>

#include <jsoc_main.h>

#include "expmax.h"
#include "expfit.h"

#define CODE_NAME 		"limbfit"
#define CODE_VERSION 	"V0.2r2" 
#define CODE_DATE 		"Wed Jul 21 09:32:15 PDT 2010" 

//#define dsin	"hmi.lev1c_nrt[]"
//#define dsout	"su_scholl.limbfit"
#define LOG_DIR	"~/LOGS/"
#define TMP_DIR	"~/TMP/"

#define LOGMSG1	"LIMBFIT"

#define ERR_USAGE 							-1
#define ERR_MALLOC_FAILED 					-11
#define ERR_LIMBFIT_FAILED 					-21
#define ERR_DISK_OUTSIDE_IMAGE 				-31
#define ERR_SIZE_ANN_TOO_BIG 				-32
#define ERR_GSL_FIN_MIN_SET_FAILED 			-41
#define ERR_GSL_FIN_MIN_PRO_FAILED 			-42
#define ERR_GSL_GAUSSFIT_SET_FAILED 		-51
#define ERR_GSL_GAUSSFIT_FDFSOLVER_FAILED 	-52
#define ERR_SPECIAL							-100
#define ERR_DRMS_WR 						-200

#define ANNULUS_WIDTH 200					// -1
#define MAX_SIZE_ANN_VARS 4000000			// -2  ! must be the same value than JPT in fortran code !
#define NUM_LDF 180							// -3   //a verifier
#define NUM_RADIAL_BINS 64					// -4 	
#define NUM_AB_BINS 180						// -5
#define LO_LIMIT 24.0						// -6  ! the sum of these 2 must be equal to ANNULUS_WIDTH 
#define UP_LIMIT 24.0						// -7
#define INC_X -4.0							// -8
#define INC_Y -4.0							// -9
#define NUM_FITPNTS 16						// -10	//  2*NUM_FITPNTS<NUM_RADIAL_BINS
#define GUESS_RANGE 8						// -11
#define NB_ITER 4							// -12
#define BAD_PIXEL_VALUE -2147483648.0
#define	SKIPGC 1							// skip the guess estimation, use X0/YO_LF

typedef struct {
	float			*data;			// image to analyze
	int				img_sz0;
	int				img_sz1;
	FILE			*opf;
	char			*dsout;
	double			ix;
	double			iy;
	double			ir;
	int				camera;
	int				fid;
	unsigned int	fsn;
	int				hpltid;
	char			*tobs;
	float			sat_rot;
} LIMBFIT_INPUT;

typedef struct {		// output files content
		
	// general keywords
	char		proc_date[50];
	char		code_name[10];
	char		code_version[10];		
	int			numext;
	float		cenx;
	float		ceny;
	double		radius;
	double		cmean;
	float		max_limb;
	int			quality;
	int			error;
		
	// result data
	float*		fits_ldfs_data; 	// main data
	float*		fits_fulldfs; 		// extension #6
	float*		fits_alpha_beta;   	// extension #1
//	double*		fits_params;   		// extension #2
	double*		fits_as;   			// extension #1
	double*		fits_es;   			// extension #2
	double*		fits_r;    			// extension #3
	float*		fits_ab_dataa;		// extension #4
	float*		fits_ab_datab;		// extension #4
	double*		fits_ip;    		// extension #5

	// info to describe extension dimensions
	long		fits_ldfs_naxis1;		//	ldf_nrow
	long		fits_ldfs_naxis2;		//	ldf_ncol
//	long		fits_fldfs_naxis1;		// 	fldf_nrow
//	long 		fits_fldfs_naxis2;		//	fldf_ncol
//	long 		fits_alpha_beta_naxis1;	//	alpha_beta_nrow
//	long		fits_alpha_beta_naxis2;	//	alpha_beta_ncol
//	long 		fits_params_naxis1;		//	params_nrow
//	long		fits_params_naxis2;		//	params_ncol
	long 		fits_as_nrows;		//	aeris_nrow
	long		fits_as_tfields;	//	aes_ncol
	long 		fits_es_nrows;		//	aeris_nrow
	long		fits_es_tfields;	//	aes_ncol
	long 		fits_r_nrows;		//	aeris_nrow
	long		fits_r_tfields;		//	1
	long 		fits_ab_nrows;		//	ab_nrow
	long		fits_ab_tfields;	//	ab_ncol
	long 		fits_ip_nrows;		//	aeris_nrow
	long		fits_ip_tfields;	//	1
	long		fits_fldfs_nrows;	//	fulldf_nrow
	long		fits_fldfs_tfields;	// 	for fits_fullldfs


	// processing parameters to save
	int			ann_wd;
	long		mxszannv;
	int			nb_ldf;
	int			nb_rdb;
	int			nb_abb;
	double		up_limit;
	double		lo_limit;
	double		inc_x;
	double		inc_y;
	int			nfitpnts;
	int			nb_iter;
	int			skipgc;

} LIMBFIT_OUTPUT;

int gaussfit(double y[], double t[],double sigma[], double A[6], double erro[6], long N, int debug, FILE *fd);
double fin_min(double A[6], double m, int range, int debug, FILE *fd);
int limbfit(LIMBFIT_INPUT *info,LIMBFIT_OUTPUT *results,int debug);
void limb_(float *anls, long *jk, float *cmx, float *cmy, float *r, int *nitr, int *ncut, int *nang, 
			int *nprf, float* rprf, float* lprf, int *nreg, float *rsi, float *rso, float *dx, float *dy, 
			int *jreg, int *jang, int *jprf, float* alph, float* beta, int *ifail, float* b0, int *centyp); 
int do_one_limbfit(LIMBFIT_INPUT *info,DRMS_Record_t *record_in,int debug);
