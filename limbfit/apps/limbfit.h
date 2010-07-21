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

#include "expmax.h"
#include "expfit.h"

#define CODE_NAME 		"limbfit"
#define CODE_VERSION 	"V0.1r8" 
#define CODE_DATE 		"Fri Apr  9 20:20:36 HST 2010" 

#define ERR_MALLOC_FAILED -11
#define ERR_LIMBFIT_FAILED -21
#define ERR_DISK_OUTSIDE_IMAGE -31
#define ERR_SIZE_ANN_TOO_BIG -32
// more ERR_LIMBFIT_ failed codes...
#define ERR_GSL_FIN_MIN_SET_FAILED -41
#define ERR_GSL_FIN_MIN_PRO_FAILED -42
#define ERR_GSL_GAUSSFIT_SET_FAILED -51
#define ERR_GSL_GAUSSFIT_FDFSOLVER_FAILED -52
// more ERR_GSL_FAILED failed codes...

#ifdef HMI

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
#define BAD_PIXEL_VALUE -2147483648.0
// more to add from the C and Fortran codes TBD

#elif  MDI

#define ANNULUS_WIDTH 50
#define MAX_SIZE_ANN_VARS 300000
#define NUM_LDF 180
#define NUM_RADIAL_BINS 16
#define NUM_AB_BINS 180
#define LO_LIMIT 6.0
#define UP_LIMIT 6.0
#define INC_X -1.0
#define INC_Y -1.0
#define NUM_FITPNTS 8
#define GUESS_RANGE 2
#define BAD_PIXEL_VALUE -32767.0
// more to add from the C and Fortran code TBD

#endif


typedef struct {
	// input paramters to the limbfit routine
#ifdef HMI
	float		*data;			// image to analyze
#elif  MDI
	short		*data;			// image to analyze
#endif
	int 		annw;			// annulus width
	long		mxsz_annv;		// maximum size if annulus variables
	int			nb_ldf;			// number of output LDFs
	int			nb_rdb;			// number of radial bins
	int			nb_abb;			// number of bins for alha and beta
	float		lol;			// lower limit to use in center-finder
	float		upl;			// upper limit to use in center-finder
	float 		incx;			// x increment to step by in center-finder
	float 		incy;			// y increment to step by in center-finder
	int			nb_fit_pnts;	// number of points used in the fitting function
	int			guess_rng;		// range for guess (fin_min)
	int			img_szx;
	int			img_szy;
	
	// info for tuning 
	double		ix;
	double		iy;
	double		ir;
	int			iter;
	FILE		*fd;
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
	int			quality;
	int			error;
		
	// result data
	float*		fits_ldfs_data; 	// main data
	float*		fits_fulldfs; 		// extension #6
	double*		fits_as;   			// extension #1
	double*		fits_es;   			// extension #2
	double*		fits_r;    			// extension #3
	float*		fits_ab_dataa;		// extension #4
	float*		fits_ab_datab;		// extension #4
	float*		fits_ab_datab1;		// extension #4
	float*		fits_ab_datab2;		// extension #4
	double*		fits_ip;    		// extension #5

	// info to describe extension dimensions
	long		fits_ext0_naxis;	// 	for fits_ldfs_data
	long		fits_ext0_naxis1;	//	ldf_nrow
	long		fits_ext0_naxis2;	//	ldf_ncol
	long 		fits_ext1_nrows;	//	aeris_nrow
	long		fits_ext1_tfields;	//	aes_ncol
	long 		fits_ext2_nrows;	//	aeris_nrow
	long		fits_ext2_tfields;	//	aes_ncol
	long 		fits_ext3_nrows;	//	aeris_nrow
	long		fits_ext3_tfields;	//	1
	long 		fits_ext4_nrows;	//	ab_nrow
	long		fits_ext4_tfields;	//	ab_ncol
	long 		fits_ext5_nrows;	//	aeris_nrow
	long		fits_ext5_tfields;	//	1
	long		fits_ext6_nrows;	//	fulldf_nrow
	long		fits_ext6_tfields;	// 	for fits_fullldfs

	// info for tunning
	double		gxc;	
	double		gyc;
	double		gr;
	float		*annulus;
	long		annulus_jk;
	float		max_limb;

} LIMBFIT_OUTPUT;

int gaussfit(double y[], double t[],double sigma[], double A[6], double erro[6], long N, int debug, FILE *fd);
double fin_min(double A[6], double m, int range, int debug, FILE *fd);
int limbfit(LIMBFIT_INPUT *info,LIMBFIT_OUTPUT *results,int skpic,int debug);
void limb_(float *anls, long *jk, float *cmx, float *cmy, float *r, int *nitr, int *ncut, int *nang, 
			int *nprf, float* rprf, float* lprf, int *nreg, float *rsi, float *rso, float *dx, float *dy, 
			int *jreg, int *jang, int *jprf, float* alph, float* beta, int *ifail, float* b0, int *centyp); 
