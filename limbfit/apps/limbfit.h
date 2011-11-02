/* I.Scholl "Sat Oct 22 09:44:08 HST 2011" 
*/

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

#include "fitsio.h"

#include <jsoc_main.h>

#include "expmax.h"
#include "expfit.h"

#define CODE_NAME 		"limbfit"
#define CODE_VERSION 	"V1.14r0" 
#define CODE_DATE 		"Sat Oct 22 09:44:08 HST 2011" 
#define LOGMSG1			"LIMBFIT"
#define	JSD_NAME		"su_scholl.hmi_lf.jsd"

//#define dsin	"hmi.lev1c_nrt[]"
//#define dsout	"su_scholl.limbfit"
//#define LOG_DIR	"~/LOGS/"
//#define TMP_DIR	"~/TMP/"

#define NUMRECLEV1 12

/*
drms_open_records
drms_create_records
drms_array_create
drms_segment_write
drms_segment_write_from_file
drms_set_key_string for the final status of the current processed record (because if I can't write the final status of the processing even the record will be in a incoherent state...)
*/
//---------------------------------------------------ERRORS
//GENERAL FAILURES -> ABORT
#define ERR_USAGE 							-1
#define ERR_MALLOC_FAILED 					-11
#define ERR_SPECIAL							-100
#define ERR_DRMS_WRITE 						-200
#define ERR_DRMS_READ 						-201
#define WAR_DRMS_NORECORD 					201

//GENERAL FAILURES -> write errors
#define ERR_DRMS_WRITE_KW 					-300
#define ERR_DRMS_READ_MISSING_DATA 			-301
#define ERR_DRMS_READ_MISSING_KW 			-302
#define ERR_DRMS_READ_MISSING_XYR_LF		-303
#define ERR_FITSIO		 					-350
#define ERR_FITSIO2		 					-351 // for mini FITS on error
#define ERR_NR_STACK_TOO_SMALL 				-352

//LIMBFIT FAILED -> write errors
#define ERR_LIMBFIT_FAILED 					-501
#define ERR_LIMBFIT_FIT_FAILED 				-502 // error on exit from Marcelo's fitting routines
#define ERR_LIMBFIT_FLDF_FAILED				-503
#define ERR_DISK_OUTSIDE_IMAGE 				-511
#define ERR_SIZE_ANN_TOO_BIG 				-512
#define ERR_GSL_FINMIN_SET_FAILED 			-541
#define ERR_GSL_FINMIN_NOMIN_FAILED			-542
#define ERR_GSL_FINMIN_PRO_FAILED 			-543
#define ERR_GSL_GAUSSFIT_SET_FAILED 		-551
#define ERR_GSL_GAUSSFIT_FDFSOLVER_FAILED 	-552

//----------------------------- Processing status codes per record
#define PROCSTAT_OK						"OK"
#define PROCSTAT_NOK					"NOK"
#define PROCSTAT_NO_LF_FAILED			"LF_FAILED"
#define PROCSTAT_NO_LF_MISSVALS 		"NO_LF_MISSVALS"
#define PROCSTAT_NO_LF_DARKIMG 			"NO_LF_DARKIMG"
#define PROCSTAT_NO_LF_OPENLOOP 		"NO_LF_OPENLOOP"
#define PROCSTAT_NO_LF_DB_READ_PB 		"NO_LF_DB_READ_PB"
#define PROCSTAT_NO_LF_XYR_LF_MISSING 	"NO_LF_XYR_LF_MISSING"
#define PROCSTAT_NO_LF_DB_WRITE_PB 		"NO_LF_DB_WRITE_PB"
#define PROCSTAT_NO_LF_FITS_WRITE_PB 	"NO_LF_FITS_WRITE_PB"

//---------------------------------- LIMBFIT PARAMETERS
#define ANNULUS_WIDTH 200					// 
#define MAX_SIZE_ANN_VARS 4000000			// ! must be the same value than JPT in fortran code !
#define NUM_LDF 180							// 
#define NUM_RADIAL_BINS 64					// 
#define NUM_AB_BINS 180						// 
#define LO_LIMIT 32.0						// ! the sum of these 2 must be equal to ANNULUS_WIDTH 
#define UP_LIMIT 32.0						// 
#define INC_X -4.0							// 
#define INC_Y -4.0							// 
#define NUM_FITPNTS 16						// 2*NUM_FITPNTS<NUM_RADIAL_BINS
#define GUESS_RANGE 8						//
#define NB_ITER 2							//
#define BAD_PIXEL_VALUE -2147483648.0
//#define	SKIPGC 1							// skip the guess estimation, use X0/YO_LF
//#define	IFAC 0								// skip the center calculation, use X0/YO_LF
#define	AHI 70000.0							// 

// alternate parameters for low LDF threshold - needed for roll analysis
#define	AHI2 	  30000.0
#define LO_LIMIT2 24.0	
#define UP_LIMIT2 24.0
#define NB_ITER2  1

//------------------------------------------------------

typedef struct {
	float			*data;			// image to analyze
	int				img_sz0;
	int				img_sz1;
	int				cc;
	int				spe;
	int				iter;
	double			ix;
	double			iy;
	double			ir;
} LIMBFIT_INPUT;

typedef struct {		// output files content
		
	// general keywords
	int			numext;
	float		cenx;
	float		ceny;
	double		radius;
	double		cmean;
	float		max_limb;
	int			quality;
	int			error1;
	int			error2;
		
	// result data
	float*		fits_ldfs_data; 		// main table / segment
	//float*		fits_fulldfs; 			// extension #2
	float*		fits_alpha_beta;  	 	// extension #0
	double*		fits_params;   			// extension #1

	// info to describe extension dimensions
	long		fits_ldfs_naxis1;		//	ldf_nrow
	long		fits_ldfs_naxis2;		//	ldf_ncol
//	long		fits_fldfs_nrows;		// 	fldf_nrow
//	long 		fits_fldfs_tfields;		//	fldf_ncol
	long 		fits_ab_nrows;			//	alpha_beta_nrow
	long		fits_ab_tfields;		//	alpha_beta_ncol
	long 		fits_params_nrows;		//	params_nrow
	long		fits_params_tfields;	//	params_ncol

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
	int			cc;
	double		ahi;
	int 		nb_fbins;
	
	// extra for error management
	char*		dsin;
	char*		comment;
	char*		code_date;
	char*		code_version;
	char*		code_name;
} LIMBFIT_OUTPUT;

void get_sdate(char *sdate);
void lf_logmsg(char * type1, char * type2, int return_code, int status, char *message, char *code_name, FILE *opf);
void lf_logmsg4fitsio(char *log_msg,char *log_msg_code,char *kw,unsigned int fsn,int status, FILE *opf);
void close_on_error(DRMS_Record_t *record_in,DRMS_Record_t *record_out, DRMS_Array_t *data_array); //, FILE *opf);
int get_set_kw(int typ, char *kw, char *kw_txt, unsigned int fsn, DRMS_Record_t *record_in,DRMS_Record_t *record_out, fitsfile *outfptr, FILE *opf, int debug, int tbf, LIMBFIT_OUTPUT *lfr, int *status);
float median(float * tmed, int siz);
int gaussfit(double y[], double t[],double sigma[], double A[6], double erro[6], long N, int debug, FILE *fd);
double fin_min(double A[6], double m, int range, int debug, FILE *fd);
int limbfit(LIMBFIT_INPUT *info,LIMBFIT_OUTPUT *results,FILE *opf,int debug);
void limb_(float *anls, long *jk, float *cmx, float *cmy, float *r, int *nitr, int *ncut, int *nang, 
			int *nprf, float* rprf, float* lprf, int *nreg, float *rsi, float *rso, float *dx, float *dy, 
			int *jreg, int *jang, int *jprf, float* alph, float* beta, int *ifail, float* b0, int *centyp, float *lahi); 
int process_n_records(char * open_dsname, char *dsout, char *tmp_dir, FILE *opf, int cc, int spe, int iter, char *dsin, char *comment, int debug, int *status);
int do_one_limbfit(unsigned int fsn, DRMS_Record_t *record_in,DRMS_Record_t *record_out, char *tmp_dir, FILE *opf, int cc, int spe, int iter, char *dsin, char *comment, int debug, int *status);
int	write_mini_output(char * errcode, DRMS_Record_t *record_in,DRMS_Record_t *record_out,FILE *opf, int tbf, LIMBFIT_OUTPUT *lfr, int debug);
