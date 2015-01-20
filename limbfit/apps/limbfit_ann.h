/* I.Scholl 
	#define CODE_VERSION 	"V1.00" 
	#define CODE_DATE 		"Mon Sep 15 15:14:19 PDT 2014" 
*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>

#include "nrutil.h"

#include <jsoc_main.h>

#define CODE_NAME 		"limbfit_ann"
#define CODE_VERSION 	"V1r0" 
#define CODE_DATE 		"Mon Sep 15 15:14:19 PDT 2014" 
#define LOGMSG1			"LIMBFITS_ANN"
#define	JSD_NAME		"scholl_limbfit_ann.jsd"

//#define dsin	"hmi.lev1c_nrt[]"
//#define dsout	"su_scholl.limbfit_ann"
//#define LOG_DIR	"~/LOGS/"
//#define TMP_DIR	"~/TMP/"

#define NUMRECPERTRANS 960 // must be equal to Unitsize in the JSD

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
#define ERR_EXIT 							1
#define ERR_USAGE 							-2
#define ERR_MALLOC_FAILED 					-11
#define ERR_SPECIAL							-100
#define ERR_DRMS_WRITE 						-200
#define ERR_DRMS_READ 						-201
#define ERR_DRMS_ARRAY_CREATE 				-202
#define ERR_DRMS_SEGMENT_WRITE 				-203
#define ERR_DRMS_SEGMENT_LOOKUPNUM			-204
#define WAR_DRMS_NORECORD 					201
#define WAR_DRMS_FREE_ARRAY					202
#define DEBUG_MSG 							999
#define	VOID								0

//GENERAL FAILURES -> write errors
#define ERR_DRMS_WRITE_KW 					-300
#define ERR_DRMS_READ_MISSING_DATA 			-301
#define ERR_DRMS_READ_MISSING_KW 			-302
#define ERR_DRMS_READ_MISSING_XYR_LA		-303
#define ERR_NR_STACK_TOO_SMALL 				-352

//LIMBFIT FAILED -> write errors
#define ERR_LIMBANN_FAILED 					-501
#define ERR_DISK_OUTSIDE_IMAGE 				-511
#define ERR_SIZE_ANN_TOO_BIG 				-512

//----------------------------- Processing status codes per record
#define PROCSTAT_OK						"OK"
#define PROCSTAT_NOK					"NOK"
#define PROCSTAT_NO_LA_FAILED			"LF_FAILED"
#define PROCSTAT_NO_LA_MISSVALS 		"NO_LA_MISSVALS"
#define PROCSTAT_NO_LA_DARKIMG 			"NO_LA_DARKIMG"
#define PROCSTAT_NO_LA_OPENLOOP 		"NO_LA_OPENLOOP"
#define PROCSTAT_NO_LA_DB_READ_PB 		"NO_LA_DB_READ_PB"
#define PROCSTAT_NO_LA_XYR_LF_MISSING 	"NO_LA_XYR_LF_MISSING"
#define PROCSTAT_NO_LA_DB_WRITE_PB 		"NO_LA_DB_WRITE_PB"

//---------------------------------- LIMBFIT PARAMETERS
#define EQNANVAL -2147483648	// NAN equivalent
#define IMG_SIZE 16777216		// = 4096*4096
#define CENTX 2048.0
#define CENTY 2048.0
#define R_MAX 1985.0
#define R_MIN 1825.0
#define NAXIS_ROW 4096
#define NAXIS_COL 4096
//------------------------------------------------------

typedef struct { // input content
	int		*mask;	// mask image
	int		*pf_mask;		//
	int		*pl_mask;		//
	unsigned int fsn;
	double		ix;
	double		iy;
} LIMBFIT_INPUT;

typedef struct {	// output files content
		
	// general keywords
	float		cmean;
	int			quality;
	float		cenx;
	float		ceny;
	float		r_min;
	float		r_max;
			
	// extra for error management
	char*		dsin;
	char*		comment;
	char*		code_date;
	char*		code_version;
	char*		code_name;
	char*		bld_vers;
	char*		series_name;	
	
	// not to save
	FILE 		*opf;
	char*		tmp_dir;
	char*		dsout;
	int			debug;
	
} LIMBFIT_OUTPUT;

// C functions
void	close_on_error(DRMS_Record_t *record_in,DRMS_Record_t *record_out, DRMS_Array_t *data_array); //, FILE *opf);
int		do_one_limbfit(DRMS_Record_t *record_in,DRMS_Record_t *record_out, 
					LIMBFIT_INPUT *input, LIMBFIT_OUTPUT *results, int *status);
void	get_sdate(char *sdate);
void	lf_logmsg(char * type1, char * type2, int return_code, int status, char *message, char *code_name, FILE *opf);
int		process_n_records_fsn(char * open_dsname, LIMBFIT_INPUT *lfv, LIMBFIT_OUTPUT *lfr, int *status);
int		process_all_records_smpl(char * open_dsname, LIMBFIT_INPUT *lfv, LIMBFIT_OUTPUT *lfr, int *status);
int		write_mini_output(char * errcode, DRMS_Record_t *record_in,DRMS_Record_t *record_out,int tbf, LIMBFIT_OUTPUT *lfr);
int		write_lf_keywords(char * errcode, DRMS_Record_t *record_out, LIMBFIT_OUTPUT *results, int pass);

