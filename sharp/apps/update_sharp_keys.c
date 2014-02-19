/*                                                                      
 *      MODULE NAME: update_sharp_keys.c
 *
 *	DESCRIPTION: This module recalculates SHARP keywords. This is accomplished by 
 *      cloning a record, recalculating keywords of choice (i.e., user input), 
 *      and pointing to the same segments as the old record. Associated error keys
 *      are computed by default. 
 * 
 *      This module accounts for versioning by appending to the keyword CODEVER7 and HISTORY:
 *      CODEVER7 will contain multiple lines: the production build of sharp.c, the 
 *      production build of include file sw_functions.c, and the production build
 *      of update_sharp_keys.c. HISTORY will include a human-readable sentence about which
 *      keywords were updated. 
 * 
 *      This module does not produce segments.
 *
 *      INPUTS     : -- DRMS SHARP series
 *                   -- DRMS SHARP CEA series
 *                   -- HARPNUM
 *                   -- comma separated list of keywords to recalculate
 *
 *	AUTHOR     : Monica Bobra
 *
 *	Version    :   v0.0	Jun 14 2013
 *
 *	EXAMPLE    :
 *      update_sharp_keys sharpseries=hmi.sharp_720s sharpceaseries=hmi.sharp_cea_720s //
 *      HARPNUM=1 keylist=USFLUX,TOTPOT
 *
 */

/* Include files */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>
#include "jsoc_main.h"
#include "astro.h"
#include "fstats.h"
#include "cartography.c"
#include "fresize.h"
#include "finterpolate.h"
#include "img2helioVector.c"
#include "copy_me_keys.c"
#include "errorprop.c"
#include "sw_functions.c"
#include <mkl_blas.h>
#include <mkl_service.h>
#include <mkl_lapack.h>
#include <mkl_vml_functions.h>
#include <omp.h>
 
/* Define values */                     
#define PI              (M_PI)
#define RADSINDEG       (PI/180.)
#define ARRLENGTH(ARR) (sizeof(ARR) / sizeof(ARR[0]))
#define DIE(msg) {fflush(stdout); fprintf(stderr,"%s, status=%d\n", msg, status); return(status);}
#define SHOW(msg) {printf("%s", msg); fflush(stdout);}

char *module_name = "update_sharp_keys";  /* Module name */
char *version_id  = "2013 Jun 14";        /* Version number */

ModuleArgs_t module_args[] =
{
	{ARG_STRING, "sharpseries",     NULL, "Input/output Sharp dataseries"},
	{ARG_STRING, "sharpceaseries",  NULL, "Input/output Sharp CEA dataseries"},
	{ARG_INT,    "HARPNUM",         NULL, "HARP number"},
	{ARG_STRING, "keylist",         NULL, "comma separated list of keywords to update"},
	{ARG_END}
};


int DoIt(void)                    
{
        int errbufstat=setvbuf(stderr, NULL, _IONBF, BUFSIZ);
        int outbufstat=setvbuf(stdout, NULL, _IONBF, BUFSIZ);

	int status = DRMS_SUCCESS;
	int nrecs, irec;

        /* Keywords */
	float mean_vf; 
	float absFlux;
        float count_mask;
	float mean_hf;
	float mean_gamma;
	float mean_derivative_btotal;
	float mean_derivative_bh; 
	float mean_derivative_bz;
	float mean_jz;
	float us_i;
	float mean_alpha;
	float mean_ih;
	float total_us_ih;
	float total_abs_ih;
	float totaljz;
	float totpot;
	float meanpot;
	float area_w_shear_gt_45;
	float meanshear_angle;
	float area_w_shear_gt_45h;
	float meanshear_angleh; 
        float mean_derivative_btotal_err;
        float mean_vf_err;
        float mean_gamma_err;
        float mean_derivative_bh_err;
        float mean_derivative_bz_err;
        float mean_jz_err;
        float us_i_err; 
        float mean_alpha_err;
        float mean_ih_err;
        float total_us_ih_err;
        float total_abs_ih_err;
        float totaljz_err;
        float meanpot_err;
        float totpot_err;
        float meanshear_angle_err;
	char *sharpseries = (char *) params_get_str(&cmdparams, "sharpseries");
	char *sharpceaseries = (char *) params_get_str(&cmdparams, "sharpceaseries");

	int harpnum = params_get_int(&cmdparams, "HARPNUM");
	char sharpquery[256];
	char sharpceaquery[256];

	sprintf(sharpceaquery, "%s[%d]\n", sharpceaseries,harpnum);
	printf(sharpceaquery, "%s[%d]\n", sharpceaseries,harpnum);

	sprintf(sharpquery, "%s[%d]\n", sharpseries,harpnum);
	printf(sharpquery, "%s[%d]\n", sharpseries,harpnum);

	DRMS_RecordSet_t *sharpinrecset = drms_open_records(drms_env, sharpquery, &status);
	if (status != DRMS_SUCCESS)
		DIE("Problem opening sharp recordset.");
	if (sharpinrecset->n == 0)
		DIE("Sharp recordset contains no records.");
	DRMS_RecordSet_t *sharpoutrecset = drms_clone_records(sharpinrecset, DRMS_PERMANENT, DRMS_SHARE_SEGMENTS, &status);
	if (status != DRMS_SUCCESS)
		DIE("Problem cloning sharp records.");

	DRMS_RecordSet_t *sharpceainrecset = drms_open_records(drms_env, sharpceaquery, &status);
	if (status != DRMS_SUCCESS)
		DIE("Problem opening sharp cea recordset.");
	if (sharpceainrecset->n == 0)
		DIE("Sharp cea recordset contains no records.");
	DRMS_RecordSet_t *sharpceaoutrecset = drms_clone_records(sharpceainrecset, DRMS_PERMANENT, DRMS_SHARE_SEGMENTS, &status);
	if (status != DRMS_SUCCESS)
		DIE("Problem cloning sharp cea records.");

	char *keylist = (char *) params_get_str(&cmdparams, "keylist");

        // Flags to indicate which keyword will be recalculated
	int meanvfflag  = (strstr(keylist,"USFLUX")  != NULL);  // generalize so that lowercase is acceptable
	int totpotflag  = (strstr(keylist,"TOTPOT")  != NULL);  
        int meangamflag = (strstr(keylist,"MEANGAM") != NULL);
        int meangbtflag = (strstr(keylist,"MEANGBT") != NULL);
        int meangbzflag = (strstr(keylist,"MEANGBZ") != NULL);
        int meangbhflag = (strstr(keylist,"MEANGBH") != NULL);
        int meanjzdflag = (strstr(keylist,"MEANJZD") != NULL);
        int totusjzflag = (strstr(keylist,"TOTUSJZ") != NULL);
        int meanalpflag = (strstr(keylist,"MEANALP") != NULL);
        int meanjzhflag = (strstr(keylist,"MEANJZH") != NULL);
        int totusjhflag = (strstr(keylist,"TOTUSJH") != NULL);
        int absnjzhflag = (strstr(keylist,"ABSNJZH") != NULL);
        int savncppflag = (strstr(keylist,"SAVNCPP") != NULL);
        int meanpotflag = (strstr(keylist,"MEANPOT") != NULL);
        int meanshrflag = (strstr(keylist,"MEANSHR") != NULL);
        int shrgt45flag = (strstr(keylist,"SHRGT45") != NULL);
	DRMS_Record_t *sharpinrec = sharpinrecset->records[0];
	DRMS_Record_t *sharpceainrec = sharpceainrecset->records[0];
	DRMS_Segment_t *inseg = drms_segment_lookup(sharpceainrec, "Br");
	int nx = inseg->axis[0];
	int ny = inseg->axis[1];
	int nxny = nx * ny;
	int dims[2] = {nx, ny};

	// Temp arrays 	
	float *bh      = (float *) (malloc(nxny * sizeof(float)));
	float *bt      = (float *) (malloc(nxny * sizeof(float)));
	float *jz      = (float *) (malloc(nxny * sizeof(float)));
	float *bpx     = (float *) (malloc(nxny * sizeof(float)));
	float *bpy     = (float *) (malloc(nxny * sizeof(float)));
	float *bpz     = (float *) (malloc(nxny * sizeof(float)));
	float *derx    = (float *) (malloc(nxny * sizeof(float)));
	float *dery    = (float *) (malloc(nxny * sizeof(float)));
	float *derx_bt = (float *) (malloc(nxny * sizeof(float)));
	float *dery_bt = (float *) (malloc(nxny * sizeof(float)));
	float *derx_bh = (float *) (malloc(nxny * sizeof(float)));
	float *dery_bh = (float *) (malloc(nxny * sizeof(float)));
	float *derx_bz = (float *) (malloc(nxny * sizeof(float)));
	float *dery_bz = (float *) (malloc(nxny * sizeof(float)));
	float *bt_err  = (float *) (malloc(nxny * sizeof(float)));
	float *bh_err  = (float *) (malloc(nxny * sizeof(float)));
        float *jz_err  = (float *) (malloc(nxny * sizeof(float)));
        float *jz_err_squared = (float *) (malloc(nxny * sizeof(float)));
        float *jz_rms_err = (float *) (malloc(nxny * sizeof(float)));
        float *jz_err_squared_smooth = (float *) (malloc(nxny * sizeof(float)));
	float *jz_smooth = (float *) (malloc(nxny * sizeof(float)));

        // ephemeris variables
	float  cdelt1_orig, cdelt1, dsun_obs, imcrpix1, imcrpix2, crpix1, crpix2;
	double rsun_ref, rsun_obs;

        // outrecords 
	DRMS_Record_t *sharpoutrec, *sharpceaoutrec;
	nrecs=sharpoutrecset->n;
        printf("nrecs=%d\n",nrecs);

        // prepare to set CODEVER7 (CVS Version of the SHARP module)
	char *cvsinfo0;
	char *history0;
	char *cvsinfo1 = strdup("$Id: update_sharp_keys.c,v 1.5 2014/02/19 14:59:25 arta Exp $");
	char *cvsinfo2 = sw_functions_version();
	char *cvsinfoall = (char *)malloc(2048);
        char historyofthemodule[2048];
        cvsinfo0    = drms_getkey_string(sharpinrec, "CODEVER7", &status);
        strcpy(cvsinfoall,cvsinfo0);
        strcat(cvsinfoall,"\n");
        strcat(cvsinfoall,cvsinfo2);

        // prepare to set HISTORY (History of the data)
	char timebuf[1024];
	float UNIX_epoch = -220924792.000; 	
	sprint_time(timebuf, (double)time(NULL) + UNIX_epoch, "ISO", 0);
        sprintf(historyofthemodule,"The following keywords were re-computed on %s: %s.",timebuf,keylist); 
	printf("historyofthemodule=%s\n",historyofthemodule);
	history0    = drms_getkey_string(sharpinrec, "HISTORY", &status);
	strcat(historyofthemodule,"\n");
        strcat(historyofthemodule,history0);
 
        // no longer need inrecs
 	drms_close_records(sharpinrecset, DRMS_FREE_RECORD);
	drms_close_records(sharpceainrecset, DRMS_FREE_RECORD);
        
	for (irec=0;irec<nrecs;irec++)
	{
   	   // Get emphemeris
  	   sharpoutrec = sharpoutrecset->records[irec];
  	   sharpceaoutrec = sharpceaoutrecset->records[irec];
	   cdelt1_orig = drms_getkey_float(sharpceaoutrec, "CDELT1",   &status);
	   dsun_obs    = drms_getkey_float(sharpoutrec, "DSUN_OBS",   &status);
	   rsun_ref    = drms_getkey_double(sharpoutrec, "RSUN_REF", &status);
	   rsun_obs    = drms_getkey_double(sharpoutrec, "RSUN_OBS", &status);
	   imcrpix1    = drms_getkey_float(sharpoutrec, "IMCRPIX1", &status);
	   imcrpix2    = drms_getkey_float(sharpoutrec, "IMCRPIX2", &status);
	   crpix1      = drms_getkey_float(sharpoutrec, "CRPIX1", &status);
	   crpix2      = drms_getkey_float(sharpoutrec, "CRPIX2", &status);
           cdelt1      = (atan((rsun_ref*cdelt1_orig*RADSINDEG)/(dsun_obs)))*(1/RADSINDEG)*(3600.);

	   // set CODEVER7 and HISTORY
	   drms_setkey_string(sharpoutrec, "CODEVER7", cvsinfoall);
           drms_setkey_string(sharpoutrec,"HISTORY",historyofthemodule);           
	   drms_setkey_string(sharpceaoutrec, "CODEVER7", cvsinfoall);
           drms_setkey_string(sharpceaoutrec,"HISTORY",historyofthemodule);

	   // Get bx, by, bz, mask	
	   // Use HARP (Turmon) bitmap as a threshold on spaceweather quantities
	   DRMS_Segment_t *bitmaskSeg = drms_segment_lookup(sharpceaoutrec, "bitmap");
	   DRMS_Array_t *bitmaskArray = drms_segment_read(bitmaskSeg, DRMS_TYPE_INT, &status);
	   int *bitmask = (int *) bitmaskArray->data;	 // get the previously made mask array	

	   //Use conf_disambig map as a threshold on spaceweather quantities 
	   DRMS_Segment_t *maskSeg = drms_segment_lookup(sharpceaoutrec, "conf_disambig");       
	   DRMS_Array_t *maskArray = drms_segment_read(maskSeg, DRMS_TYPE_INT, &status);	
           if (status != DRMS_SUCCESS)
              DIE("No CONF_DISAMBIG segment.");
	   int *mask = (int *) maskArray->data;		// get the previously made mask array	

	   DRMS_Segment_t *bxSeg = drms_segment_lookup(sharpceaoutrec, "Bp");
	   DRMS_Array_t *bxArray = drms_segment_read(bxSeg, DRMS_TYPE_FLOAT, &status);
	   float *bx = (float *) bxArray->data;		// bx
	
	   DRMS_Segment_t *bySeg = drms_segment_lookup(sharpceaoutrec, "Bt");
	   DRMS_Array_t *byArray = drms_segment_read(bySeg, DRMS_TYPE_FLOAT, &status);
	   float *by = (float *) byArray->data;		// by
	   for (int i = 0; i < nxny; i++) by[i] *= -1;
	
	   DRMS_Segment_t *bzSeg = drms_segment_lookup(sharpceaoutrec, "Br");
	   DRMS_Array_t *bzArray = drms_segment_read(bzSeg, DRMS_TYPE_FLOAT, &status);
	   float *bz = (float *) bzArray->data;		// bz

   	   DRMS_Segment_t *bz_errSeg = drms_segment_lookup(sharpceaoutrec, "Br_err");
	   DRMS_Array_t *bz_errArray = drms_segment_read(bz_errSeg, DRMS_TYPE_FLOAT, &status);
	   float *bz_err = (float *) bz_errArray->data;		// bz_err

	   DRMS_Segment_t *by_errSeg = drms_segment_lookup(sharpceaoutrec, "Bt_err");
	   DRMS_Array_t *by_errArray = drms_segment_read(by_errSeg, DRMS_TYPE_FLOAT, &status);
	   float *by_err = (float *) by_errArray->data;		// by_err

	   DRMS_Segment_t *bx_errSeg = drms_segment_lookup(sharpceaoutrec, "Bp_err");
	   DRMS_Array_t *bx_errArray = drms_segment_read(bx_errSeg, DRMS_TYPE_FLOAT, &status);
	   float *bx_err = (float *) bx_errArray->data;		// bx_err

           /***** USFLUX, Example Function 1 ************************************sdfdsf*/ 
	   if (meanvfflag)
	   {
              // Compute unsigned flux 
              if (computeAbsFlux(bz_err, bz , dims, &absFlux, &mean_vf,  &mean_vf_err, 
                                 &count_mask, mask, bitmask, cdelt1, rsun_ref, rsun_obs))
              {
	         mean_vf           = DRMS_MISSING_FLOAT;
                 mean_vf_err       = DRMS_MISSING_FLOAT;
                 count_mask        = DRMS_MISSING_INT;
	      }

           drms_setkey_float(sharpoutrec, "USFLUX",  mean_vf);
           drms_setkey_float(sharpoutrec, "CMASK",   count_mask); 	
           drms_setkey_float(sharpoutrec, "ERRVF",   mean_vf_err);

           drms_setkey_float(sharpceaoutrec, "USFLUX",  mean_vf);
           drms_setkey_float(sharpceaoutrec, "CMASK",   count_mask); 	
           drms_setkey_float(sharpceaoutrec, "ERRVF",   mean_vf_err);
           }
   
           /***** MEANPOT and TOTPOT, Example Function 13  (Requires Keiji's code) **/
       
   	   if (totpotflag || meanpotflag)
	   {
              // First compute potential field	
              for (int i = 0; i < nxny; i++) bpz[i] = bz[i];
	      greenpot(bpx, bpy, bpz, nx, ny);
 
              // Then compute energy   
	      if (computeFreeEnergy(bx_err, by_err, bx, by, bpx, bpy, dims, 
                                    &meanpot, &meanpot_err, &totpot, &totpot_err, 
			            mask, bitmask, cdelt1, rsun_ref, rsun_obs)) 
              {
	         meanpot           = DRMS_MISSING_FLOAT; 
	         totpot            = DRMS_MISSING_FLOAT;
                 meanpot_err       = DRMS_MISSING_FLOAT;
                 totpot_err        = DRMS_MISSING_FLOAT;
	      }
	   
              // Set sharp keys          
              drms_setkey_float(sharpoutrec, "MEANPOT",  meanpot); 
              drms_setkey_float(sharpoutrec, "TOTPOT",   totpot); 	
              drms_setkey_float(sharpoutrec, "ERRMPOT",  meanpot_err); 
              drms_setkey_float(sharpoutrec, "ERRTPOT",  totpot_err);  

              // Set sharp cea keys
              drms_setkey_float(sharpceaoutrec, "MEANPOT",  meanpot); 
              drms_setkey_float(sharpceaoutrec, "TOTPOT",   totpot); 	
              drms_setkey_float(sharpceaoutrec, "ERRMPOT",  meanpot_err); 
              drms_setkey_float(sharpceaoutrec, "ERRTPOT",  totpot_err);
           }
       
           /***** MEANGAM, Example Function 3 ************************************/
   
           if (meangamflag)
           {
              // First compute horizontal field
              computeBh(bx_err, by_err, bh_err, bx, by, bz, bh, dims, &mean_hf, mask, bitmask);

              // Then compute inclination angle, gamma
              if (computeGamma(bz_err, bh_err, bx, by, bz, bh, dims, &mean_gamma, &mean_gamma_err,mask, bitmask))
	      {	
                 mean_gamma                 =  DRMS_MISSING_FLOAT;
                 mean_gamma_err             =  DRMS_MISSING_FLOAT;
              }

              // Set sharp keys    
              drms_setkey_float(sharpoutrec, "MEANGAM",  mean_gamma); 
              drms_setkey_float(sharpoutrec, "ERRGAM",   mean_gamma_err); 	

              // Set sharp cea keys
              drms_setkey_float(sharpceaoutrec, "MEANGAM", mean_gamma); 
              drms_setkey_float(sharpceaoutrec, "ERRGAM",  mean_gamma_err); 	
           }

           /***** MEANGBT, Example Function 5 (Requires Function 4) *************/

           if (meangbtflag)
           {
              // First, compute Bt 
              computeB_total(bx_err, by_err, bz_err, bt_err, bx, by, bz, bt, dims, mask, bitmask);

              // Then take the derivative of Bt
	      if (computeBtotalderivative(bt, dims, &mean_derivative_btotal, mask, bitmask, derx_bt, dery_bt, bt_err, 
                                          &mean_derivative_btotal_err))
              {
	         mean_derivative_btotal     = DRMS_MISSING_FLOAT;
	         mean_derivative_btotal_err = DRMS_MISSING_FLOAT;
              }

              // Set sharp keys    
              drms_setkey_float(sharpoutrec, "MEANGBT",  mean_derivative_btotal); 
              drms_setkey_float(sharpoutrec, "ERRBT",  mean_derivative_btotal_err); 
	
              // Set sharp cea keys  
              drms_setkey_float(sharpceaoutrec, "MEANGBT", mean_derivative_btotal); 
              drms_setkey_float(sharpceaoutrec, "ERRBT", mean_derivative_btotal_err); 	
           }
   
           /***** MEANGBH, Example Function 6 (Requires Function 2) *************/

           if (meangbhflag)
           {
              // First, compute Bh
   	      computeBh(bx_err, by_err, bh_err, bx, by, bz, bh, dims, &mean_hf, mask, bitmask);

              // Then take the derivative of Bh
	      if (computeBhderivative(bh, bh_err, dims, &mean_derivative_bh, &mean_derivative_bh_err, mask, bitmask, derx_bh, dery_bh))
              {
	         mean_derivative_bh       = DRMS_MISSING_FLOAT;
                 mean_derivative_bh_err   = DRMS_MISSING_FLOAT;
	      }

              // Set sharp keys   
              drms_setkey_float(sharpoutrec, "MEANGBH", mean_derivative_bh); 
              drms_setkey_float(sharpoutrec, "ERRBH", mean_derivative_bh_err); 
	
              // Set sharp cea keys  
              drms_setkey_float(sharpceaoutrec, "MEANGBH", mean_derivative_bh); 
              drms_setkey_float(sharpceaoutrec, "ERRBH", mean_derivative_bh_err); 	
           }

           /***** MEANGBZ, Example Function 7 ************************************/

           if (meangbzflag)
           {
              // Compute Bz derivative
              if (computeBzderivative(bz, bz_err, dims, &mean_derivative_bz, &mean_derivative_bz_err, mask, bitmask, derx_bz, dery_bz))
              {
	         mean_derivative_bz     = DRMS_MISSING_FLOAT; 
                 mean_derivative_bz_err = DRMS_MISSING_FLOAT; 
              }
	
              // Set sharp keys   
              drms_setkey_float(sharpoutrec, "MEANGBZ",  mean_derivative_bz); 
              drms_setkey_float(sharpoutrec, "ERRBZ",  mean_derivative_bz_err); 
	
              // Set sharp cea keys   
              drms_setkey_float(sharpceaoutrec, "MEANGBZ", mean_derivative_bz); 
              drms_setkey_float(sharpceaoutrec, "ERRBZ", mean_derivative_bz_err); 	
           }

           /***** MEANJZD and TOTUSJZ, Example Function 9 (Requires Function 8) ***/

           if (meanjzdflag || totusjzflag)
           {
              // First, compute Jz on all pixels
	      computeJz(bx_err, by_err, bx, by, dims, jz, jz_err, jz_err_squared, mask, bitmask, cdelt1, rsun_ref, rsun_obs, 
                     derx, dery);
   
              // Then take the sums of the appropriate pixels
              if(computeJzsmooth(bx, by, dims, jz, jz_smooth, jz_err, jz_rms_err, jz_err_squared_smooth, &mean_jz,
                       &mean_jz_err, &us_i, &us_i_err, mask, bitmask, cdelt1,
                       rsun_ref, rsun_obs, derx, dery))


              {
                 mean_jz                = DRMS_MISSING_FLOAT;
	         us_i                   = DRMS_MISSING_FLOAT;
                 mean_jz_err            = DRMS_MISSING_FLOAT;
                 us_i_err               = DRMS_MISSING_FLOAT;
	      }

              // Set sharp keys   
              drms_setkey_float(sharpoutrec, "MEANJZD",  mean_jz); 
              drms_setkey_float(sharpoutrec, "TOTUSJZ",  us_i); 	
              drms_setkey_float(sharpoutrec, "ERRJZ",  mean_jz_err); 
              drms_setkey_float(sharpoutrec, "ERRUSI",  us_i_err); 	

              // Set sharp cea keys   
              drms_setkey_float(sharpceaoutrec, "MEANJZD", mean_jz); 
              drms_setkey_float(sharpceaoutrec, "TOTUSJZ", us_i); 	
              drms_setkey_float(sharpceaoutrec, "ERRJZ", mean_jz_err); 
              drms_setkey_float(sharpceaoutrec, "ERRUSI", us_i_err); 	
           }

           /***** MEANALP, Example Function 10 (Requires Function 8)*********/

           if (meanalpflag)
           {
              // First, compute Jz on all pixels
   	      computeJz(bx_err, by_err, bx, by, dims, jz, jz_err, jz_err_squared, mask, bitmask, cdelt1, rsun_ref, rsun_obs, 
                        derx, dery);

              // Then compute alpha quantities on select pixels
	      if (computeAlpha(jz_err, bz_err, bz, dims, jz, jz_smooth, &mean_alpha, &mean_alpha_err, 
                               mask, bitmask, cdelt1, rsun_ref, rsun_obs))
              {
	         mean_alpha             = DRMS_MISSING_FLOAT;
                 mean_alpha_err         = DRMS_MISSING_FLOAT;
              }

              // Set sharp keys   	
              drms_setkey_float(sharpoutrec, "MEANALP", mean_alpha); 
              drms_setkey_float(sharpoutrec, "ERRALP", mean_alpha_err); 	

              // Set sharp cea keys   
              drms_setkey_float(sharpceaoutrec, "MEANALP", mean_alpha); 
              drms_setkey_float(sharpceaoutrec, "ERRALP", mean_alpha_err); 	

           }

           /***** MEANJZH, TOTUSJH, ABSNJZH, Example Function 11 (Requires Function 8) ***/

           if (meanjzhflag || totusjhflag || absnjzhflag)
           {
              // First, compute Jz on all pixels
    	      computeJz(bx_err, by_err, bx, by, dims, jz, jz_err, jz_err_squared, mask, bitmask, cdelt1, rsun_ref, rsun_obs, 
                        derx, dery);

              // Then compute helicity quantities on select pixels
           if (computeHelicity(jz_err, jz_rms_err, bz_err, bz, dims, jz, &mean_ih, &mean_ih_err, &total_us_ih, 
                               &total_abs_ih, &total_us_ih_err, &total_abs_ih_err, mask, bitmask, cdelt1, rsun_ref, rsun_obs))
              {  
	         mean_ih                = DRMS_MISSING_FLOAT; 
	         total_us_ih            = DRMS_MISSING_FLOAT;
  	         total_abs_ih           = DRMS_MISSING_FLOAT;
                 mean_ih_err            = DRMS_MISSING_FLOAT;
                 total_us_ih_err        = DRMS_MISSING_FLOAT;
                 total_abs_ih_err       = DRMS_MISSING_FLOAT;
	      }

              // Set sharp keys 
              drms_setkey_float(sharpoutrec, "MEANJZH",  mean_ih); 
              drms_setkey_float(sharpoutrec, "TOTUSJH",  total_us_ih); 	
              drms_setkey_float(sharpoutrec, "ABSNJZH",  total_abs_ih); 
              drms_setkey_float(sharpoutrec, "ERRMIH",  mean_ih_err); 	
              drms_setkey_float(sharpoutrec, "ERRTUI",  total_us_ih_err); 
              drms_setkey_float(sharpoutrec, "ERRTAI",  total_abs_ih_err); 	

              // Set sharp cea keys 
              drms_setkey_float(sharpceaoutrec, "MEANJZH", mean_ih); 
              drms_setkey_float(sharpceaoutrec, "TOTUSJH", total_us_ih); 	
              drms_setkey_float(sharpceaoutrec, "ABSNJZH", total_abs_ih); 
              drms_setkey_float(sharpceaoutrec, "ERRMIH", mean_ih_err); 	
              drms_setkey_float(sharpceaoutrec, "ERRTUI", total_us_ih_err); 
              drms_setkey_float(sharpceaoutrec, "ERRTAI", total_abs_ih_err); 	
           }

           /***** SAVNCPP, Example Function 12 (Requires Function 8) *******************/
   
           if (savncppflag)
           {
              // First, compute Jz on all pixels
  	      computeJz(bx_err, by_err, bx, by, dims, jz, jz_err, jz_err_squared, mask, bitmask, cdelt1, rsun_ref, rsun_obs, 
                        derx, dery);

              // Then compute sums of Jz on select pixels
  	      if (computeSumAbsPerPolarity(jz_err, bz_err, bz, jz, dims, &totaljz, &totaljz_err, 
				           mask, bitmask, cdelt1, rsun_ref, rsun_obs))
              {
	         totaljz                = DRMS_MISSING_FLOAT;
                 totaljz_err            = DRMS_MISSING_FLOAT;
	      }

              // Set sharp keys 	
              drms_setkey_float(sharpoutrec, "SAVNCPP",  totaljz); 
              drms_setkey_float(sharpoutrec, "ERRJHT",  totaljz_err); 	

              // Set sharp cea keys 
              drms_setkey_float(sharpceaoutrec, "SAVNCPP", totaljz); 
              drms_setkey_float(sharpceaoutrec, "ERRJHT", totaljz_err); 	
           }

           /***** MEANSHR and SHRGT45, Example Function 14 (Requires Keiji's code) **********/

           if (meanshrflag || shrgt45flag)
           {
              // First compute potential field
	      for (int i = 0; i < nxny; i++) bpz[i] = bz[i];
	      greenpot(bpx, bpy, bpz, nx, ny); 

              // Then compute shear angles
	      if (computeShearAngle(bx_err, by_err, bh_err, bx, by, bz, bpx, bpy, bpz, dims, 
		   		    &meanshear_angle, &meanshear_angle_err, &area_w_shear_gt_45,  
				    mask, bitmask)) 
              {
	          meanshear_angle    = DRMS_MISSING_FLOAT; // If fail, fill in NaN
		  area_w_shear_gt_45 = DRMS_MISSING_FLOAT;
                  meanshear_angle_err= DRMS_MISSING_FLOAT;
	      }

              // Set sharp keys 	
              drms_setkey_float(sharpoutrec, "MEANSHR",  meanshear_angle); 
              drms_setkey_float(sharpoutrec, "SHRGT45",  area_w_shear_gt_45); 	
              drms_setkey_float(sharpoutrec, "ERRMSHA",  meanshear_angle_err); 	

              // Set sharp cea keys 
              drms_setkey_float(sharpceaoutrec, "MEANSHR", meanshear_angle); 
              drms_setkey_float(sharpceaoutrec, "SHRGT45", area_w_shear_gt_45); 	
              drms_setkey_float(sharpceaoutrec, "ERRMSHA", meanshear_angle_err); 	  
           }

           /******************************* END FLAGS **********************************/

	   drms_free_array(bitmaskArray);		
	   drms_free_array(maskArray);
	   drms_free_array(bxArray);           
	   drms_free_array(byArray);
	   drms_free_array(bzArray);
	   drms_free_array(bx_errArray);           
	   drms_free_array(by_errArray);
	   drms_free_array(bz_errArray);
	} //endfor

	drms_close_records(sharpoutrecset, DRMS_INSERT_RECORD);
	drms_close_records(sharpceaoutrecset, DRMS_INSERT_RECORD);
   
        //used for select calculations
	free(bh); free(bt); free(jz); 
	free(bpx); free(bpy); free(bpz); 
	free(derx); free(dery);	 
	free(derx_bt); free(dery_bt); 	
	free(derx_bz); free(dery_bz);	 
	free(derx_bh); free(dery_bh); 
	free(bt_err); free(bh_err);  free(jz_err); 
        free(jz_err_squared); free(jz_rms_err);
	free(cvsinfoall);
 
return 0;    
}	// DoIt


