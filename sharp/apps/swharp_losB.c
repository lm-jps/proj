/*
 *  swharp_losB.
 *
 *  Created by Xudong Sun on 8/22/11.
 *  Modified to -- include ALL spaceweather keywords by Monica Bobra 25 Aug 2011
 *              -- include potential field calculation 
 *              -- run only on los data 5 March 2011
 *  Bz arrays
 *  Write out abs(B) as data segment and a few keywords as SW index
 *
 *  Use:
 *  First use the mmap module to create the in= and mask= parameters:
 *  /home/xudong/cvs/JSOC/bin/linux_x86_64/mmap "in=hmi.M_720s[2012.02.01_03:48:00_TAI]" "harp=su_turmon.Mharpv9_720s[1][2012.02.01_03:48:00_TAI]" "out=su_mbobra.test_mmap" "map=Postel" -a -e -v
 *  /home/xudong/cvs/JSOC/bin/linux_x86_64/mmap "harp=su_turmon.Mharpv9_720s[1][2012.02.01_03:48:00_TAI]" -z  "out=su_mbobra.test_mmap_bitmap" "map=Postel" "segment=bitmap" -a -v
 *
 *  then run this module:
 *  swharp_losB  "in=su_mbobra.test_mmap[][2011.10.06_23:59:60_TAI/1d]" /
 *  "mask=su_mbobra.test_mmap_bitmap[][2011.10.06_23:59:60_TAI/1d]" "out=su_mbobra.swharp_test_v1" "dzvalue=0.001"
 */


#include <jsoc_main.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sharp_functions.c"

#define DIE(msg) {fflush(stdout); fprintf(stderr, "%s, status=%d\n", msg, status); return(status);}
#define SHOW(msg) {printf("%s", msg); fflush(stdout);}
#define IN_FILES	3   /* Number of input files */
#define PI (3.141592653589793)
#define CMPERPIX (0.504277*696000000.0*100.)/(943.) /* cm/pixel = (CDELT1*RSUN_REF*100./RSUN_OBS) */

/* CMSQUARED = CMX*CMY
   CMX = CDELT1*RSUN_REF*(#of pixels in the x-direction)*100/RSUN_OBS
   CMY = CDELT2*RSUN_REF*(#of pixels in the y-direction)*100/RSUN_OBS */


/* declaring all the functions */

int computeAbsFlux(float *bz, int *dims, float *absFlux, float *mean_vf_ptr, int *mask);
int computeBh(float *bpx, float *bpy, float *bz, float *bh, int *dims, float *mean_hf_ptr, int *mask);
int computeGamma(float *bpx, float *bpy, float *bz, float *bh, int *dims, float *mean_gamma_ptr, int *mask);
int readFits(char *filename, float **image, int *dims);
int writeFits(char *filename, float *image, int *dims);
int computeB_total(float *bpx, float *bpy, float *bz, float *bt, int *dims, int *mask);
int computeBtotalderivative(float *bt, int *dims, float *mean_derivative_btotal_ptr, int *mask);
int computeBhderivative(float *bh, int *dims, float *mean_derivative_bh_ptr, int *mask);
int computeBzderivative(float *bz, int *dims, float *mean_derivative_bz_ptr, int *mask);
void greenpot(float *bx, float *by, float *bz, int nnx, int nny); 

char *module_name = "swharp_losB";	/* Module name */

ModuleArgs_t module_args[] =
{
    {ARG_STRING, "in", NULL, "Input vec mag recordset."},
    {ARG_STRING, "mask", NULL, "Input bitmap recordset."},
    {ARG_STRING, "out", NULL, "Output series."},
    {ARG_FLOAT,  "dzvalue", NULL, "Monopole depth."},
    {ARG_END}
};


int DoIt(void)
{
	
    int status = DRMS_SUCCESS;
	
	char *inQuery, *outQuery;		// input series query string
        char *maskQuery;		// mask series query string
	DRMS_RecordSet_t *inRecSet, *outRecSet, *maskRecSet;
	DRMS_Record_t *inRec, *outRec, *maskRec;
        DRMS_Segment_t *inSegBx, *inSegBy, *inSegBz, *outSeg, *maskSeg;
	DRMS_Array_t *inArrayBx, *inArrayBy, *inArrayBz, *outArray, *maskArray;	
        float *bx, *by, *bz, *outData, *bh, *bt, *jz, *bpx, *bpy, *bpz;
        int *mask;
        int dims[2], nxny, nx, ny;		// dimensions;  NAXIS1 = dims[0] which is the number of columns.
        float mean_vf = 0.0; 
        float absFlux = 0.0;
        float mean_hf = 0.0;
        float mean_gamma = 0.0;
        float mean_derivative_btotal = 0.0;
        float mean_derivative_bh = 0.0; 
        float mean_derivative_bz = 0.0;
        float mean_jz = 0.0;
        float us_i = 0.0;
        float mean_alpha = 0.0;
        float mean_ih = 0.0;
        float total_us_ih = 0.0;
        float total_abs_ih = 0.0;
        float totaljz = 0.0;
        float totpot  =0.0;
        float meanpot = 0.0;
        float area_w_shear_gt_45 = 0.0;
        float meanshear_angle = 0.0;
        float area_w_shear_gt_45h = 0.0;
        float meanshear_angleh = 0.0;
	int nrecs, irec, i;

	/* Input */
	
	inQuery = (char *) params_get_str(&cmdparams, "in");
	inRecSet = drms_open_records(drms_env, inQuery, &status);
	if (status || inRecSet->n == 0) DIE("No input data found");
        nrecs = inRecSet->n;
	
	/* Mask */
	
	maskQuery = (char *) params_get_str(&cmdparams, "mask");
	maskRecSet = drms_open_records(drms_env, maskQuery, &status);
	if (status || maskRecSet->n == 0) DIE("No mask data found");
        if (maskRecSet->n != nrecs) DIE("Mask and Input series do not have a 1:1 match"); 

	/* Output */
	
	outQuery = (char *) params_get_str(&cmdparams, "out");
	outRecSet = drms_create_records(drms_env, nrecs, outQuery, DRMS_PERMANENT, &status);
        if (status) DIE("Output recordset not created");
	
	/* Do this for each record */
	
        for (irec = 0; irec < nrecs; irec++)
        {
		
		/* Input record and data */
		
                inRec = inRecSet->records[irec];
		printf("Input Record #%d of #%d\n", irec+1, nrecs); fflush(stdout);
		
                maskRec = maskRecSet->records[irec];
		printf("Mask Record #%d of #%d\n", irec+1, nrecs); fflush(stdout);
		
		inSegBz = drms_segment_lookupnum(inRec, 0); /* Assume this is Bz equivalent */

		maskSeg = drms_segment_lookupnum(maskRec, 0); /* This is the bitmap */

		inArrayBz = drms_segment_read(inSegBz, DRMS_TYPE_FLOAT, &status);
		if (status) DIE("No Bz data file found. \n");
		
                maskArray = drms_segment_read(maskSeg, DRMS_TYPE_INT, &status);
		if (status) DIE("No mask data file found. \n");

                bz   = (float *)inArrayBz->data;
                mask = (int *)maskArray->data;

		nx = dims[0] = inArrayBz->axis[0];
		ny = dims[1] = inArrayBz->axis[1];
		nxny = dims[0] * dims[1];
                if (maskArray->axis[0] != nx || maskArray->axis[1] != ny) DIE("Mask and Input series are not of the same size"); 

                /* This is to modify the data for each PROJECTION method */
                int flag;
                char* value1;
                value1 = drms_getkey_string(inRec, "PROJECTION", &status);
                flag = strcmp("LambertCylindrical",value1);
                if (flag == 0)
                {
                      int i, j;   
                      for (j = 0; j < ny; j++) 
                      {
		         for (i = 0; i < nx; i++) 
		         {
		            by[j * nx + i] = - by[j * nx + i];
		         }	
                      }
		}
		
		/* Output data */

		outRec = outRecSet->records[irec];
                drms_setlink_static(outRec, "SRCLINK",  inRec->recnum);

                /*===========================================*/    
                /* Malloc some arrays     */

                bh   = (float *)malloc(nx*ny*sizeof(float));
                bt   = (float *)malloc(nx*ny*sizeof(float));
                jz   = (float *)malloc(nx*ny*sizeof(float));
                bpx  = (float *)malloc(nx*ny*sizeof(float));
                bpy  = (float *)malloc(nx*ny*sizeof(float));
                bpz  = (float *)malloc(nx*ny*sizeof(float));

                /*===========================================*/   
                /* SW Keyword computation */
                
                if (computeAbsFlux(bz, dims, &absFlux, &mean_vf, mask)) 
                {
                   absFlux = 0.0 / 0.0;		// If fail, fill in NaN
                   mean_vf = 0.0 / 0.0;
                }
		drms_setkey_float(outRec, "USFLUX", mean_vf);

                for (i=0 ;i<nxny; i++){bpz[i]=bz[i];}
                greenpot(bpx, bpy, bpz, nx, ny);
			
                computeBh(bpx, bpy, bz, bh, dims, &mean_hf, mask);	
               
                if (computeGamma(bpx, bpy, bz, bh, dims, &mean_gamma, mask)) mean_gamma = 0.0 / 0.0;	
                drms_setkey_float(outRec, "MEANGAM", mean_gamma);

                computeB_total(bpx, bpy, bz, bt, dims, mask);	
                
                if (computeBtotalderivative(bt, dims, &mean_derivative_btotal, mask)) mean_derivative_btotal = 0.0 / 0.0;
                drms_setkey_float(outRec, "MEANGBT", mean_derivative_btotal);

                if (computeBhderivative(bh, dims, &mean_derivative_bh, mask)) mean_derivative_bh = 0.0 / 0.0;
                drms_setkey_float(outRec, "MEANGBH", mean_derivative_bh);
                
                if (computeBzderivative(bz, dims, &mean_derivative_bz, mask)) mean_derivative_bz = 0.0 / 0.0; // If fail, fill in NaN
                drms_setkey_float(outRec, "MEANGBZ", mean_derivative_bz);

                /*===========================================*/  
		/* Set non-SW keywords */
		
		drms_copykey(outRec, inRec, "T_REC");
		drms_copykey(outRec, inRec, "HARPNUM");
		
		/* Clean up */
    		drms_free_array(inArrayBz);
		drms_free_array(maskArray);

    }

	drms_close_records(inRecSet, DRMS_FREE_RECORD);
        drms_close_records(outRecSet, DRMS_INSERT_RECORD);
	
    return 0;

}
