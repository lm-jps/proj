/*
 *  swharp_vectorB.
 *
 *  Created by Xudong Sun on 8/22/11.
 *  Modified to -- include ALL spaceweather keywords by Monica Bobra 25 Aug 2011
 *              -- include potential field calculation 
 *              -- run only on los data 5 March 2011
 *  Bz arrays
 *  Write out abs(B) as data segment and a few keywords as SW index
 *
 *  Use:
 *  First use the bmap module to create the in= and mask= parameters.
 *
 *  then run this module:
 *  swharp_vectorB  "in=su_mbobra.test_mmap_me[][]" /
 *  "mask=su_mbobra.test_mmap_bitmap_me[][]" "out=su_mbobra.swharp_test_v2" "dzvalue=0.001"
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

/* declaring all the functions */

int computeAbsFlux(float *bz, int *dims, float *absFlux, float *mean_vf_ptr, int *mask);
int computeBh(float *bx, float *by, float *bz, float *bh, int *dims, float *mean_hf_ptr, int *mask);
int computeGamma(float *bx, float *by, float *bz, float *bh, int *dims, float *mean_gamma_ptr, int *mask);
int readFits(char *filename, float **image, int *dims);
int writeFits(char *filename, float *image, int *dims);
int computeB_total(float *bx, float *by, float *bz, float *bt, int *dims, int *mask);
int computeBtotalderivative(float *bt, int *dims, float *mean_derivative_btotal_ptr, int *mask);
int computeBhderivative(float *bh, int *dims, float *mean_derivative_bh_ptr, int *mask);
int computeBzderivative(float *bz, int *dims, float *mean_derivative_bz_ptr, int *mask);
int computeJz(float *by, float *bx, int *dims, float *jz, float *mean_jz_ptr, float *us_i_ptr, int *mask);
int computeAlpha(float *bz, int *dims, float *jz, float *mean_alpha_ptr, int *mask);
int computeHelicity(float *bz, int *dims, float *jz, float *mean_ih_ptr, float *total_us_ih_ptr, float *total_abs_ih_ptr, int *mask);
int computeSumAbsPerPolarity(float *bz, float *jz, int *dims, float *totaljzptr, int *mask);
void greenpot(float *bx, float *by, float *bz, int nnx, int nny); 
int computeFreeEnergy(float *bx, float *by, float *bpx, float *bpy, int *dims, float *meanpotptr, float *totpotptr, int *mask);
int computeShearAngle(float *bx, float *by, float *bz, float *bpx, float *bpy, float *bpz, int *dims, float *meanshear_angleptr, float *area_w_shear_gt_45ptr, float *meanshear_anglehptr, float *area_w_shear_gt_45hptr, int *mask);

char *module_name = "swharp_vectorB";	/* Module name */

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
 
                inSegBx = drms_segment_lookupnum(inRec, 0); /* Assume this is Bx equivalent */
                inSegBy = drms_segment_lookupnum(inRec, 1); /* Assume this is By equivalent */
		inSegBz = drms_segment_lookupnum(inRec, 2); /* Assume this is Bz equivalent */

		maskSeg = drms_segment_lookupnum(maskRec, 0); /* This is the bitmap */

                inArrayBx = drms_segment_read(inSegBx, DRMS_TYPE_FLOAT, &status);
                if (status) DIE("No Bx data file found. \n");
                inArrayBy = drms_segment_read(inSegBy, DRMS_TYPE_FLOAT, &status);
                if (status) DIE("No By data file found. \n");
		inArrayBz = drms_segment_read(inSegBz, DRMS_TYPE_FLOAT, &status);
		if (status) DIE("No Bz data file found. \n");
		
                maskArray = drms_segment_read(maskSeg, DRMS_TYPE_INT, &status);
		if (status) DIE("No mask data file found. \n");

                bx   = (float *)inArrayBx->data;
                by   = (float *)inArrayBy->data;
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
			
                computeBh(bx, by, bz, bh, dims, &mean_hf, mask);	
               
                if (computeGamma(bx, by, bz, bh, dims, &mean_gamma, mask)) mean_gamma = 0.0 / 0.0;	
                drms_setkey_float(outRec, "MEANGAM", mean_gamma);

                computeB_total(bx, by, bz, bt, dims, mask);	
                
                if (computeBtotalderivative(bt, dims, &mean_derivative_btotal, mask)) mean_derivative_btotal = 0.0 / 0.0;
                drms_setkey_float(outRec, "MEANGBT", mean_derivative_btotal);

                if (computeBhderivative(bh, dims, &mean_derivative_bh, mask)) mean_derivative_bh = 0.0 / 0.0;
                drms_setkey_float(outRec, "MEANGBH", mean_derivative_bh);
                
                if (computeBzderivative(bz, dims, &mean_derivative_bz, mask)) mean_derivative_bz = 0.0 / 0.0; // If fail, fill in NaN
                drms_setkey_float(outRec, "MEANGBZ", mean_derivative_bz);

                if(computeJz(bx, by, dims, jz, &mean_jz, &us_i, mask))
                {  
                   mean_jz = 0.0 / 0.0;
                   us_i = 0.0 /0.0; 
                } 
                drms_setkey_float(outRec, "MEANJZD", mean_jz);
                drms_setkey_float(outRec, "TOTUSJZ", us_i);

                if (computeAlpha(bz, dims, jz, &mean_alpha, mask)) mean_alpha = 0.0 / 0.0; 
                drms_setkey_float(outRec, "MEANALP", mean_alpha);

                if (computeHelicity(bz, dims, jz, &mean_ih, &total_us_ih, &total_abs_ih, mask)) 
                {  
                   mean_ih     = 0.0/0.0;
                   total_us_ih = 0.0/0.0;
                   total_abs_ih= 0.0/0.0;
                } 
                drms_setkey_float(outRec, "MEANJZH", mean_ih);
                drms_setkey_float(outRec, "TOTUSJH", total_us_ih);
                drms_setkey_float(outRec, "ABSNJZH", total_abs_ih);

                if (computeSumAbsPerPolarity(bz, jz, dims, &totaljz, mask)) totaljz = 0.0 / 0.0;
                drms_setkey_float(outRec, "SAVNCPP", totaljz);

                if (computeFreeEnergy(bx, by, bpx, bpy, dims, &meanpot, &totpot, mask)) 
                {
                   meanpot = 0.0 / 0.0; // If fail, fill in NaN
                   totpot = 0.0 / 0.0;
                }
                drms_setkey_float(outRec, "MEANPOT", meanpot);
                drms_setkey_float(outRec, "TOTPOT", totpot);

                if (computeShearAngle(bx, by, bz, bpx, bpy, bpz, dims, &meanshear_angle, &area_w_shear_gt_45, &meanshear_angleh, &area_w_shear_gt_45h, mask)) 
                {
                   meanshear_angle = 0.0 / 0.0; // If fail, fill in NaN
                   area_w_shear_gt_45 = 0.0/0.0;
                   meanshear_angleh = 0.0 / 0.0; // If fail, fill in NaN
                   area_w_shear_gt_45h = 0.0/0.0;
                }
                printf("meanshear_angle=%f, area_w_shear_gt_45=%f, meanshear_angleh=%f, area_w_shear_gt_45h=%f\n",meanshear_angle,area_w_shear_gt_45,meanshear_angleh,area_w_shear_gt_45h);
                drms_setkey_float(outRec, "MEANSHR", meanshear_angle);
                drms_setkey_float(outRec, "SHRGT45", area_w_shear_gt_45);
                //drms_setkey_float(outRec, "MEANSHRH", meanshear_angleh);
                //drms_setkey_float(outRec, "SHRGT45H", area_w_shear_gt_45h);

 
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
