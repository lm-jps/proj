/*
 * Module name:		cp_polarfield.c
 *
 * Description:
 *  This module copies su_xudong.polmf_f into hmi.meanpf_720s
 *  The original, double 180x4 segment is splitted into
 *  four 180-element, Rice compressed segments
 *
 * Written by:		Xudong Sun (xudongs@stanford.edu)
 *
 * Version:
 *  v1.0		Aug 27 2015
 *
 * Issues:
 *  v1.0    No error checking as of now
 *
 *
 * Example:
 *  cp_polarfield "in=su_xudong.polmf_f[2011.11.01]"
 *
 */

#include <jsoc_main.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* ############# Macros ############# */

#define SHOW(msg) {printf("%s", msg); fflush(stdout);}
#define DIE(msg) {fflush(stdout); fprintf(stderr,"%s (status=%d)\n", msg, status); return(status);}

#define FREE_ARR(arr) {if (arr) {free(arr);}}
#define DRMS_FREE_ARR(arr) {if (arr) {drms_free_array(arr);}}

/* ################################################# */
/* ################## Main Module ################## */
/* ################################################# */

char *module_name = "cp_polarfield";

ModuleArgs_t module_args[] =
{
    {ARG_STRING, "in", "", "Input query"},
    {ARG_STRING, "out", "hmi.meanpf_720s", "Output query"},
    {ARG_END}
};

int DoIt(void)
{
    
    int status = DRMS_SUCCESS;
    
    /* ============== */
    /* Get parameters */
    /* ============== */
    
    char *inQuery = (char *)params_get_str(&cmdparams, "in");
    char *outQuery = (char *)params_get_str(&cmdparams, "out");
    
    /* ============= */
    /* Input records */
    /* ============= */
    
    DRMS_RecordSet_t *inRS = NULL;
    inRS = drms_open_records(drms_env, inQuery, &status);
    if (status || inRS->n == 0) DIE("Input records error.");
    
    int nrecs = inRS->n;
    
    /* ============== */
    /* Output records */
    /* ============== */
    
    DRMS_RecordSet_t *outRS = drms_create_records(drms_env, nrecs, outQuery, DRMS_PERMANENT, &status);
    if (status) DIE("Output records not created.");
    
    /* ============ */
    /* Loop records */
    /* ============ */
    
    printf("Start, %d records in total\n", nrecs);
    
    for (int irec = 0; irec < nrecs; irec++) {
        
        // Input
        
        DRMS_Record_t *inRec = inRS->records[irec];
        TIME t_rec = drms_getkey_time(inRec, "T_REC", &status);
        
        char t_rec_str[100];
        sprint_time(t_rec_str, t_rec, "TAI", 0);
        printf("Record #%d, [%s]\n", irec, t_rec_str);

        DRMS_Segment_t *inSeg = drms_segment_lookupnum(inRec, 0);
        DRMS_Array_t *inArray = NULL;
        inArray = drms_segment_read(inSeg, DRMS_TYPE_DOUBLE, &status);
        if (status) {
            SHOW("Input array read error, record skipped.\n");
            if (inArray) drms_free_array(inArray);
            continue;
        }
        
        int nx = inArray->axis[0], ny = inArray->axis[1];       // 180x4
        int outDims[1] = {nx};
        
        double *inData = (double *) (inArray->data);
        double *mf_br = (double *) (malloc(nx * sizeof(double)));
        double *mf_bl = (double *) (malloc(nx * sizeof(double)));
        double *w = (double *) (malloc(nx * sizeof(double)));
        double *num = (double *) (malloc(nx * sizeof(double)));
        for (int i = 0; i < nx; i++) {
        	mf_br[i] = inData[i];
        	mf_bl[i] = inData[nx + i];
        	w[i] = inData[nx * 2 + i];
        	num[i] = inData[nx * 3 + i];
        }
        
        // Output
        
        DRMS_Record_t *outRec = outRS->records[irec];
        
        status = drms_copykeys(outRec, inRec, 0, 0);
        
        DRMS_Segment_t *outSeg_br = drms_segment_lookup(outRec, "mf_br");
        DRMS_Array_t *outArray_br = drms_array_create(DRMS_TYPE_DOUBLE, 1, outDims, mf_br, &status);
        status = drms_segment_write(outSeg_br, outArray_br, 0);
        
        DRMS_Segment_t *outSeg_bl = drms_segment_lookup(outRec, "mf_bl");
        DRMS_Array_t *outArray_bl = drms_array_create(DRMS_TYPE_DOUBLE, 1, outDims, mf_bl, &status);
        status = drms_segment_write(outSeg_bl, outArray_bl, 0);
        
        DRMS_Segment_t *outSeg_w = drms_segment_lookup(outRec, "w");
        DRMS_Array_t *outArray_w = drms_array_create(DRMS_TYPE_DOUBLE, 1, outDims, w, &status);
        status = drms_segment_write(outSeg_w, outArray_w, 0);
        
        DRMS_Segment_t *outSeg_num = drms_segment_lookup(outRec, "num");
        DRMS_Array_t *outArray_num = drms_array_create(DRMS_TYPE_DOUBLE, 1, outDims, num, &status);
        status = drms_segment_write(outSeg_num, outArray_num, 0);
        
        // Clean up
        DRMS_FREE_ARR(inArray);
        DRMS_FREE_ARR(outArray_br); DRMS_FREE_ARR(outArray_bl);
        DRMS_FREE_ARR(outArray_w); DRMS_FREE_ARR(outArray_num);
        
    }   // irec
    
    /* ============= */
    
    drms_close_records(inRS, DRMS_FREE_RECORD);
    drms_close_records(outRS, DRMS_INSERT_RECORD);
    
    return DRMS_SUCCESS;
    
}   // DoIt
