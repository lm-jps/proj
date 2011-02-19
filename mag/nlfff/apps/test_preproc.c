/*
 * Module name:		test_localpf.c
 *
 * Description:		Local PF using Wiegelmann's Green function code
 *
 * Calling:		
 *
 * Original source:	C NLFFF model by Thomas Wiegelmann (weigelmann@linmpi.mg.de)
 *
 * Written by:		Xudong Sun (xudongs@stanford.edu)
 *
 * Version:
 *			v1.0		Oct 13 2010
 *
 * Issues:
 *			v1.0
 *
 * Example:
 *			test_preproc "in=su_xudong.nlfff_test_in[2010]" "out=su_xudong.nlfff_preproc" "test=1"
 *			test_preproc "in=su_xudong.mercator_vec_S5[2010.07.01_11:48:00_TAI]" "out=su_xudong.nlfff_preproc" "test=0"
 */


#include <jsoc_main.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>

#define PI	(M_PI)
#define	DTOR	(PI / 180.)

#define ARRLENGTH(ARR) (sizeof(ARR) / sizeof(ARR[0]))
#define DIE(msg) {fflush(stdout); fprintf(stderr, "%s, status=%d\n", msg, status); return(status);}
#define SHOW(msg) {printf("%s", msg); fflush(stdout);}

#define Macro
#define MCENTERGRAD(f,id) ((f[i+id]-f[i-id])/(2*h))
#define MLEFTGRAD(f,id)   ((-3*f[i]+4*f[i+id]-f[i+2*id])/(2*h))
#define MRIGHTGRAD(f,id)  ((+3*f[i]-4*f[i-id]+f[i-2*id])/(2*h))
#define GRADX(f,i) ((ix>0 && ix<nx-1) ? (MCENTERGRAD(f,nynz)) : ((ix==0) ? (MLEFTGRAD(f,nynz)) : ((ix==nx-1) ? (MRIGHTGRAD(f,nynz)) : (0.0))))
#define GRADY(f,i) ((iy>0 && iy<ny-1) ? (MCENTERGRAD(f,nz)) : ((iy==0) ? (MLEFTGRAD(f,nz)) : ((iy==ny-1) ? (MRIGHTGRAD(f,nz)) : (0.0))))
#define GRADZ(f,i) ((iz>0 && iz<nz-1) ? (MCENTERGRAD(f,1)) : ((iz==0) ? (MLEFTGRAD(f,1)) : ((iz==nz-1) ? (MRIGHTGRAD(f,1)) : (0.0))))

// Working part
// ======== OpenMP =========
#ifdef _OPENMP 
#include <omp.h>
#endif 

#include "rebin.c"
#include "preproc.c"
#include "green.c"
#include "relax.c"


/* Timing by T. P. Larson */

double getwalltime(void)
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec * 1000.0 + tv.tv_usec/1000.0;
}

double getcputime(double *utime, double *stime)
{
	
	struct rusage ru;
	getrusage(RUSAGE_SELF, &ru);
	*utime = ru.ru_utime.tv_sec * 1000.0 + ru.ru_utime.tv_usec/1000.0;
	*stime = ru.ru_stime.tv_sec * 1000.0 + ru.ru_stime.tv_usec/1000.0;
	return *utime + *stime;
}




char *module_name = "test_preproc";	/* Module name */

ModuleArgs_t module_args[] =
{
    {ARG_STRING, "in", NULL, "Input data series."},
    {ARG_STRING, "out", NULL,  "Output data series."},
    {ARG_DOUBLE, "mu1", "1.0", "For vector magnetogram preprocessing."},
    {ARG_DOUBLE, "mu2", "1.0", "For vector magnetogram preprocessing."},
    {ARG_DOUBLE, "mu3", "0.001", "For vector magnetogram preprocessing."},
    {ARG_DOUBLE, "mu4", "0.01", "For vector magnetogram preprocessing."},
    {ARG_INT, "maxit", "10000", "Maximum itertation number."},
    {ARG_INT, "test", "0", "1 for test data with 3 layer input cube rather than 3 images."},
    {ARG_INT, "VERB", "2", "Level of verbosity: 0=errors/warnings; 1=status; 2=all"},
    {ARG_STRING, "COMMENT", " ", "Comments."}, 
    {ARG_END}
};


/* ################## Main Module ################## */

int DoIt(void)
{
    int status = DRMS_SUCCESS;
	
#ifdef _OPENMP
    printf("Compiled by OpenMP\n");
#endif
	
    char *inQuery, *outQuery;
    DRMS_RecordSet_t *inRS, *outRS;
    int irec, nrecs;
    DRMS_Record_t *inRec, *outRec;
    DRMS_Segment_t *inSegBx, *inSegBy, *inSegBz;
    DRMS_Segment_t *outSegBx, *outSegBy, *outSegBz;
    DRMS_Array_t *inArrayBx, *inArrayBy, *inArrayBz;
    DRMS_Array_t *outArrayBx, *outArrayBy, *outArrayBz;
    float *inDataBx, *inDataBy, *inDataBz;
    float *outDataBx, *outDataBy, *outDataBz;
    
    double *bx0, *by0, *bz0;
    double *bx0_now, *by0_now, *bz0_now;
    double *Bx, *By, *Bz, *Bx_tmp, *By_tmp, *Bz_tmp;
	
    int verbflag, prepflag;
    int outDims[3];
	
    int i, j, k, itest, dpt, dpt0;
    int nx, ny, nz, nxny, nynz, nxnynz;
    int nd, nd_now, maxit;
    double mu1, mu2, mu3, mu4;
    int multi, multi_now;
    int nx_now, ny_now, nz_now, nxny_now, nynz_now, nxnynz_now;
    int test;
    char *comment;
	
    // Time measuring
    double wt0, wt1, wt;
    double ut0, ut1, ut;
    double st0, st1, st;
    double ct0, ct1, ct;
    wt0 = getwalltime();
    ct0 = getcputime(&ut0, &st0);
	
    /* Get parameters */
    inQuery = (char *)params_get_str(&cmdparams, "in");
    outQuery = (char *)params_get_str(&cmdparams, "out");
    comment = (char *)params_get_str(&cmdparams, "COMMENT");
    verbflag = params_get_int(&cmdparams, "VERB");
    maxit = params_get_int(&cmdparams, "maxit");
    mu1 = params_get_double(&cmdparams, "mu1");
    mu2 = params_get_double(&cmdparams, "mu2");
    mu3 = params_get_double(&cmdparams, "mu3");
    mu4 = params_get_double(&cmdparams, "mu4");
    test = params_get_int(&cmdparams, "test");
    multi = 1;
	
    /* Open input */
    inRS = drms_open_records(drms_env, inQuery, &status);
    if (status || inRS->n == 0) DIE("No input data found");
    nrecs = inRS->n;
	
    /* Create output */
    outRS = drms_create_records(drms_env, nrecs, outQuery, DRMS_PERMANENT, &status);
    if (status) DIE("Output recordset not created");
	
    /* Do this for each record */
    for (irec = 0; irec < nrecs; irec++)
    {
		
        if (verbflag) {
            wt1 = getwalltime();
            ct1 = getcputime(&ut1, &st1);
            printf("processing record %d...\n", irec);
        }
        
        /* Input record and data */
        inRec = inRS->records[irec];
        if (!test) {
            inSegBx = drms_segment_lookup(inRec, "Bx");
            inSegBy = drms_segment_lookup(inRec, "By");
            inSegBz = drms_segment_lookup(inRec, "Bz");
            inArrayBx = drms_segment_read(inSegBx, DRMS_TYPE_FLOAT, &status);
            if (status) DIE("No Bx data file found. \n");
            inArrayBy = drms_segment_read(inSegBy, DRMS_TYPE_FLOAT, &status);
            if (status) DIE("No By data file found. \n");
            inArrayBz = drms_segment_read(inSegBz, DRMS_TYPE_FLOAT, &status);
            if (status) DIE("No Bz data file found. \n");
            inDataBx = (float *)inArrayBx->data;
            inDataBy = (float *)inArrayBy->data;
            inDataBz = (float *)inArrayBz->data;
            nx = inArrayBx->axis[0]; ny = inArrayBx->axis[1];
            if (inArrayBy->axis[0] != nx || inArrayBz->axis[0] != nx || 
                inArrayBy->axis[1] != ny || inArrayBz->axis[1] != ny)
                DIE("Dimension error. \n");
            nxny = nx * ny; // nynz = ny * nz; nxnynz = nxny * nz;
        } else {
			inSegBz = drms_segment_lookupnum(inRec, 0);
			inArrayBz = drms_segment_read(inSegBz, DRMS_TYPE_FLOAT, &status);
			if (status) DIE("No data file found. \n");
			if (inArrayBz->naxis != 3) DIE("Wrong data dimension. \n");
			nx = inArrayBz->axis[0]; ny = inArrayBz->axis[1];
			nxny = nx * ny; // nynz = ny * nz; nxnynz = nxny * nz;
			inDataBx = (float *)inArrayBz->data;
			inDataBy = inDataBx + nxny;
			inDataBz = inDataBx + 2 * nxny;
        }
		
        /* Output data */
        outDims[0] = nx; outDims[1] = ny;
		
        outArrayBx = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, NULL, &status);
        outArrayBy = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, NULL, &status);
        outArrayBz = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, NULL, &status);
        outDataBx = (float *)outArrayBx->data;
        outDataBy = (float *)outArrayBy->data;
        outDataBz = (float *)outArrayBz->data;
        
		
        /* ======================== */
        /* This is the working part.*/
        /* ======================== */
		
        // Copy in magnetogram--------------- CHECKED AGAINST WIEGELMANN'S CODE, COLUMN FIRST
		Bx = (double *)(calloc(nxny, sizeof(double)));
        By = (double *)(calloc(nxny, sizeof(double)));
        Bz = (double *)(calloc(nxny, sizeof(double)));
        dpt = 0;
        for (i = 0; i < nx; i++)	// See green()
			for (j = 0; j < ny; j++) {
				dpt0 = j * nx + i;
				Bx[dpt] = inDataBx[dpt0];
				By[dpt] = inDataBy[dpt0];
				Bz[dpt] = inDataBz[dpt0];
				dpt++;
			}
		
		// Main
		preproc(Bx, By, Bz, nx, ny, 
				mu1, mu2, mu3, mu4, 0.0, maxit, verbflag);
		
		// Copy out result -------------- CHECKED AGAINST WIEGELMANN'S CODE
		for (j = 0; j < ny; j++) {
		for (i = 0; i < nx; i++) {
			outDataBx[j * nx + i] = Bx[i * ny + j];
			outDataBy[j * nx + i] = By[i * ny + j];
			outDataBz[j * nx + i] = Bz[i * ny + j];
		}}
		
        /* Output record */
        outRec = outRS->records[irec];
        outSegBx = drms_segment_lookupnum(outRec, 0);
        outSegBy = drms_segment_lookupnum(outRec, 1);
        outSegBz = drms_segment_lookupnum(outRec, 2);
		
        for (i = 0; i < 2; i++) {	// For variable dimensions
            outSegBx->axis[i] = outArrayBx->axis[i];
            outSegBy->axis[i] = outArrayBy->axis[i];
            outSegBz->axis[i] = outArrayBz->axis[i];
        }
		
        outArrayBx->parent_segment = outSegBx;
        outArrayBy->parent_segment = outSegBy;
        outArrayBz->parent_segment = outSegBz;
		
        /* Result writing */
        status = drms_segment_write(outSegBx, outArrayBx, 0);
        if (status) DIE("Problem writing file Bx");
        drms_free_array(outArrayBx);
        status = drms_segment_write(outSegBy, outArrayBy, 0);
        if (status) DIE("Problem writing file By");
        drms_free_array(outArrayBy);
        status = drms_segment_write(outSegBz, outArrayBz, 0);
        if (status) DIE("Problem writing file Bz");
        drms_free_array(outArrayBz);
		
        /* Set keywords */
    	drms_setkey_string(outRec, "BLD_VERS", jsoc_version);
        drms_setkey_time(outRec, "DATE", CURRENT_SYSTEM_TIME);
        if (test) {
            drms_copykey(outRec, inRec, "NUM");
        } else {
            drms_copykey(outRec, inRec, "T_REC");
            drms_copykey(outRec, inRec, "T_OBS");
            drms_copykey(outRec, inRec, "PNUM");
            drms_setkey_string(outRec, "COMMENT", comment);
        }
		
		drms_setkey_int(outRec, "NUM", 100);
		
        /* Clean up */
        if (!test) {
            drms_free_array(inArrayBx);
            drms_free_array(inArrayBy);
            drms_free_array(inArrayBz);
        } else {
            drms_free_array(inArrayBz);
        }
		free(Bx); free(By); free(Bz);
		
        /* Time measure */
        if (verbflag) {
            wt = getwalltime();
            ct = getcputime(&ut, &st);
            printf("record %d done, %.2f ms wall time, %.2f ms cpu time\n", 
				   irec, wt - wt1, ct - ct1);
        }
		
    }
	
    drms_close_records(inRS, DRMS_FREE_RECORD);
    drms_close_records(outRS, DRMS_INSERT_RECORD);
	
    if (verbflag) {
        wt = getwalltime();
        ct = getcputime(&ut, &st);
        printf("total time spent: %.2f ms wall time, %.2f ms cpu time\n", 
			   wt - wt0, ct - ct0);
    }
	
    return(DRMS_SUCCESS);
}

