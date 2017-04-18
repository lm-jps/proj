/*
 *  qmap4pfss.c
 *
 *  This module computes squashing degree Q in a PFSS model
 *  on several 2D slices, i.e. pre-defined grid.
 *  The PFSS field vectors Br, Bt, Bp are read in as 3 data cubes.
 *  This module then calls a Fortran subroutine, mapfl_wrapper()
 *  that controls working Fortran functions. The subroutine
 *  returns Q-maps and coronal hole maps.
 *  
 *  Author:
 *      Xudong Sun, Zoran Mikic
 *
 *  Version:
 *      v0.0    Mar 29 2016
 *      v0.1    Apr 18 2017
 *
 *  Notes:
 *      v0.0
 *          Input/output grid hard-coded for now
 *      v0.1
 *          Changed output coronal hole to lowest layer only
 *          Store logQ now in compressed form
 *
 *
 *  Example calls:
 *      qmap4pfss "in=su_xudong.pfss[2099][0]" "out=su_xudong.qmap_pfss" -v -c "np=1441" "nt=721"

 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>
#include "jsoc_main.h"
#include <omp.h>
#include "radius.h"

#define PI              (M_PI)
#define TWOPI           (2.*M_PI)
#define RADSINDEG		(PI/180.)
#define RAD2ARCSEC		(648000./M_PI)
#define SECINDAY		(86400.)
#define FOURK			(4096)
#define FOURK2          (16777216)

#define DIE(msg) {fflush(stdout); fprintf(stderr,"%s, status=%d\n", msg, status); return(status);}
#define SHOW(msg) {printf("%s", msg); fflush(stdout);}

#define kNotSpecified "Not Specified"

#define ARRLENGTH(ARR) (sizeof(ARR) / sizeof(ARR[0]))
#define NR_OUT      (ARRLENGTH(r_out))

// FORTRAN function

extern void mapfl_wrapper_(double *bp, double *bt, double *br,
                           double *p, double *t, double *r,
                           int *np, int *nt, int *nr,
                           double *p_q, double *t_q, double *r_q,
                           int *np_q, int *nt_q, double *r_ch,
                           int *do_chmap, int *vb, int *do_cubic,
                           double *qmap, double *chmap);

//int ptest = 325, ttest = 363, rtest = 9;		// test pixel address of phi,theta, starting from 0

/* ========================================================================================================== */

char *module_name = "qmap4pfss";

ModuleArgs_t module_args[] =
{
    {ARG_STRING, "in", kNotSpecified, "Input PFSS series."},
    {ARG_STRING, "out", kNotSpecified, "Output Q-map series."},
    {ARG_FLAG, "v", "", "Verbose mode for Fortran."},
    {ARG_INT, "np", "1441", "Output map gridpoints in longitude."},
    {ARG_INT, "nt", "721", "Output map gridpoints in longitude."},
    {ARG_FLAG, "c", "", "Use cubic interpolation."},
    {ARG_END}
};

int DoIt(void)
{
    
    int status = DRMS_SUCCESS;

    /* Get data series */
    
    char *inQuery = (char *) params_get_str(&cmdparams, "in");
    char *outQuery = (char *) params_get_str(&cmdparams, "out");
    int vb = params_isflagset(&cmdparams, "v");
    int do_cubic = params_isflagset(&cmdparams, "c");
    int np_q = params_get_int(&cmdparams, "np");
    int nt_q = params_get_int(&cmdparams, "nt");
    
    /* Input Data */
    
    DRMS_RecordSet_t *inRS = drms_open_records(drms_env, inQuery, &status);
    int nrecs = inRS->n;
    if (status || nrecs < 1) {
        DIE("Input data series error");
    }
    
    /* Output */
    
    DRMS_RecordSet_t *outRS = drms_create_records(drms_env, nrecs, outQuery, DRMS_PERMANENT, &status);
    if (status) {
        DIE("Output data series error");
    }

    /* Loop */
    
    for (int irec = 0; irec < nrecs; irec++) {
        
        // Input/output
        
        DRMS_Record_t *inRec = inRS->records[irec];
        DRMS_Record_t *outRec = outRS->records[irec];
        
        // Read input data
        
        DRMS_Segment_t *inSeg_bp = drms_segment_lookup(inRec, "Bp");
        DRMS_Segment_t *inSeg_bt = drms_segment_lookup(inRec, "Bt");
        DRMS_Segment_t *inSeg_br = drms_segment_lookup(inRec, "Br");
        
        DRMS_Array_t *inArray_bp = drms_segment_read(inSeg_bp, DRMS_TYPE_DOUBLE, &status);
        if (status) { SHOW("input error"); continue; }
        DRMS_Array_t *inArray_bt = drms_segment_read(inSeg_bt, DRMS_TYPE_DOUBLE, &status);
        if (status) { SHOW("input error"); continue; }
        DRMS_Array_t *inArray_br = drms_segment_read(inSeg_br, DRMS_TYPE_DOUBLE, &status);
        if (status) { SHOW("input error"); continue; }
        
        int np = inArray_br->axis[0], nt = inArray_br->axis[1], nr = inArray_br->axis[2];       // no error check
        printf("np=%d, nt=%d, nr=%d\n", np, nt, nr);
        int npix = np * nt * nr;
        
        double *bp_jsoc = (double *) inArray_bp->data;
        double *bt_jsoc = (double *) inArray_bt->data;
        double *br_jsoc = (double *) inArray_br->data;
        
        double *bp = (double *) (malloc(npix * sizeof(double)));
        double *bt = (double *) (malloc(npix * sizeof(double)));
        double *br = (double *) (malloc(npix * sizeof(double)));
        
        // Coordinate swap
        int idx = 0;
        int npnt = np * nt;
        for (int idx_p = 0; idx_p < np; idx_p++) {
            for (int idx_t = 0; idx_t < nt; idx_t++) {
                for (int idx_r = 0; idx_r < nr; idx_r++) {
                    bp[idx] = bp_jsoc[idx_r*npnt+(nt-1-idx_t)*np+idx_p];
                    bt[idx] = bt_jsoc[idx_r*npnt+(nt-1-idx_t)*np+idx_p];
                    br[idx++] = br_jsoc[idx_r*npnt+(nt-1-idx_t)*np+idx_p];
                }
            }
        }
        
        // Test order of elements
        /*
        printf("%f,%f,%f,\n", br_jsoc[0], br_jsoc[1], br_jsoc[2]);      // along p
        printf("%f,%f,%f,\n", br_jsoc[0], br_jsoc[np], br_jsoc[2*np]);  // along t
        printf("%f,%f,%f,\n", br_jsoc[0], br_jsoc[np*nt], br_jsoc[np*nt*2]);    // along r
         */
        
        // Input grid
        
        double *p = (double *) (malloc(np * sizeof(double)));
        double *t = (double *) (malloc(nt * sizeof(double)));
        double *r = (double *) (malloc(nr * sizeof(double)));
        for (idx = 0; idx < np; idx++) { p[idx] = idx * TWOPI / (np - 1); }
        for (idx = 0; idx < nt; idx++) { t[idx] = idx * PI / (nt - 1); }
        for (idx = 0; idx < nr; idx++) { r[idx] = r_pfss[idx]; }
        
        /*
        idx = ttest * np_q + ptest;
        printf("ph=%19.16f\n",p_q[idx]);
        printf("th=%19.16f\n",t_q[idx]);
         */
        
        // Output grid (2D)
        
        int nr_q = NR_OUT;
        int npnt_q = np_q * nt_q;
        printf("np_q=%d, nt_q=%d, nr_q=%d, npnt_q=%d\n", np_q, nt_q, nr_q, npnt_q);
        double *p_q = (double *) (malloc(npnt_q * sizeof(double)));
        double *t_q = (double *) (malloc(npnt_q * sizeof(double)));
        double *r_q = (double *) (malloc(npnt_q * sizeof(double)));
        for (int idx_t = 0; idx_t < nt_q; idx_t++) {
            for (int idx_p = 0; idx_p < np_q; idx_p++) {
                idx = idx_t * np_q + idx_p;
                p_q[idx] = idx_p * TWOPI / (np_q - 1);
                t_q[idx] = idx_t * PI / (nt_q - 1);
            }
        }
        
        // Output data
        
        double *qmap_cube = (double *) (malloc(npnt_q * nr_q * sizeof(double)));
        double *chmap = (double *) (malloc(npnt_q * sizeof(double)));
        
        // Loop over slices

        for (int ir = 0; ir < nr_q; ir++) {
//        for (int ir = rtest; ir < rtest+1; ir++) {
            
            for (idx = 0; idx < npnt_q; idx++) {
                r_q[idx] = r_out[ir];
            }
            
//            printf("r=%19.16f\n",r_q[ttest * np_q + ptest]);
            
            
            double r_ch = r_out[ir];
            double *qmap = qmap_cube + ir * npnt_q;     // layer #ir
            
            int do_chmap = (ir == 0) ? 1 : 0;       // CH map for bottom layer only
            
            if (ir == 0) {
            	for (idx = 0; idx < npnt_q; idx++) {
            		chmap[idx] = 1;
            	}
            }
            
            // Call Fortran function
            // Output qmap and chmap readily rotated into C convention
            
            mapfl_wrapper_(bp, bt, br, p, t, r, &np, &nt, &nr,
                           p_q, t_q, r_q, &np_q, &nt_q, &r_ch,
                           &do_chmap, &vb, &do_cubic,
                           qmap, chmap);
            
            // test
//            printf("%f,%f,%f\n",qmap[0],qmap[1],qmap[2]);
//            printf("%f,%f,%f\n",qmap[0],qmap[np_q],qmap[np_q*2]);
            
        }
        
        // Write out
        
        int npix_q = np_q * nt_q * nr_q;
        for (int idx = 0; idx < npix_q; idx++) {
            qmap_cube[idx] = log10(qmap_cube[idx]);
        }
        DRMS_Segment_t *outSeg_q = drms_segment_lookup(outRec, "logq");
        int dims_q[3] = {np_q, nt_q, nr_q};
        DRMS_Array_t *outArray_q = drms_array_create(DRMS_TYPE_DOUBLE, 3, dims_q, qmap_cube, &status);
        for (idx = 0; idx < 3; idx++) {
            outSeg_q->axis[idx] = outArray_q->axis[idx];
        }
        outArray_q->israw = 0;		// always compressed
        outArray_q->bzero = outSeg_q->bzero;
        outArray_q->bscale = outSeg_q->bscale;
        status = drms_segment_write(outSeg_q, outArray_q, 0);
        
        DRMS_Segment_t *outSeg_ch = drms_segment_lookup(outRec, "chmap");
        int dims_ch[2] = {np_q, nt_q};
        DRMS_Array_t *outArray_ch = drms_array_create(DRMS_TYPE_DOUBLE, 2, dims_ch, chmap, &status);
        for (idx = 0; idx < 2; idx++) {
            outSeg_ch->axis[idx] = outArray_ch->axis[idx];
        }
        outArray_ch->israw = 0;		// always compressed
        outArray_ch->bzero = outSeg_ch->bzero;
        outArray_ch->bscale = outSeg_ch->bscale;
        status = drms_segment_write(outSeg_ch, outArray_ch, 0);
        
        // Keywords
        
        drms_copykeys(outRec, inRec, 0, 0);     // copy all keys
        drms_setkey_int(outRec, "CUBIC", do_cubic);
        drms_setkey_float(outRec, "R_CHMAP", r_out[0]);
        
        // Clean up
        drms_free_array(inArray_bp); drms_free_array(inArray_bt); drms_free_array(inArray_br);
        drms_free_array(outArray_q); drms_free_array(outArray_ch);
        free(bp); free(bt); free(br);
        free(p); free(t); free(r);
        free(p_q); free(t_q); free(r_q);
        
    }
    
    drms_close_records(inRS, DRMS_FREE_RECORD);
    drms_close_records(outRS, DRMS_INSERT_RECORD);
 
    //
    
    return DRMS_SUCCESS;
    
}
