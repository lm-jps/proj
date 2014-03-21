/*
 *  sharp.c
 *
 *	This module computes the Lorentz force from vector field maps
 *  Outputs include three images Fx, Fy, Fz
 *  and the sum over the whole image
 *  as well as the sum within the active pixels
 *
 *	Author:
 *		Xudong Sun
 *
 *	Version:
 *		v0.0	Aug 15 2013
 *
 *	Notes:
 *		v0.0
 *
 *	Example:
 *	lorentz "in=hmi.sharp_cea_720s[377][2011.02.15_00:00]" "out=su_xudong.lorentz"
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>
#include "jsoc_main.h"
#include "astro.h"
#include "copy_me_keys.c"

#define ARRLENGTH(ARR) (sizeof(ARR) / sizeof(ARR[0]))

#define DIE(msg) {fflush(stdout); fprintf(stderr,"%s, status=%d\n", msg, status); return(status);}
#define SHOW(msg) {printf("%s", msg); fflush(stdout);}

#define kNotSpecified "Not Specified"

#define PIXSIZE ((double)3.64e7)

/* Data segment names, in Bx, By, Bz order */
char *vectorSegs[] = {"Bp", "Bt", "Br"};

/* Multiplier */
float vectorMulti[] = {1.,-1.,1.};

/* Lorentz segment names */
char *lorentzSegs[] = {"Fx", "Fy", "Fz"};

/* Mask name, thresholding */
char *maskSegs[] = {"bitmap", "conf_disambig"};
int maskThresh[] = {32, 61};

/* Set all keywords, no error checking for now */
void setKeys(DRMS_Record_t *outRec, DRMS_Record_t *inRec);

/* ====================================================== */

char *module_name = "lorentz";

ModuleArgs_t module_args[] =
{
	{ARG_STRING, "in", kNotSpecified, "Input sharp series."},
	{ARG_STRING, "out", kNotSpecified, "Output Lorentz force series."},
    {ARG_INT, "mask", "1", "Mask for weak field: 0 for bitmap, 1 for conf_disambig"},
	{ARG_END}
};

int DoIt(void)
{

    int status = DRMS_SUCCESS;
	int nrecs, irec;
    
    char *inQuery, *outQuery;
    DRMS_RecordSet_t *inRS = NULL, *outRS = NULL;       // Mar 20 2014
    
    /* Get parameters */
    
	inQuery = (char *) params_get_str(&cmdparams, "in");
	outQuery = (char *) params_get_str(&cmdparams, "out");
    int mask_id = params_get_int(&cmdparams, "mask");
    mask_id = mask_id ? 1 : 0;
    
    char *maskSegName = maskSegs[mask_id];
    int thresh = maskThresh[mask_id];
    printf("mask=%d, name=%s, thresh=%d\n", mask_id, maskSegName, thresh);
    
    /* Open input data */
    
    inRS = drms_open_records(drms_env, inQuery, &status);
    if (status || inRS->n == 0) DIE("No input data found");
    nrecs = inRS->n;
    
    /* Output data */

    outRS = drms_create_records(drms_env, nrecs, outQuery, DRMS_PERMANENT, &status);
    if (status) {
        DIE("Output error");
    }
    
    /* Loop through each record */
	
	for (irec = 0; irec < nrecs; irec++) {
        
        /* Info */
        
        DRMS_Record_t *inRec = inRS->records[irec];
        
        TIME t_rec = drms_getkey_time(inRec, "T_REC", &status);
        int harpnum = drms_getkey_time(inRec, "HARPNUM", &status);
        
        char *inRecInfo = (char *) malloc(100 * sizeof(char));
        char *t_rec_str = (char *) malloc(30 * sizeof(char));
        sprint_time(t_rec_str, t_rec, "TAI", 0);
        sprintf(inRecInfo, "%s[%d][%s]\n", inRec->seriesinfo->seriesname, harpnum, t_rec_str);
        
        printf("Record #%d of %d, %s\n", irec+1, nrecs, inRecInfo); fflush(stdout);
        
        /* Output record */
        
        DRMS_Record_t *outRec = outRS->records[irec];
        
        /* Input arrays */
        
        DRMS_Segment_t *maskSeg = drms_segment_lookup(inRec, maskSegName);
        DRMS_Array_t *maskArray = drms_segment_read(maskSeg, DRMS_TYPE_CHAR, &status);
        if (status) {
            SHOW("Error while reading mask, skip record\n");
            drms_free_array(maskArray);
        }
        
        DRMS_Segment_t *bxSeg = drms_segment_lookup(inRec, vectorSegs[0]);
        DRMS_Array_t *bxArray = drms_segment_read(bxSeg, DRMS_TYPE_FLOAT, &status);
        if (status) {
            SHOW("Error while reading Bx, skip record\n");
            drms_free_array(maskArray);
            drms_free_array(bxArray);
        }
        
        DRMS_Segment_t *bySeg = drms_segment_lookup(inRec, vectorSegs[1]);
        DRMS_Array_t *byArray = drms_segment_read(bySeg, DRMS_TYPE_FLOAT, &status);
        if (status) {
            SHOW("Error while reading By, skip record\n");
            drms_free_array(maskArray);
            drms_free_array(bxArray); drms_free_array(byArray);
        }
        
        DRMS_Segment_t *bzSeg = drms_segment_lookup(inRec, vectorSegs[2]);
        DRMS_Array_t *bzArray = drms_segment_read(bzSeg, DRMS_TYPE_FLOAT, &status);
        if (status) {
            SHOW("Error while reading Bz, skip record\n");
            drms_free_array(maskArray);
            drms_free_array(bxArray); drms_free_array(byArray); drms_free_array(bzArray);
        }
        
        int nx = bzArray->axis[0], ny = bzArray->axis[1];
        int nxny = nx * ny;
        
        char *mask = (char *) maskArray->data;
        float *bx = (float *) bxArray->data;
        float *by = (float *) byArray->data;
        float *bz = (float *) bzArray->data;
        
        /* Output arrays */
        
        int outDims[2] = {nx, ny};
        
        DRMS_Array_t *fxArray = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, NULL, &status);
        if (status) {
            SHOW("Error while creating Fx, skip record\n");
            drms_free_array(fxArray);
            drms_free_array(maskArray);
            drms_free_array(bxArray); drms_free_array(byArray); drms_free_array(bzArray);
        }

        DRMS_Array_t *fyArray = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, NULL, &status);
        if (status) {
            SHOW("Error while creating Fy, skip record\n");
            drms_free_array(fxArray); drms_free_array(fyArray);
            drms_free_array(maskArray);
            drms_free_array(bxArray); drms_free_array(byArray); drms_free_array(bzArray);
        }
        
        DRMS_Array_t *fzArray = drms_array_create(DRMS_TYPE_FLOAT, 2, outDims, NULL, &status);
        if (status) {
            SHOW("Error while creating Fz, skip record\n");
            drms_free_array(fxArray); drms_free_array(fyArray); drms_free_array(fzArray);
            drms_free_array(maskArray);
            drms_free_array(bxArray); drms_free_array(byArray); drms_free_array(bzArray);
        }
        
        float *fx = (float *) fxArray->data;
        float *fy = (float *) fyArray->data;
        float *fz = (float *) fzArray->data;
        
        /* Do the work */
        
        double totfx = 0, totfy = 0, totfz = 0;
        double totfx1 = 0, totfy1 = 0, totfz1 = 0;
        double bsq = 0, totbsq = 0, totbsq1 = 0;
        double epsx = 0, epsy = 0, epsz = 0;
        double epsx1 = 0, epsy1 = 0, epsz1 = 0;
        
        double area = pow(PIXSIZE, 2.0);
        double k_h = -1.0 * area / (4. * M_PI) / 1.0e20;
        double k_z = area / (8. * M_PI) / 1.0e20;
        
        for (int i = 0; i < nxny; i++) {
            
            fx[i] = (bx[i] * vectorMulti[0]) * (bz[i] * vectorMulti[2]) * k_h;
            fy[i] = (by[i] * vectorMulti[1]) * (bz[i] * vectorMulti[2]) * k_h;
            fz[i] = (bx[i] * bx[i] + by[i] * by[i] - bz[i] * bz[i]) * k_z;
            bsq = bx[i] * bx[i] + by[i] * by[i] + bz[i] * bz[i];
            
            totfx += fx[i]; totfy += fy[i]; totfz += fz[i];
            totbsq += bsq;
            
            if (mask[i] > thresh) {
                totfx1 += fx[i]; totfy1 += fy[i]; totfz1 += fz[i];
                totbsq1 += bsq;
            }
            
        }
        
        epsx = (totfx / k_h) / totbsq;
        epsy = (totfy / k_h) / totbsq;
        epsz = (totfz / k_z) / totbsq;
        epsx1 = (totfx1 / k_h) / totbsq1;
        epsy1 = (totfy1 / k_h) / totbsq1;
        epsz1 = (totfz1 / k_z) / totbsq1;
        
        /* Write image */
        
        DRMS_Segment_t *fxSeg = drms_segment_lookup(outRec, lorentzSegs[0]);
        DRMS_Segment_t *fySeg = drms_segment_lookup(outRec, lorentzSegs[1]);
        DRMS_Segment_t *fzSeg = drms_segment_lookup(outRec, lorentzSegs[2]);
        
        fxSeg->axis[0] = outDims[0]; fxSeg->axis[1] = outDims[1];
        fySeg->axis[0] = outDims[0]; fySeg->axis[1] = outDims[1];
        fzSeg->axis[0] = outDims[0]; fzSeg->axis[1] = outDims[1];
        
        fxArray->israw = 0;
        fyArray->israw = 0;
        fzArray->israw = 0;
        
        fxArray->bzero = fxSeg->bzero; fxArray->bscale = fxSeg->bscale;
        fyArray->bzero = fySeg->bzero; fyArray->bscale = fySeg->bscale;
        fzArray->bzero = fzSeg->bzero; fzArray->bscale = fzSeg->bscale;
        
        status = drms_segment_write(fxSeg, fxArray, 0);
        if (status) {
            SHOW("Error while writing Fx, skip record\n");
            drms_free_array(fxArray); drms_free_array(fyArray); drms_free_array(fzArray);
            drms_free_array(maskArray);
            drms_free_array(bxArray); drms_free_array(byArray); drms_free_array(bzArray);
        }
        
        status = drms_segment_write(fySeg, fyArray, 0);
        if (status) {
            SHOW("Error while writing Fy, skip record\n");
            drms_free_array(fxArray); drms_free_array(fyArray); drms_free_array(fzArray);
            drms_free_array(maskArray);
            drms_free_array(bxArray); drms_free_array(byArray); drms_free_array(bzArray);
        }
        
        status = drms_segment_write(fzSeg, fzArray, 0);
        if (status) {
            SHOW("Error while writing Fz, skip record\n");
            drms_free_array(fxArray); drms_free_array(fyArray); drms_free_array(fzArray);
            drms_free_array(maskArray);
            drms_free_array(bxArray); drms_free_array(byArray); drms_free_array(bzArray);
        }
        
        /* Links */
        
        DRMS_Link_t *sharpLink = hcon_lookup_lower(&outRec->links, "SHARP");
        if (sharpLink) drms_link_set("SHARP", outRec, inRec);
        
        /* Keywords */
        
        setKeys(outRec, inRec);     // All other keys
        
        drms_copykey(outRec, inRec, "T_REC");       // Prime keys
        drms_copykey(outRec, inRec, "HARPNUM");
        
        drms_setkey_double(outRec, "TOTFX", totfx); // indices
        drms_setkey_double(outRec, "TOTFY", totfy);
        drms_setkey_double(outRec, "TOTFZ", totfz);
        drms_setkey_double(outRec, "TOTBSQ", totbsq);
        drms_setkey_double(outRec, "TOTFX1", totfx1);
        drms_setkey_double(outRec, "TOTFY1", totfy1);
        drms_setkey_double(outRec, "TOTFZ1", totfz1);
        drms_setkey_double(outRec, "TOTBSQ1", totbsq1);
        drms_setkey_double(outRec, "EPSX", epsx);
        drms_setkey_double(outRec, "EPSY", epsy);
        drms_setkey_double(outRec, "EPSZ", epsz);
        drms_setkey_double(outRec, "EPSX1", epsx1);
        drms_setkey_double(outRec, "EPSY1", epsy1);
        drms_setkey_double(outRec, "EPSZ1", epsz1);
        drms_setkey_double(outRec, "AREA", area);
        drms_setkey_string(outRec, "MASKNAME", maskSegName);        // Feb 14
        drms_setkey_int(outRec, "MASKTHRS", thresh);        
        
        /* Clean up */
        
        drms_free_array(fxArray); drms_free_array(fyArray); drms_free_array(fzArray);
        drms_free_array(maskArray);
        drms_free_array(bxArray); drms_free_array(byArray); drms_free_array(bzArray);
        free(inRecInfo); free(t_rec_str);

    }   // irec
    
    drms_close_records(outRS, DRMS_INSERT_RECORD);
    drms_close_records(inRS, DRMS_FREE_RECORD);
    
    return DRMS_SUCCESS;
    
}   // DoIt



/*
 * Set all keywords, no error checking for now
 *
 */

void setKeys(DRMS_Record_t *outRec, DRMS_Record_t *inRec)
{

	copy_me_keys(inRec, outRec);
	copy_patch_keys(inRec, outRec);
	copy_geo_keys(inRec, outRec);
	copy_ambig_keys(inRec, outRec);
	
	drms_setkey_string(outRec, "BUNIT_000", "1e20 dyne");
	drms_setkey_string(outRec, "BUNIT_001", "1e20 dyne");
	drms_setkey_string(outRec, "BUNIT_002", "1e20 dyne");
	
	char timebuf[1024];
	float UNIX_epoch = -220924792.000; /* 1970.01.01_00:00:00_UTC */
	double val;
	int status = DRMS_SUCCESS;
	
	sprint_time(timebuf, (double)time(NULL) + UNIX_epoch, "ISO", 0);
	drms_setkey_string(outRec, "DATE", timebuf);
	
    char *cvsinfo = strdup("$Id: lorentz.c,v 1.2 2014/03/21 00:21:05 xudong Exp $");
    drms_setkey_string(outRec, "LOR_VERS", cvsinfo);        // Mar 12
	
};
