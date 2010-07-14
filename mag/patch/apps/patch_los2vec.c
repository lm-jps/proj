/*
 * Module name:
 *		patch_los2vec.c
 *
 * Description:
 *		Read in patch byte mask from LoS patch and keywords
 *		Copy them into corresponding vector patch
 *
 * Input:
 *		LoS patch (su_xudong.test_patch_720s)
 *		Vector magnetogram (su_keiji.vmagf_2d_720s_3_fltprf)
 *
 * Output:
 *		Vector patch (su_xudong.test_patch_v)
 *
 * Written by:
 *		Xudong Sun (xudongs@stanford.edu)
 *
 * Version:
 *		v1.0		Apr 15 2010
 *
 * Notes:
 *		v1.0
 *		1. Adapted from test_patch_v
 *		2. There needs to be one vector mag record in the time slot specified
 *		by the LoS patch query. Otherwise the code stops.
 *
 * Example:
 *		patch_los2vec "plos=su_xudong.test_patch_720s[2010.03.29_12:00:00_TAI][2]" "vmag=su_keiji.vmagf_2d_720s_3_fltprf" "pvec=su_xudong.test_patch_vec"
 *
 */


#include <jsoc_main.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>

#define DIE(msg) {fflush(stdout); fprintf(stderr, "%s, status=%d\n", msg, status); return(status);}

#define DTOR	(M_PI/180.)

// Macros for WCS transformations.  assume crpix1, crpix2 = CRPIX1, CRPIX2, sina,cosa = sin and cos of CROTA2 resp.
// and crvalx and crvaly are CRVAL1 and CRVAL2, cdelt = CDELT1 == CDELT2, then
// PIX_X and PIX_Y are CCD pixel addresses, WX and WY are arc-sec W and N on the Sun from disk center.
#define PIX_X(wx,wy) ((((wx-crvalx)*cosa + (wy-crvaly)*sina)/cdelt)+crpix1)
#define PIX_Y(wx,wy) ((((wy-crvaly)*cosa - (wx-crvalx)*sina)/cdelt)+crpix2)
#define WX(pix_x,pix_y) (((pix_x-crpix1)*cosa - (pix_y-crpix2)*sina)*cdelt+crvalx)
#define WY(pix_x,pix_y) (((pix_y-crpix2)*cosa + (pix_x-crpix1)*sina)*cdelt+crvaly)


/* ################## Timing from T. Larson ################## */

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
  *utime = ru.ru_utime.tv_sec * 1000.0 + ru.ru_utime.tv_usec / 1000.0;
  *stime = ru.ru_stime.tv_sec * 1000.0 + ru.ru_stime.tv_usec / 1000.0;
  return *utime + *stime;
}



/* ################## Main Module ################## */


char *module_name = "patch_los2vec";	/* Module name */

//char *vmag = "su_keiji.vmagf_2d_720s";	/* Vector magnetogram series */

ModuleArgs_t module_args[] =
{
    {ARG_STRING, "plos", NULL, "Input LoS patches."},
    {ARG_STRING, "vmag", NULL, "Input vector magnetograms."},
    {ARG_STRING, "pvec", NULL, "Output vector patches."},
    {ARG_INT, "VERB", "1", "Level of verbosity: 0=errors/warnings; 1=all messages"},
    {ARG_END}
};


int DoIt(void)
{
    int status = DRMS_SUCCESS;
    char *plosQuery, *vmag, *vmagQuery, *pvecQuery;
    DRMS_RecordSet_t *plosRS, *vmagRS, *pvecRS;
    int irec, nrecs, imag, nmag;
    DRMS_Record_t *plosRec, *vmagRec, *pvecRec;
    DRMS_Segment_t *plosSeg, *pvecSeg;
    DRMS_Array_t *plosArray, *pvecArray;
    DRMS_Link_t *srcLink;
    char *plosBitmap, *pvecBitmap;

    int i, verbflag, outDims[2];
    int n0, n1, mapsz;
    int pnum;
    TIME t_rec, dt;
    char *trec_str;

    // Time measuring
    double wt0, wt1, wt;
    double ut0, ut1, ut;
    double st0, st1, st;
    double ct0, ct1, ct;
    wt0 = getwalltime();
    ct0 = getcputime(&ut0, &st0);

    /* Get parameters */
    plosQuery = (char *)params_get_str(&cmdparams, "plos");
    vmag = (char *)params_get_str(&cmdparams, "vmag");
    pvecQuery = (char *)params_get_str(&cmdparams, "pvec");
    verbflag = params_get_int(&cmdparams, "VERB");

    /* Open LoS patch input */
    plosRS = drms_open_records(drms_env, plosQuery, &status);
    if (status || plosRS->n == 0) DIE("No magnetic input data found");
    nrecs = plosRS->n;

    /* Do this for each record */
    for (irec = 0; irec < nrecs; irec++)
    {
    
        if (verbflag) {
            wt1 = getwalltime();
            ct1 = getcputime(&ut1, &st1);
            printf("processing record %d...\n", irec);
        }
        
        /* Input record and data */
        plosRec = plosRS->records[irec];
        t_rec = drms_getkey_time(plosRec, "T_REC", &status);
        
        /* Search for available vector mag record */
        vmagQuery = (char *)malloc(100 * sizeof(char));
        trec_str = (char *)malloc(30 * sizeof(char));
        sprint_time(trec_str, t_rec, "TAI", 0);
        sprintf(vmagQuery, "%s[%s]", vmag, trec_str); printf("%s, %s\n", trec_str, vmagQuery);
        vmagRS = drms_open_records(drms_env, vmagQuery, &status);
        if (status || vmagRS->n != 1) DIE("No magnetic input data found");
        vmagRec = vmagRS->records[0];
        printf("#%d: t_vmag = %lf; t_plos = %lf\n", irec, drms_getkey_time(vmagRec, "T_REC", &status), t_rec);
        
        /* Data */
        plosSeg = drms_segment_lookup(plosRec, "bitmap");
        plosArray = drms_segment_read(plosSeg, DRMS_TYPE_CHAR, &status);
        if (status) {
            printf("No bitmap found for patch #%d. \n", irec);
            fflush(stdout);
            drms_free_array(plosArray);
            continue;
        }
        plosBitmap = (char *)plosArray->data;
        n0 = plosArray->axis[0]; n1 = plosArray->axis[1];
        mapsz = n0 * n1;
        
        /* Output array */
        outDims[0] = n0; outDims[1] = n1;
	    pvecArray = drms_array_create(DRMS_TYPE_CHAR, 2, outDims, NULL, &status);
        pvecBitmap = (char *)pvecArray->data;
        for (i = 0; i < mapsz; i++) {
            pvecBitmap[i] = plosBitmap[i];
        }
        
        /* Output record */
        pvecRec = drms_create_record(drms_env, pvecQuery, DRMS_PERMANENT, &status);
        if (status)
            DIE("Output record not created");
        pvecSeg = drms_segment_lookup(pvecRec, "bitmap");
        pvecSeg->axis[0] = pvecArray->axis[0];
        pvecSeg->axis[1] = pvecArray->axis[1];
        pvecArray->parent_segment = pvecSeg;
        
        /* Links */
        srcLink = hcon_lookup_lower(&pvecRec->links, "MDATA");
        if (srcLink) {
        	drms_setlink_static(pvecRec, "MDATA", vmagRec->recnum);
        }
        
        
        /*********************************/
        /* Ephemeris derivation, by Phil */
        /*********************************/
        // For CRVALn
        float crvalx = drms_getkey_float(vmagRec, "CRVAL1", &status);	// center of solar disc in arcsec
        float crvaly = drms_getkey_float(vmagRec, "CRVAL2", &status);
        float cdelt = drms_getkey_float(vmagRec, "CDELT1", &status);	// in arcsec, assumimg dx=dy
        float crpix1 = drms_getkey_float(vmagRec, "CRPIX1", &status);	// disk center in ccd, original
        float crpix2 = drms_getkey_float(vmagRec, "CRPIX2", &status);
        float crota2 = drms_getkey_float(vmagRec, "CROTA2", &status);	// rotation
        float sina = sin(crota2 * DTOR); 
        float cosa = cos(crota2 * DTOR);
        float pcx = drms_getkey_float(plosRec, "CRPIX1", &status);
        float pcy = drms_getkey_float(plosRec, "CRPIX2", &status);
        
        /* Keywords */
        drms_setkey_string(pvecRec, "BLD_VERS", jsoc_version);
        drms_setkey_time(pvecRec, "DATE", CURRENT_SYSTEM_TIME);
        // Prime key
        drms_copykey(pvecRec, vmagRec, "T_REC");
        drms_copykey(pvecRec, plosRec, "PNUM");
        // Geometry
        drms_copykey(pvecRec, plosRec, "HWIDTH1");
        drms_copykey(pvecRec, plosRec, "HWIDTH2");
        drms_copykey(pvecRec, plosRec, "CRPIX1");
        drms_copykey(pvecRec, plosRec, "CRPIX2");
        drms_copykey(pvecRec, vmagRec, "CDELT1");
        drms_copykey(pvecRec, vmagRec, "CDELT2");
        drms_setkey_float(pvecRec, "CRVAL1", WX(pcx,pcy));
        drms_setkey_float(pvecRec, "CRVAL2", WY(pcx,pcy));
        // Flags
        drms_copykey(pvecRec, plosRec, "PATCHNUM");
        drms_copykey(pvecRec, plosRec, "MASK");
        // Stats
        drms_copykey(pvecRec, plosRec, "TOTVALS");
        drms_copykey(pvecRec, plosRec, "DATAVALS");
        drms_copykey(pvecRec, plosRec, "MISSVALS");
        
        /* Finishing up */
        status = drms_segment_write(pvecSeg, pvecArray, 0);
        if (status) DIE("Problem writing file"); 
        drms_close_record(pvecRec, DRMS_INSERT_RECORD);

        /* Clean up */
        drms_free_array(plosArray);
        drms_free_array(pvecArray);
        free(trec_str); trec_str = NULL;
        free(vmagQuery); vmagQuery = NULL;
        drms_close_records(vmagRS, DRMS_FREE_RECORD);
        
        /* Time measure */
        if (verbflag) {
            wt = getwalltime();
            ct = getcputime(&ut, &st);
            printf("record %d done, %.2f ms wall time, %.2f ms cpu time\n", 
                     irec, wt - wt1, ct - ct1);
        }
        
    }
    
    drms_close_records(plosRS, DRMS_FREE_RECORD);

    if (verbflag) {
        wt = getwalltime();
        ct = getcputime(&ut, &st);
        printf("total time spent: %.2f ms wall time, %.2f ms cpu time\n", 
                wt - wt0, ct - ct0);
    }

    return(DRMS_SUCCESS);

}  // End of module
