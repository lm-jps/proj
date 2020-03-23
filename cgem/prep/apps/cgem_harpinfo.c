/*
 *  cgem_harpinfo.c
 *
 *  This module takes a SHARP number and returns the following information:
 *  + Overall success (1) or failure (0). Can encode more info later
 *  + HARP number, NOAA number (0 if unavailable)
 *  + T_OBS of first & last frame where patch center is within 70 longitude
 *  + T_OBS of reference frame, where patch center is closest to CM
 *  + Carrington lon/lat of patch center of reference frame; rounded to integer times 0.03
 *  + Minimal frame size in native coordinate that encompass all frames: max(CRSIZEn) with 50 pixel padding on each side
 *  + Recommended map size in Plate Caree: round(LONDTMAX-LONDTMIN)/0.03
 *
 *  Author:
 *      Xudong Sun
 *
 *  Version:
 *      v0.0 Dec 27 2019
 *      v0.1 Mar 23 2020
 *
 *  Notes:
 *      v0.1
 *      Declare success even no unique NOAA AR is found, as we decide
 *      CGEMNUM=HARPNUM+100000
 *
 *  Example Calls:
 *      > cgem_harpinfo "harpnum=377" -h -i    # h for readable format
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>
#include "jsoc_main.h"

#define PI              (M_PI)
#define RADSINDEG       (PI/180.)
#define RAD2ARCSEC      (648000./M_PI)
#define SECINDAY        (86400.)
#define FOURK           (4096)
#define FOURK2          (16777216)
#define NMAXSTR         (256)       // max string length
#define LON_THRESH      (70.)       // longitude threshold
#define NPAD            (50)        // padding columns on each size
#define MPAD            (50)        // padding rows on each size
#define DPC             (0.03)      // Plate-Carree map resolution

#define ARRLENGTH(ARR) (sizeof(ARR) / sizeof(ARR[0]))

// Some other things
#ifndef MIN
#define MIN(a,b) (((a)<(b)) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) (((a)>(b)) ? (a) : (b))
#endif

#define DIE(msg) {fflush(stdout); fprintf(stderr,"%s, status=%d\n", msg, status); return(status);}
#define SHOW(msg) {printf("%s", msg); fflush(stdout);}

#define kNotSpecified "Not Specified"
#define dpath    "/home/jsoc/cvs/Development/JSOC"

// Necessary information
struct recInfo {
    TIME t_rec;
//    double disk_lonc;           // Carrington lon of disk center
    double lonc, latc;          // Stonyhurst lon/lat of patch center
    double dlon, dlat;          // Size of patch in Plate Caree
    int crsize1, crsize2;       // Cutout sizes in HARP
};

/* ========================================================================================================== */

char *module_name = "cgem_harpinfo";

ModuleArgs_t module_args[] =
{
    {ARG_INT, "harpnum", "-1", "HARP number."},
    {ARG_FLAG, "h", "", "Output with readable comments."},
    {ARG_FLAG, "i", "", "Output additional info."},
    {ARG_END}
};

int DoIt(void)
{
    
    int status = DRMS_SUCCESS;
    
    /* Get data series */
    
    int harpnum = params_get_int(&cmdparams, "harpnum");
    char mharpQuery[NMAXSTR];
    sprintf(mharpQuery, "hmi.Mharp_720s[%d]", harpnum);
    int help = params_isflagset(&cmdparams, "h");
    int info = params_isflagset(&cmdparams, "i");
    
    /* Output variables */
    
    int result = 0;                 // 1 for success, 0 for failure
    int noaa_ar = DRMS_MISSING_INT, noaa_num = 0;      // NOAA numbers
    char *noaa_ars;
    TIME t_start = DRMS_MISSING_TIME, t_end = DRMS_MISSING_TIME, t_ref = DRMS_MISSING_TIME;     // Times
    char tstr_start[NMAXSTR] = "", tstr_end[NMAXSTR] = "", tstr_ref[NMAXSTR] = "";
    double lonc_start = DRMS_MISSING_DOUBLE, latc_start = DRMS_MISSING_DOUBLE;
    double lonc_end = DRMS_MISSING_DOUBLE, latc_end = DRMS_MISSING_DOUBLE;
    double lonc_ref = DRMS_MISSING_DOUBLE, latc_ref = DRMS_MISSING_DOUBLE;          // Reference Carrington lon/lat
    int cols = 0, rows = 0;                 // Patch size for cut out
    int cols_pc = 0, rows_pc = 0;             // Patch size for Plate Carree map
    int nGoodRecs = 0;
    TIME t_frst1 = DRMS_MISSING_TIME, t_last1 = DRMS_MISSING_TIME;      // start/end time with strong field
    double dlon_ref, dlat_ref;          // Patch size in degree
    
    /* Input data */
    
    DRMS_RecordSet_t *mharpRS = drms_open_records(drms_env, mharpQuery, &status);        // Input data series
    int nrecs = mharpRS->n;
    if (status || nrecs == 0 || !mharpRS) { printf("%d\n", result); return(status);}
//    printf("%s has %d records.\n", mharpQuery, nrecs); fflush(stdout);
    
    /* Start */
    
    result = 1;     // assume success; check against all criteria below
    
    struct recInfo *recInfoArr = (struct recInfo *) (malloc(nrecs * sizeof(struct recInfo)));

    // Loop over record set to retrieve necessary info
    // No error checking
    
    for (int irec = 0; irec < nrecs; irec++) {
        
        DRMS_Record_t *mharpRec = mharpRS->records[irec];
        
        recInfoArr[irec].t_rec = drms_getkey_time(mharpRec, "T_REC", &status);
        
        double latdtmin = drms_getkey_double(mharpRec, "LATDTMIN", &status);
        double latdtmax = drms_getkey_double(mharpRec, "LATDTMAX", &status);
        double londtmin = drms_getkey_double(mharpRec, "LONDTMIN", &status);
        double londtmax = drms_getkey_double(mharpRec, "LONDTMAX", &status);
        recInfoArr[irec].latc = (latdtmin + latdtmax) / 2.;
        recInfoArr[irec].lonc = (londtmin + londtmax) / 2.;
        recInfoArr[irec].dlat = latdtmax - latdtmin;
        recInfoArr[irec].dlon = londtmax - londtmin;
        
        recInfoArr[irec].crsize1 = drms_getkey_int(mharpRec, "CRSIZE1", &status);
        recInfoArr[irec].crsize2 = drms_getkey_int(mharpRec, "CRSIZE2", &status);
        
    }
    
    noaa_ar = drms_getkey_int(mharpRS->records[0], "NOAA_AR", &status);
    noaa_num = drms_getkey_int(mharpRS->records[0], "NOAA_NUM", &status);
    noaa_ars = drms_getkey_string(mharpRS->records[0], "NOAA_ARS", &status);
    
    t_frst1 = drms_getkey_time(mharpRS->records[0], "T_FRST1", &status);
    t_last1 = drms_getkey_time(mharpRS->records[0], "T_LAST1", &status);
        
    // Process info
    
    int index_start = nrecs - 1, index_end = 0, index_ref = 0;
    
    for (int irec = 0; irec < nrecs; irec++) {
        
        if (fabs(recInfoArr[irec].lonc) > LON_THRESH ||
            recInfoArr[irec].t_rec < t_frst1 ||
            recInfoArr[irec].t_rec > t_last1) { continue; }
        
        nGoodRecs++;
        if (recInfoArr[irec].lonc < recInfoArr[index_start].lonc) { index_start = irec; }
        if (recInfoArr[irec].lonc > recInfoArr[index_end].lonc) { index_end = irec; }
        if (fabs(recInfoArr[irec].lonc) < fabs(recInfoArr[index_ref].lonc)) { index_ref = irec; }
        if (recInfoArr[irec].crsize1 > cols) { cols = recInfoArr[irec].crsize1; }
        if (recInfoArr[irec].crsize2 > rows) { rows = recInfoArr[irec].crsize2; }
        
    }
    
//    if (nGoodRecs == 0 || isnan(noaa_ar) || noaa_num == 0) { result = 0; }
    if (nGoodRecs == 0) { result = 0; }     // Mar 23 2020
    
    lonc_start = recInfoArr[index_start].lonc;
    latc_start = recInfoArr[index_start].latc;
    lonc_end = recInfoArr[index_end].lonc;
    latc_end = recInfoArr[index_end].latc;
    lonc_ref = round(recInfoArr[index_ref].lonc / DPC) * DPC;
    latc_ref = round(recInfoArr[index_ref].latc / DPC) * DPC;
    
    sprint_time(tstr_start, recInfoArr[index_start].t_rec, "TAI", 0);
    sprint_time(tstr_end, recInfoArr[index_end].t_rec, "TAI", 0);
    sprint_time(tstr_ref, recInfoArr[index_ref].t_rec, "TAI", 0);
    
    cols += (NPAD * 2); rows += (MPAD *2);
    
    dlat_ref = recInfoArr[index_ref].dlat;
    dlon_ref = recInfoArr[index_ref].dlon;
    rows_pc = round(dlat_ref / DPC);
    cols_pc = round(dlon_ref / DPC);
    
    
    /* Print */
    
    if (help) {     // print comments and arguments
        
        printf("success: %d\n", result);
        printf("harpnum: %d\n", harpnum);
        printf("noaa_ar: %d\n", noaa_ar);        // 0 if unavailable
        printf("t_start: %s\n", tstr_start);
        printf("t_end: %s\n", tstr_end);
        printf("t_ref: %s\n", tstr_ref);
        printf("lonref: %f\n", lonc_ref);
        printf("latref: %f\n", latc_ref);
        printf("cols: %d\n", cols);
        printf("rows: %d\n", rows);
        printf("cols_pc: %d\n", cols_pc);
        printf("rows_pc: %d\n", rows_pc);
        
    } else {        // print arguments only
        
        printf("%d\n", result);
        printf("%d\n", harpnum);
        printf("%d\n", noaa_ar);
        printf("%s\n", tstr_start);
        printf("%s\n", tstr_end);
        printf("%s\n", tstr_ref);
        printf("%f\n", lonc_ref);
        printf("%f\n", latc_ref);
        printf("%d\n", cols);
        printf("%d\n", rows);
        printf("%d\n", cols_pc);
        printf("%d\n", rows_pc);
        
    }
    
    if (info) {     // print additional info
        
        printf("nframes: %d\n", nGoodRecs);
        printf("noaa_num: %d\n", noaa_num);
        printf("noaa_ars: %s\n", noaa_ars);
        printf("lonc_start: %f\n", lonc_start);
        printf("lonc_end: %f\n", lonc_end);
        printf("latc_start: %f\n", latc_start);
        printf("latc_end: %f\n", latc_end);
        printf("dlon: %f\n", dlon_ref);
        printf("dlat: %f\n", dlat_ref);
        
    }
    
    /* Clean up */
    
    free(recInfoArr);
    drms_close_records(mharpRS, DRMS_FREE_RECORD);

    return 0;

}   // DoIt
