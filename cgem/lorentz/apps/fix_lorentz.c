/*
 * Name:    fix_lorentz
 *
 * Purpose:
 *          Populate a few keywords from hmi.sharp_cea_720s to cgem.Lorentz
 *
 * Written by:
 *          Xudong Sun (xudong@sun.stanford.edu)
 *
 * Usage:
 *          fix_lorentz "in=cgem.Lorentz[377][2011.02.15_01/2h]"
 *
 * Version:
 *          v1.0    2014 May 16
 *
 * Notes:
 *          Take advantage of the fact that hmi.sharp_cea_720s is dynamically linked
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>
#include "jsoc_main.h"

#define DIE(msg) {fflush(stdout); fprintf(stderr,"%s, status=%d\n", msg, status); return(status);}
#define SHOW(msg) {printf("%s", msg); fflush(stdout);}
//#define UNIX_EPOCH (-220924792.000) // 1970.01.01_00:00:00_UTC

char *module_name = "fix_lorentz";

ModuleArgs_t module_args[] =
{
	{ARG_STRING, "in", NULL, "Input/output lorentz series"},
    {ARG_END}
};

int DoIt(void)
{
    
    char *cvsinfo = strdup("$Id: fix_lorentz.c,v 1.1 2014/05/17 01:28:49 xudong Exp $");

    int status = DRMS_SUCCESS;
    
    char *inQuery = (char *) params_get_str(&cmdparams, "in");

    // Input
    
    DRMS_RecordSet_t *inRS = drms_open_records(drms_env, inQuery, &status);
    int nrecs = inRS->n;
    // record set name???
    if (status || (nrecs == 0)) DIE("Error open data series");
    
    // Clone as output
    
    DRMS_RecordSet_t *outRS = drms_clone_records(inRS, DRMS_PERMANENT, DRMS_SHARE_SEGMENTS, &status);    // share segments
    if (status) DIE("Error cloning data series");
    
    // No longer needed
    
    drms_close_records(inRS, DRMS_FREE_RECORD);
    
    // Count for success
    
    int iGood = 0;
    
    for (int irec = 0; irec < nrecs; irec++) {
        
        // Record
        
        DRMS_Record_t *outRec = outRS->records[irec];
        
        TIME t_rec = drms_getkey_time(outRec, "T_REC", &status);
        char t_rec_str[100];
        sprint_time(t_rec_str, t_rec, "TAI", 0);
        
        printf("Record #%d of %d, T_REC=%s...\n", irec + 1, nrecs, t_rec_str);
        
        // Get linked hmi.sharp_cea_720s record

        DRMS_Record_t *sharpRec = drms_link_follow(outRec, "SHARP", &status);
        if (status) {
            SHOW("Error retrieving sharp record, skipped\n");
            continue;
        }
        
        // Copy keys
        
        drms_copykey(outRec, sharpRec, "NOAA_AR");
        drms_copykey(outRec, sharpRec, "NOAA_NUM");
        drms_copykey(outRec, sharpRec, "NOAA_ARS");
        drms_copykey(outRec, sharpRec, "AMBPATCH");
        
        char timebuf[512], histbuf[2014];
        double val;
        sprint_time(timebuf, (double)time(NULL) + UNIX_EPOCH, "ISO", 0);
        drms_setkey_string(outRec, "DATE", timebuf);
        
        sprintf(histbuf, "NOAA keywords fixed with %s", cvsinfo);
        drms_setkey_string(outRec, "HISTORY", histbuf);
        
        //
        
        drms_close_record(sharpRec, DRMS_FREE_RECORD);
        iGood++;
        
    }
    
    // Close RS
    
    drms_close_records(outRS, DRMS_INSERT_RECORD);
    
    printf("%d of %d records fixed.\n", iGood, nrecs);
    
    return DRMS_SUCCESS;
    
}