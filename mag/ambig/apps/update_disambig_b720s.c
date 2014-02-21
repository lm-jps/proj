/*
 * Name:    update_disambig_b720s
 *
 * Purpose:
 *          Fix the second/third bit in the disambig segment of the full disk disambiguation series,
 *          by default hmi.B_720s, processed with CODEVER5 before v1.17. It looks at conf_disambig, look for
 *          pixels of 60, and make all three bits agreeing with the lowest bit. All other info are cloned over,
 *          using the DRMS_COPY_SEGMENTS mode. HISTORY and DATE keywords are changed.
 *
 * Written by:
 *          Xudong Sun (xudong@sun.stanford.edu)
 *
 * Usage:
 *          update_disambig_b720s b="hmi_test.B_720s[2010.05.01_18:12]"
 *
 * Version:
 *          v1.0    2014 Feb 20
 *
 * Notes:
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

#define DIE(msg) {fflush(stdout); fprintf(stderr,"%s, status=%d\n", msg, status); return(status);}
#define SHOW(msg) {printf("%s", msg); fflush(stdout);}
//#define UNIX_EPOCH (-220924792.000) // 1970.01.01_00:00:00_UTC

char *module_name = "update_disambig_b720s";
char *version_id  = "2014 Feb 14";

ModuleArgs_t module_args[] =
{
	{ARG_STRING, "b",     NULL, "Input/output B series"},
    {ARG_END}
};

int DoIt(void)
{
    
    char *cvsinfo = strdup("$Id: update_disambig_b720s.c,v 1.1 2014/02/21 20:55:02 xudong Exp $");
    
    int status = DRMS_SUCCESS;
    
    char *bQuery = (char *) params_get_str(&cmdparams, "b");
    
    // Input
    
    DRMS_RecordSet_t *inRS = drms_open_records(drms_env, bQuery, &status);
    int nrecs = inRS->n;
    if (status || (nrecs == 0)) DIE("Error open data series");
    
    // Clone as output
    
    DRMS_RecordSet_t *outRS = drms_clone_records(inRS, DRMS_PERMANENT, DRMS_COPY_SEGMENTS, &status);    // copy segments
    if (status) DIE("Error cloning data series");
    
    // No longer needed
    
    drms_close_records(inRS, DRMS_FREE_RECORD);
    
    // Loop
    
    for (int irec = 0; irec < nrecs; irec++) {
        
        // Record
        
//        DRMS_Record_t *inRec = inRS->records[irec];
        DRMS_Record_t *outRec = outRS->records[irec];
        
        TIME t_rec = drms_getkey_time(outRec, "T_REC", &status);
        char t_rec_str[100];
		    sprint_time(t_rec_str, t_rec, "TAI", 0);
        
        printf("Record #%d of %d, T_REC=%s...\n", irec + 1, nrecs, t_rec_str);
        
        // Get conf_ambig and disambig
        
        DRMS_Segment_t *confSeg = drms_segment_lookup(outRec, "conf_disambig");
        DRMS_Array_t *confArray = drms_segment_read(confSeg, DRMS_TYPE_CHAR, &status);
        if (status) {
            printf("Conf_disambig read error, skipping record #%d\n", irec + 1);
            continue;
        }
        char *conf = (char *) confArray->data;
        
        DRMS_Segment_t *disambSeg = drms_segment_lookup(outRec, "disambig");
        DRMS_Array_t *disambArray = drms_segment_read(disambSeg, DRMS_TYPE_CHAR, &status);
        if (status) {
            printf("Disambig read error, skipping record #%d\n", irec + 1);
            drms_free_array(confArray);
            continue;
        }
        char *disamb = (char *) disambArray->data;
        
        // Create a new disambig array
        
        int nx = disambArray->axis[0], ny = disambArray->axis[1], nxny = nx * ny;
        char *disamb_new = (char *) (malloc(nxny * sizeof(char)));
        
        for (int i = 0; i < nxny; i++) {
            disamb_new[i] = disamb[i];
            if (conf[i] > 85 || conf[i] < 55) continue;     // for 60 only
            disamb_new[i] = (disamb[i] % 2) ? 7 : 0;            // odd=>111, even=>000
        }
        
        // Clean up #1
        
        drms_free_array(confArray);
        drms_free_array(disambArray);
        
        // Output new disambig
        
        int outDims[2] = {nx, ny};
        
        DRMS_Segment_t *outSeg = disambSeg;         // Write to the same segment. Can we?
        DRMS_Array_t *outArray = drms_array_create(DRMS_TYPE_CHAR, 2, outDims, disamb_new, &status);
        if (status) {
            printf("Error while creating disambig, skip record #%d\n", irec + 1);
            continue;
        }
        outArray->israw = 0;
        outArray->bzero = outSeg->bzero;
        outArray->bscale = outSeg->bscale;
        
        status = drms_segment_write(outSeg, outArray, 0);
        if (status) {
            printf("Error while writing disambig, skip record #%d\n", irec + 1);
            drms_free_array(outArray);
            continue;
        }
        
        // Keywords
        
        char timebuf[512], histbuf[2014];
        double val;
        sprint_time(timebuf, (double)time(NULL) + UNIX_EPOCH, "ISO", 0);
        drms_setkey_string(outRec, "DATE", timebuf);
        
        sprintf(histbuf, "disambig fixed with %s", cvsinfo);
        drms_setkey_string(outRec, "HISTORY", histbuf);
        
        // Clean up #2
        
        drms_free_array(outArray);
        
    }
    
    // Close RS
    
    drms_close_records(outRS, DRMS_INSERT_RECORD);
    
    return DRMS_SUCCESS;
    
}
