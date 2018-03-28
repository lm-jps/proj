/*
 * Module name:
 *		m2lostarp.c
 *
 * Description:
 *		Convert Mtarp into lostarp series
 *		Copy bitmap, create links to los segments and copy keywords
 *		Adapted from patch_los2vec.c by X. Sun
 *
 * Input:
 *		MDI LoS patch (mdi.Mtarp)
 *		MDI Los magnetogram (mdi.fd_M_96m_lev182)
 *
 * Output:
 *		LoS patch (mdi.lostarp_96m)
 *
 * Written by:
 *		Monica Bobra and Xudong Sun
 *
 * Version:
 *              v0.0 9 February 2018 Monica Bobra
 * Notes:
 *		v0.0 No explicit notes
 *
 * Example:
 *              > m2lostarp "mtarp=mdi.Mtarp[][2010.10.14_20:48:00_TAI]" "los=mdi.fd_M_96m_lev182" "lostarp=mdi.lostarp_96m"
 */


#include <jsoc_main.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "copy_me_keys.c"
#include "fstats.h"

#define DIE(msg) {fflush(stdout); fprintf(stderr, "%s, status=%d\n", msg, status); return(status);}
#define SHOW(msg) {printf("%s", msg); fflush(stdout);}

#define DTOR	(M_PI/180.)

#ifndef MIN
#define MIN(a,b) (((a)<(b)) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) (((a)>(b)) ? (a) : (b))
#endif

// Macros for WCS transformations.  assume crpix1, crpix2 = CRPIX1, CRPIX2, sina,cosa = sin and cos of CROTA2 resp.
// and crvalx and crvaly are CRVAL1 and CRVAL2, cdelt = CDELT1 == CDELT2, then
// PIX_X and PIX_Y are CCD pixel addresses, WX and WY are arc-sec W and N on the Sun from disk center.
#define PIX_X(wx,wy) ((((wx-crvalx)*cosa + (wy-crvaly)*sina)/cdelt)+crpix1)
#define PIX_Y(wx,wy) ((((wy-crvaly)*cosa - (wx-crvalx)*sina)/cdelt)+crpix2)
#define WX(pix_x,pix_y) (((pix_x-crpix1)*cosa - (pix_y-crpix2)*sina)*cdelt+crvalx)
#define WY(pix_x,pix_y) (((pix_y-crpix2)*cosa + (pix_x-crpix1)*sina)*cdelt+crvaly)


/* ################## Main Module ################## */


char *module_name = "m2lostarp";	/* Module name */

ModuleArgs_t module_args[] =
{
    {ARG_STRING, "mtarp", NULL, "Input LoS HARP series."},
    {ARG_STRING, "los", NULL, "Input ME inversion series, series name only."},
    {ARG_STRING, "lostarp", NULL, "Output ME HARP series."},
    //{ARG_STRING, "mharp", NULL, "Input LoS HARP series."}, mharp->mtarp
    //{ARG_STRING, "me", NULL, "Input ME inversion series, series name only."}, me->los
    //{ARG_STRING, "meharp", NULL, "Output ME HARP series."}, meharp->lostarp
    {ARG_INT, "buffx", "0", "Buffer around the identified patch, as needed."},
    {ARG_INT, "buffy", "0", "Buffer around the identified patch, as needed."},
    {ARG_END}
};

int DoIt(void)
{
    int status = DRMS_SUCCESS;
	
	char *mHarpQuery, *meQuery, *meHarpQuery, *meSeriesName;
	
	DRMS_RecordSet_t *mHarpRS=NULL, *meRS=NULL, *meHarpRS=NULL;
	DRMS_Record_t *mHarpRec=NULL, *meRec=NULL, *meHarpRec=NULL;
	DRMS_Segment_t *mHarpSeg=NULL, *meHarpSeg=NULL;
	DRMS_Array_t *mHarpArray=NULL, *meHarpArray=NULL;
	DRMS_Link_t *meLink=NULL;
	DRMS_Link_t *pLink=NULL, *pRec=NULL;

	char *mBitmap=NULL, *meBitmap=NULL;
	
	int irec, nrecs, imag, nmag;
	int nx, ny, mapsz;		// input bitmap size
	int outDims[2], buffx, buffy;	// output size, buffer size
	int val;	
	TIME t_rec;
	char *t_rec_str=NULL;
	
	/* Input */
	
	mHarpQuery = (char *) params_get_str(&cmdparams, "mtarp");
	meSeriesName = (char *) params_get_str(&cmdparams, "los");
	meHarpQuery = (char *) params_get_str(&cmdparams, "lostarp");
	
        printf("mHarpQuery %s\n", mHarpQuery);
        printf("meSeriesName %s\n", meSeriesName);
        printf("meHarpQuery %s\n", meHarpQuery);

	buffx = params_get_int(&cmdparams, "buffx"); buffx = MAX(buffx,0);
	buffy = params_get_int(&cmdparams, "buffy"); buffy = MAX(buffy,0);
	
	/* Access Mharp input */
	
	mHarpRS = drms_open_records(drms_env, mHarpQuery, &status);
        if (status || mHarpRS->n == 0) DIE("No Mtarp input data found");
        nrecs = mHarpRS->n;
	
	/* Loop through each record */
	
	for (irec = 0; irec < nrecs; irec++) {
		
		/* Mharp record */
		
                mHarpRec = mHarpRS->records[irec];
                t_rec = drms_getkey_time(mHarpRec, "T_REC", &status);
		printf("Record #%d of #%d\n", irec+1, nrecs); fflush(stdout);
		
		/* ME record */ 
		
		meQuery = (char *) malloc(300 * sizeof(char));
		t_rec_str = (char *) malloc(100 * sizeof(char));
		
		sprint_time(t_rec_str, t_rec, "TAI", 0);
		sprintf(meQuery, "%s[%s]\n", meSeriesName, t_rec_str); // printf("%s\n", trec_str);

                printf(meQuery, "%s\n", meQuery); fflush(stdout);
		
		meRS = drms_open_records(drms_env, meQuery, &status);

		if (status || meRS->n < 1) {
		        //SHOW("No Los input data found, skip current record\n");
			free(t_rec_str); free(meQuery);
			if (meRS) drms_close_records(meRS, DRMS_FREE_RECORD);
			t_rec_str = NULL; meQuery = NULL;
			meRS = NULL;
			DIE("No MDI LOS input data found, skip current record\n");		// changed on Dec 17 2012, exit when no ME is available
			continue;
		}
		
		for (int ii = 0; ii < meRS->n; ii++) {
			meRec = meRS->records[ii];
			
			/* Data */
			
			mHarpSeg = drms_segment_lookup(mHarpRec, "bitmap");
			mHarpArray = drms_segment_read(mHarpSeg, DRMS_TYPE_CHAR, &status);
			if (status) {
				SHOW("No bitmap found, skip current record\n");
				drms_free_array(mHarpArray);
				continue;
			}
			mBitmap = (char *) mHarpArray->data;
			nx = mHarpArray->axis[0]; ny = mHarpArray->axis[1];
			mapsz = nx * ny;
			
			/* Output array */
			
			outDims[0] = nx + 2 * buffx; outDims[1] = ny + 2 * buffy;
			meHarpArray = drms_array_create(DRMS_TYPE_CHAR, 2, outDims, NULL, &status);
			meBitmap = (char *) meHarpArray->data;
			
			for (int j = 0; j < ny; j++) {
				for (int i = 0; i < nx; i++) {
					meBitmap[(j + buffy) * outDims[0] + (i + buffx)] = mBitmap[j * nx + i];
				}
			}

			/* Output record */
			
			meHarpRec = drms_create_record(drms_env, meHarpQuery, DRMS_PERMANENT, &status);
			if (status) {
				SHOW("Output error, skip current record\n");
				drms_free_array(mHarpArray); drms_free_array(meHarpArray);
				continue;
			}
			meHarpSeg = drms_segment_lookup(meHarpRec, "bitmap");
			meHarpSeg->axis[0] = outDims[0];
			meHarpSeg->axis[1] = outDims[1];
			meHarpArray->parent_segment = meHarpSeg;
			
			/* Links */

			meLink = hcon_lookup_lower(&meHarpRec->links, "MDATA");
			if (meLink) {
				drms_setlink_static(meHarpRec, "MDATA", meRec->recnum);
			}

			pLink = hcon_lookup_lower(&meHarpRec->links, "PDATA");
			if (pLink) {
				drms_setlink_static(meHarpRec, "PDATA", meRec->recnum);
			}
			       
			/* Keywords */
                        // not setting any statistics keywords becuase it doesn't make sense for a bitmap
			drms_copykey(meHarpRec, meRec, "T_REC");
	                val = drms_getkey_double(mHarpRec, "HARPNUM", &status);
                        drms_setkey_double(meHarpRec, "TARPNUM", val);	
			drms_copykey(meHarpRec, mHarpRec, "HARPNUM");
			drms_copykey(meHarpRec, meRec, "RUNNUM");
			copy_me_keys(meRec, meHarpRec);
			copy_patch_keys(mHarpRec, meHarpRec);
			copy_geo_keys(mHarpRec, meHarpRec);
			drms_copykey(meHarpRec, mHarpRec, "NPIX_SPT");
			drms_copykey(meHarpRec, mHarpRec, "T_MAXPIX");
			drms_copykey(meHarpRec, mHarpRec, "ARS_NCLN");
			drms_copykey(meHarpRec, mHarpRec, "ARS_MODL");
			drms_copykey(meHarpRec, mHarpRec, "ARS_EDGE");
			drms_copykey(meHarpRec, mHarpRec, "ARS_BETA");
			drms_copykey(meHarpRec, mHarpRec, "T_MID1");
			drms_copykey(meHarpRec, mHarpRec, "T_CMPASS");
			drms_copykey(meHarpRec, mHarpRec, "SIZE_SPT");
			drms_copykey(meHarpRec, mHarpRec, "AREA_SPT");

                        // Set time
			drms_setkey_string(meHarpRec, "BLD_VERS", jsoc_version);
			char timebuf[1024];
			float UNIX_epoch = -220924792.000; /* 1970.01.01_00:00:00_UTC */
			sprint_time(timebuf, (double)time(NULL) + UNIX_epoch, "ISO", 0);
			drms_setkey_string(meHarpRec, "DATE", timebuf);
			// Misc
			drms_setkey_string(meHarpRec, "BUNIT_000", " ");
			drms_setkey_float(meHarpRec, "CRPIX1", 
							  drms_getkey_float(mHarpRec, "CRPIX1", &status)-buffx);
			drms_setkey_float(meHarpRec, "CRPIX2", 
							  drms_getkey_float(mHarpRec, "CRPIX2", &status)-buffy);

                        // set cvs commit version into keyword CODEVER6
                        char *cvsinfo  = strdup("$Id");
                        status = drms_setkey_string(meHarpRec, "CODEVER6", cvsinfo);

			/* Write data */
			status = drms_segment_write(meHarpSeg, meHarpArray, 0);
			if (status) {
				SHOW("Output error, skip current record\n");
				drms_free_array(mHarpArray); drms_free_array(meHarpArray);
				continue;
			} 
			
			drms_free_array(meHarpArray);
			drms_close_record(meHarpRec, DRMS_INSERT_RECORD);
			
		} //ii

		/* Clean up */
        drms_free_array(mHarpArray); meHarpArray = NULL;
        free(t_rec_str); free(meQuery); t_rec_str = NULL; meQuery = NULL;
        drms_close_records(meRS, DRMS_FREE_RECORD); meRS = NULL;
		
	} //irec
	
	drms_close_records(mHarpRS, DRMS_FREE_RECORD);
	
	return 0;
	
}
