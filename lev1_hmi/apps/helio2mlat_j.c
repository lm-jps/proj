#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include "jsoc_main.h"
#include "astro.h"
#include "drms_dsdsapi.h"
#include "errlog.h"


#define PI		(M_PI)
#define kRecSetIn "in"
#define kSeriesOut "out"
#define kParamSegIn "segin"
#define kParamSegOut "segout"
#define kParamLMAX "LMAX"
#define kParamSUBMEAN "S"
#define kParamNORMALIZE "N"
#define kParamCENTLONG "C"
#define kParamZEROMISS "Z"
#define kParamLGAPOD "L"
#define kParamLGAPMIN "LGAPMIN"
#define kParamLGAPWIDT "LGAPWIDT"

#define kDefSegIn "data"
#define kDefSegOut "data"

#define kI_DREC         "I_DREC"
#define kT_REC          "T_REC"
#define kQUALITY        "QUALITY"
#define kQUAL_ALL_DATA_MISSING       (0)
#define kQUAL_DATA_AVAILABLE         (1)

typedef enum
{
   H2MStatus_Success,
   H2MStatus_MissingParameter,
} H2MStatus_t;

char *module_name = "helio2mlat";

ModuleArgs_t module_args[] = 
{
   /* Should be ARG_DATASET, but not implemented */
   {ARG_STRING, kRecSetIn,      NULL,       "Input data records",   NULL},
   /* Should be ARG_DATASERIES, but not implemented */
   {ARG_STRING, kSeriesOut,     NULL,       "Output data series",   NULL},
   {ARG_STRING, kParamSegIn,    kDefSegIn,  "Input data segment",   NULL},
   {ARG_STRING, kParamSegOut,   kDefSegOut, "Output data segment",  NULL},
   {ARG_INT,    kParamLMAX,    "249",       NULL,                   NULL},
   /* These flags used to be <name>=<value> parameters in SOI.
    * So, rather than put SUBMEAN=1 on the cmd-line, you'd 
    * put "-S" on the cmd-line. */
   {ARG_FLAG,    kParamSUBMEAN,  "0",       NULL,                   NULL},
   {ARG_FLAG,    kParamNORMALIZE, "0",      NULL,                   NULL},
   {ARG_FLAG,    kParamCENTLONG, "1",       NULL,                   NULL},
   {ARG_FLAG,    kParamZEROMISS, "0",       NULL,                   NULL},
   {ARG_FLAG,    kParamLGAPOD,   "0",       NULL,                   NULL},  /* Apodize in longitude? */
   /* Start of apodize. Degrees. */
   {ARG_DOUBLE,  kParamLGAPMIN,  "60.0",    NULL,                   "[0.0,90.0]"}, /* range unimplemented */
   /* Width of apodize. Degrees. */
   {ARG_DOUBLE,  kParamLGAPWIDT, "10.0",    NULL,                   "[0.0,90.0]"}, /* range unimplemented */
   {ARG_END,     NULL,           NULL,      NULL,                   NULL}
};

/* Segment will be empty, but there will be a record! */
static void CreateBlankRecord(DRMS_Record_t *inrec, DRMS_Record_t *outrec)
{
  /* create 'blank' data */
   drms_copykey(outrec, inrec, kI_DREC);
   drms_copykey(outrec, inrec, kT_REC);
   if (drms_copykey(outrec, inrec, kQUALITY))
   {
      drms_setkey_int(outrec, kQUALITY, kQUAL_ALL_DATA_MISSING);
   }
}

int DoIt(void)
{
   int status = DRMS_SUCCESS;
   ErrNo_t error;
   WarnNo_t warn;

   char *inRecQuery = NULL;
   char *outseries = NULL;
   char *inseries = NULL;
   DRMS_RecordSet_t *inrecs = NULL;

   double mean, norm=1.0, normx;
   double tmp;
   int subtract_mean, normalize, cent_long, zero_miss, lgapod;
   double lgapmin, lgapwidt, lgapmax, lon;

   int col;
   int requested_lmax;
   int lmax, len, lfft, mapped_lmax, map_cols, map_rows; 
   int map_cols2;
   int nfft, nmean, nok, nout, row, sn, verbose;
   int forward = 1; /* forward transform option */
   int length[2];
   int qual;

   clock_t t0;

   LogPrint_t printerr = errlog_errwrite;
   LogPrint_t printhst = errlog_hstwrite;

   requested_lmax = params_get_int(&cmdparams, kParamLMAX);
   subtract_mean = params_isflagset(&cmdparams, kParamSUBMEAN);
   normalize = params_isflagset(&cmdparams, kParamNORMALIZE);
   /* CENTLONG=1 centers the longitude Fourier transform on the center
      of the remapped image */
   cent_long = params_isflagset(&cmdparams, kParamCENTLONG);
   /* ZEROMISS=1 sets missing data to 0,
      ZEROMISS=0 fills the output row with missing */
   zero_miss = params_isflagset(&cmdparams, kParamZEROMISS);
   lgapod = params_isflagset(&cmdparams, kParamZEROMISS);
   lgapmin = params_get_double(&cmdparams, kParamLGAPMIN);
   lgapwidt = params_get_double(&cmdparams, kParamLGAPWIDT);
   lgapmax = lgapmin+lgapwidt;

   t0 = clock();
   (*printhst)("Setup %f\n", t0);

   inRecQuery = cmdparams_get_str(&cmdparams, kRecSetIn, NULL);
   outseries = cmdparams_get_str(&cmdparams, kSeriesOut, NULL);
   inseries = drms_recordset_acquireseriesname(inRecQuery);

   inrecs = drms_open_recordset(drms_env, inRecQuery, &status);
   if (inrecs)
   {
      int nRecs = inrecs->n;
      int iset;
      int nsets = drms_recordset_getnumss(inrecs);
      DRMS_Record_t *inrec = NULL;
      DRMS_Record_t *orec = NULL;

      if (drms_env->verbose)
      {
	 (*printhst)("There are %d datasets.\n", nsets);
      }

      for (iset = 0; iset < nsets; iset++)
      {
	 while ((inrec = drms_recordset_fetchnextinset(drms_env, inrecs, &iset, &status)) != NULL)
	 {
	    /* create the corresponding output record (which may be 'blank) */
	    orec = drms_create_record(drms_env, outseries, DRMS_PERMANENT, &status);

	    /* There was some test for VDS_FITS_MERGE - not sure how to check for this in DRMS */
	    qual = drms_getkey_int(inrec, kQUALITY, &status);
	    if (status || qual != kQUAL_DATA_AVAILABLE)
	    {
	       errlog_staterr(printerr, "bad quality", !status ? qual : status, inseries, inrec->recnum);
	       CreateBlankRecord(inrec, orec);
	       continue;
	    }

	    // mapped_lmax = sds_getkey_int(in_sds, "MAPMMAX");
	 }
      }
   }

   char *segnamein = cmdparams_get_str(&cmdparams, kParamSegIn, NULL);
   char *segnameout = cmdparams_get_str(&cmdparams, kParamSegOut, NULL);

   if (inseries)
   {
      free(inseries);
   }

   return error;
} /* DoIt */

//errlog_paramerr(printerr, pname, err, series, recnum);
