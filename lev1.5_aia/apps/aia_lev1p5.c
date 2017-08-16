#include <string.h>
#include "jsoc_main.h"
#include "drms.h"
#include "drms_names.h"

#define NOT_SPECIFIED "***Not Specified***"
#define DIE(msg) {fprintf(stderr,"$$$$ %s: %s\n", module_name, msg); return 1;}

#define NBINS 1048576
static int hist[NBINS];

typedef enum { WL94A, WL131A, WL171A, WL193A, WL211A, WL304A, WL335A,
               WL1600A, WL1700A, WL4500A } wl_type;
int image_magrotate(void *, int, int, int, float, float, float, float,
                    void **, int *, int *, int, int);

ModuleArgs_t module_args[] =
{
  {ARG_STRING, "dsinp", NOT_SPECIFIED, "Input series query"},
  {ARG_STRING, "dsout", NOT_SPECIFIED, "Output series"},
  {ARG_STRING, "mpt", NOT_SPECIFIED, "Master Pointing Table series"},
  {ARG_STRING, "despike", "1", "remove cosmic ray hits"},
  {ARG_STRING, "rescale", "1", "rescale to fixed plate scale"},
  {ARG_STRING, "regrid", "1", "regrid type 0: nearest neighbor, 1: bicubic"},
  {ARG_STRING, "scale_to", "0.6", "rescale to fixed plate scale"},
  {ARG_STRING, "do_stretchmarks", "0", "fill in empty pixels created"},
  {ARG_STRING, "with_keys", "0", "write segment with all FITS keywords"},
  {ARG_STRING, "requestid", NOT_SPECIFIED, "Export request id if this program was invoked by the export system."},
  {ARG_FLAG, "h", "0", "Print usage message and quit"},
  {ARG_FLAG, "v", "0", "verbose flag"},
  {ARG_END}
};

char *module_name = "aia_lev1p5";
int despike = 1, rescale = 1, regridtype = 1, verbose = 0, do_stretchmarks = 0;
int withkeys = 0;
float scale_to = 0.6, scale_by;

int nice_intro(int help)
{
  int usage = cmdparams_get_int(&cmdparams, "h", NULL) != 0;
  verbose = cmdparams_get_int(&cmdparams, "v", NULL) != 0;
  if (usage || help) {
    printf("aia_lev1p5 {-h} {-v} dsinp=series_record_spec dsout=output_series\n"
        "  -h: print this message\n"
        "  -v: verbose\n"
        "despike=0, do not despike, default=1, despike\n"
        "rescale=0, do not rescale to fixed plate scale, default=1\n"
        "regrid=0, nearest neighbor, default=1, bicubic\n"
        "scale_to=0.6, rescale to fixed plate scale, default=0.6\n"
        "do_stretchmarks=1, fill in empty pixels created, default=0\n"
        "with_keys=1, write segment with all FITS keywords, default=0\n"
        "mpt=<master pointing table series>, default=NULL (use WCS keywords)\n"
        "dsinp=<recordset query> as <series>{[record specifier]} - required\n"
        "dsout=<series> - required\n");
    return(1);
  }
  return(0);
}
   
void sprint_time_ISO (char *tstring, TIME t)
{ 
  sprint_at(tstring,t);
  tstring[4] = tstring[7] = '-';
  tstring[10] = 'T';
  tstring[19] = '\0';
} 

static void FreeRecSpecParts(char ***snames, char ***filts, int nitems)
{
    if (snames)
    {
        int iname;
        char **snameArr = *snames;
        
        if (snameArr)
        {
            for (iname = 0; iname < nitems; iname++)
            {
                char *oneSname = snameArr[iname];
                
                if (oneSname)
                {
                    free(oneSname);
                }
            }
            
            free(snameArr);
        }
        
        *snames = NULL;
    }
    
    if (filts)
    {
        int ifilt;
        char **filtArr = *filts;
        
        if (filtArr)
        {
            for (ifilt = 0; ifilt < nitems; ifilt++)
            {
                char *onefilt = filtArr[ifilt];
                
                if (onefilt)
                {
                    free(onefilt);
                }
            }
            
            free(filtArr);
        }
        
        *filts = NULL;
    }
}

int ndx2wl(int ndx)
{
  static int wl[] = { 94, 131, 171, 193, 211, 304, 335, 1600, 1700, 4500 };
  if (-1 < ndx && ndx < 10) return wl[ndx];
  return -1;
}

int wl2ndx(int wl)
{
  switch (wl) {
    case   94: return WL94A;
    case  131: return WL131A;
    case  171: return WL171A;
    case  193: return WL193A;
    case  211: return WL211A;
    case  304: return WL304A;
    case  335: return WL335A;
    case 1600: return WL1600A;
    case 1700: return WL1700A;
    case 4500: return WL4500A;
    default: return -1;
  }
}

int DoIt ()
{
  int irec, iseg, nrecs, nsegs, status, is_aia=0;
  char *dsinp, *dsout, now_str[100];
  float crpix1, crpix2, cdelt1, cdelt2, crota2;
  float mag = 1.0, dx, dy;
  DRMS_Record_t *inprec, *outrec, *mptrec=NULL;
  DRMS_RecordSet_t *inprs, *mptrs=NULL;
  DRMS_Array_t *inparr=NULL, *outarr=NULL;
  DRMS_Segment_t *inpseg, *outseg;
  const char *reqid = NULL;
  const char *mpt = NULL;
  char seriesout[DRMS_MAXSERIESNAMELEN], mpr_str[256];
  char *allvers = NULL;
  char **sets = NULL;
  DRMS_RecordSetType_t *settypes = NULL; /* a maximum doesn't make sense */
  char **snames = NULL;
  char **filts = NULL;
  int nsets = 0;
  DRMS_RecQueryInfo_t rsinfo;
  TIME tstart, tstop;
  char mptqs[256];
  const char *x0names[] = { "A_094_X0", "A_131_X0", "A_171_X0", "A_193_X0",
    "A_211_X0", "A_304_X0", "A_335_X0", "A_1600_X0", "A1700_X0", "A_4500_X0"
  };
  const char *y0names[] = { "A_094_Y0", "A_131_Y0", "A_171_Y0", "A_193_Y0",
    "A_211_Y0", "A_304_Y0", "A_335_Y0", "A_1600_Y0", "A1700_Y0", "A_4500_Y0"
  };
  if (nice_intro(0)) return(0);

  dsinp = strdup(cmdparams_get_str(&cmdparams, "dsinp", NULL));
  dsout = strdup(cmdparams_get_str(&cmdparams, "dsout", NULL));
  mpt = strdup(cmdparams_get_str(&cmdparams, "mpt", NULL));
  if (strcmp(dsinp, NOT_SPECIFIED)==0) DIE("dsinp argument is required");
  if (strcmp(dsout, NOT_SPECIFIED)==0) DIE("dsout argument is required");
  if (strcmp(mpt, NOT_SPECIFIED)==0) mpt = NULL; 
  despike = cmdparams_get_int(&cmdparams, "despike", NULL);
  rescale = cmdparams_get_int(&cmdparams, "rescale", NULL);
  scale_to = cmdparams_get_float(&cmdparams, "scale_to", NULL);
  regridtype = cmdparams_get_int(&cmdparams, "regrid", NULL);
  do_stretchmarks = cmdparams_get_int(&cmdparams, "do_stretchmarks", NULL);
  withkeys = cmdparams_get_int(&cmdparams, "with_keys", NULL);
  reqid = cmdparams_get_str(&cmdparams, "requestid", NULL);
   
    if (strcmp(reqid, NOT_SPECIFIED) == 0)
    {
        reqid = NULL;
    }
    
    /* Parse output series name. */
    if (DRMS_SUCCESS != drms_record_parserecsetspec(dsout, &allvers, &sets, &settypes, &snames, &filts, &nsets, &rsinfo))
    {
        DIE("Invalid output record-set specification.");
    }
    
    if (nsets != 1)
    {
        FreeRecSpecParts(&snames, &filts, nsets);
        DIE("aia_lev1p5 supports writing to a single output series.");
    }
    
    snprintf(seriesout, sizeof(seriesout), "%s", snames[0]);
    FreeRecSpecParts(&snames, &filts, nsets);

  if (strstr(dsinp, "aia")) is_aia = 1;
  inprs = drms_open_records(drms_env, dsinp, &status);
  if (status) DIE("cant open recordset query");
  drms_stage_records(inprs, 1, 0);
  nrecs = inprs->n;
  printf("%d records\n", nrecs);
  for (irec=0; irec<nrecs; irec++) {
    int save_rec = 1, qualmask = 0x800100f0;
    inprec = inprs->records[irec];
    /* dsout isn't the correct argument to drms_create_record(). dsout is a record-set
     * specification, but this call takes a seriesname.
     outrec = drms_create_record(drms_env, dsout, DRMS_PERMANENT, &status); */
    outrec = drms_create_record(drms_env, seriesout, DRMS_PERMANENT, &status); 
    if (status) DIE("cant create recordset");
    status = drms_copykeys(outrec, inprec, 0, kDRMS_KeyClass_Explicit);
    if (status) DIE("Error in drms_copykeys()");
    {
       int fsn, quallev0, trec_off=0;
       double tr_step;
       double gaex_obs, gaey_obs, gaez_obs, haex_obs, haey_obs, haez_obs;
       long long tr_index;
       TIME t_rec, t_obs, tr_epoch;

       sprint_time_ISO(now_str, CURRENT_SYSTEM_TIME);
       drms_setkey_string(outrec, "DATE", now_str);
       t_obs = drms_getkey_time(inprec, "T_OBS", &status);
       if (status) DIE("T_OBS not found!");
       if (is_aia) {
         if ( 0 == drms_setkey_time(outrec, "T_REC", t_obs)) {
           tr_index = drms_getkey_longlong(outrec, "T_REC_index", &status);
           if (status) DIE("T_REC_index not found!");
           tr_step =  drms_getkey_double(outrec, "T_REC_step", &status);
           if (status) DIE("T_REC_step not found!");
           tr_epoch = drms_getkey_time(outrec, "T_REC_epoch", &status);
           if (status) DIE("T_REC_epoch not found!");
           t_rec = tr_epoch + tr_index*tr_step;
           drms_setkey_time(outrec, "T_REC", t_rec);
         }
         gaex_obs = drms_getkey_double(inprec, "GAEX_OBS", &status);
         if (status) {
           gaex_obs = drms_getkey_double(inprec, "GCIEC_X", &status);
           if(!status) drms_setkey_double(outrec, "GAEX_OBS", gaex_obs);
         }
         gaey_obs = drms_getkey_double(inprec, "GAEY_OBS", &status);
         if (status) {
           gaey_obs = drms_getkey_double(inprec, "GCIEC_Y", &status);
           if(!status) drms_setkey_double(outrec, "GAEY_OBS", gaey_obs);
         }
         gaez_obs = drms_getkey_double(inprec, "GAEZ_OBS", &status);
         if (status) {
           gaez_obs = drms_getkey_double(inprec, "GCIEC_Z", &status);
           if(!status) drms_setkey_double(outrec, "GAEZ_OBS", gaez_obs);
         }
         haex_obs = drms_getkey_double(inprec, "HAEX_OBS", &status);
         if (status) {
           haex_obs = drms_getkey_double(inprec, "HCIEC_X", &status);
           if(!status) drms_setkey_double(outrec, "HAEX_OBS", haex_obs);
         }
         haey_obs = drms_getkey_double(inprec, "HAEY_OBS", &status);
         if (status) {
           haey_obs = drms_getkey_double(inprec, "HCIEC_Y", &status);
           if(!status) drms_setkey_double(outrec, "HAEY_OBS", haey_obs);
         }
         haez_obs = drms_getkey_double(inprec, "HAEZ_OBS", &status);
         if (status) {
           haez_obs = drms_getkey_double(inprec, "HCIEC_Z", &status);
           if(!status) drms_setkey_double(outrec, "HAEZ_OBS", haez_obs);
         }
       }
        
        /* Set RequestID keyword, if it exists. */
        if (reqid && drms_keyword_lookup(outrec, "RequestID", 0))
        {
            if (DRMS_SUCCESS != drms_setkey_string(outrec, "RequestID", reqid))
            {
                DIE("Unable to set RequestID in output record.");
            }
        }
        
       if (t_obs < 0) save_rec = 0;
       if (mpt) {
         int wi, wl;
         if (mptrs == NULL) {
           snprintf(mptqs, 256, "%s[?%f BETWEEN T_START and T_STOP?]", 
                   mpt, t_obs);
           mptrs = drms_open_records(drms_env, mptqs, &status);
           if (status) DIE("Can't open MPT recordset");
           if(mptrs->n) mptrec = mptrs->records[0];
           else DIE("No MPT record found for T_OBS");
           tstart = drms_getkey_time(mptrec, "T_START", &status);
           tstop = drms_getkey_time(mptrec, "T_STOP", &status);
         } else if (tstart > t_obs || t_obs > tstop) {
           drms_close_records(mptrs, DRMS_FREE_RECORD);
           snprintf(mptqs, 256, "%s[?%f BETWEEN T_START and T_STOP?]",
                   mpt, t_obs);
           mptrs = drms_open_records(drms_env, mptqs, &status);
           if (status) DIE("Can't open MPT recordset");
           if(mptrs->n) mptrec = mptrs->records[0];
           else DIE("No MPT record found for T_OBS");
           tstart = drms_getkey_time(mptrec, "T_START", &status);
           tstop = drms_getkey_time(mptrec, "T_STOP", &status);
         }
         wl = drms_getkey_int(inprec, "WAVELNTH", &status);
         wi = wl2ndx(wl);
         crpix1 = drms_getkey_float(mptrec, x0names[wi], &status);
         crpix2 = drms_getkey_float(mptrec, y0names[wi], &status);
         sprintf(mpr_str, "%s[:#%lld]", mpt, mptrec->recnum);
         drms_setkey_string(outrec, "MPO_REC", mpr_str);
       } else {
         crpix1 = drms_getkey_float(inprec, "CRPIX1", &status);
         if (status) DIE("CRPIX1 not found!");
         crpix2 = drms_getkey_float(inprec, "CRPIX2", &status);
         if (status) DIE("CRPIX2 not found!");
       }
       cdelt1 = drms_getkey_float(inprec, "CDELT1", &status);
       if (status) DIE("CDELT1 not found!");
       cdelt2 = drms_getkey_float(inprec, "CDELT2", &status);
       if (status) DIE("CDELT2 not found!");
       if (rescale) {
         mag = fabs(cdelt1 / scale_to);
         cdelt1 /= mag;
         cdelt2 /= mag;
         drms_setkey_float(outrec, "CDELT1", cdelt1);
         drms_setkey_float(outrec, "CDELT2", cdelt2);
       }
       crota2 = drms_getkey_float(inprec, "CROTA2", &status);
       if (status) DIE("CROTA2 not found!");
       if (is_aia) {
         quallev0 = drms_getkey_int(inprec, "QUALLEV0", &status);
         if (status) DIE("QUALLEV0 not found!");
         if (quallev0 & qualmask) save_rec = 0;
       }
    }
    nsegs = hcon_size(&inprec->segments);
    for (iseg=0; iseg<1; iseg++) {
      void *output_array = NULL;
      int i, ix, iy, n, m, dtyp, nx, ny, npix;
      char *filename = NULL, *inpsegname = NULL;
      float *fval, tmpval, z = 16383.5;
      double s, s2, s3, s4, ss, datamin, datamax, datamedn, datamean;
      double data_rms, dataskew, datakurt, dtmp;
      inpseg = drms_segment_lookupnum(inprec, iseg);
      inpsegname = inpseg->info->name;
      filename = inpseg->filename;
      if (0 == iseg) {
        outseg = drms_segment_lookupnum(outrec, iseg);
        if (!outseg) DIE("Cant get output segment");
        inparr = drms_segment_read(inpseg, DRMS_TYPE_FLOAT, &status);
        if (status) DIE("drms_segment_read failed!");
        dtyp = 3;
        n = inparr->axis[0]; m = inparr->axis[1];
//      dx = crpix1 - (n + 1.0)*0.5; dy = crpix2 - (m + 1.0)*0.5;
        dx = (n + 1.0)*0.5 - crpix1; dy = (m + 1.0)*0.5 - crpix2;

        fval = inparr->data;
        for (ix = 0; ix<n; ix++) {
          fval[ix] = 0.0;
          fval[(m - 1)*n + ix] = 0.0;
        }
        for (iy = 0; iy<m; iy++) {
          fval[iy*n] = 0.0;
          fval[iy*n + n - 1] = 0.0;
        }
        status = image_magrotate( inparr->data, n, m, dtyp, crota2,
                 mag, dx, dy, &output_array, &nx, &ny,
                 regridtype, do_stretchmarks);
        if (status) DIE("image_magrotate failed!");
        /* if out dimen != inp dimen, inparr->axis below is incorrect */
        if (is_aia) {
          outarr = drms_array_create(DRMS_TYPE_INT, 2, inparr->axis,
                                   NULL, &status);
        } else {
          outarr = drms_array_create(DRMS_TYPE_INT, 2, inparr->axis,
                                   NULL, &status);
        }
        if (status) DIE("drms_array_create failed!");
        s = s2 = s3 = s4 = 0.0; datamin = 9.9e9; datamax = -9.9e9, npix = 0; 
        for (i=0; i<nx*ny; i++) {
          if (is_aia) {
            if (*((float *)(output_array)+i)<0) *((float *)(output_array)+i)=0;
            if (*((float *)(output_array)+i)>z) *((float *)(output_array)+i)=z;
            *((int *)(outarr->data)+i) = *((float *)(output_array)+i);
          } else {
            *((int *)(outarr->data)+i) = *((float *)(output_array)+i);
          }
          tmpval = *((float *)(output_array)+i);
          if (finite(tmpval)) {
            npix++;
            if (tmpval < datamin) datamin = tmpval;
            if (tmpval > datamax) datamax = tmpval;
            s += tmpval;
            s2 += tmpval*tmpval;
            s3 += tmpval*tmpval*tmpval;
            s4 += tmpval*tmpval*tmpval;
          }
        }
        drms_setkey_int(outrec, "DATAMIN", (int) datamin);
        drms_setkey_int(outrec, "DATAMAX", (int) datamax);
        s /= npix;
        drms_setkey_float(outrec, "DATAMEAN", s);
        ss = s*s;
        s2 /= npix;
        s3 /= npix;
        s4 /= npix;
        if (npix > 1) {
          dtmp = npix * (s2-ss) / (npix-1);
          data_rms = sqrt(dtmp);
          drms_setkey_float(outrec, "DATARMS", (float) data_rms);
          if (dtmp > 0.0) {
            dataskew = (s3 - s * (3*s2 - 2*ss)) / (dtmp*data_rms);
            drms_setkey_float(outrec, "DATASKEW", (float) dataskew);
            datakurt = (s4 - 4*s*s3 + 3*ss*(2*s2-ss)) / (dtmp*dtmp) - 3;
            drms_setkey_float(outrec, "DATAKURT", (float) datakurt);
          }
        }
        drms_setkey_float(outrec, "CROTA2", 0.0);
        drms_setkey_float(outrec, "CRPIX1", (n + 1.0)*0.5);
        drms_setkey_float(outrec, "CRPIX2", (m + 1.0)*0.5);
        drms_setkey_float(outrec, "X0", 0.0);
        drms_setkey_float(outrec, "Y0", 0.0);
        if (withkeys) {
          drms_segment_writewithkeys(outseg, outarr, 0);
         } else {
          drms_segment_write(outseg, outarr, 0);
        }
        if (inparr) drms_free_array(inparr);
        if (output_array) free(output_array);
        if (outarr) drms_free_array(outarr);
      }
    }
    if (save_rec) drms_close_record(outrec, DRMS_INSERT_RECORD);
    else drms_close_record(outrec, DRMS_FREE_RECORD);
  }

  return 0;
}
