/**

synframe sdate='2010.08.12' stime='16:00:00' in='su_yang.fd_M12m_remap_los' out='su_yang.synframe_los' synoptic='su_yang.hmi_M12m_synop' drmethod='Snodgrass' magresoln=1 synresoln=1 MAXMISSVALS=0 xx1=30

*/

// -----------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "jsoc_main.h"
#include "/home0/yliu/cvs/JSOC/proj/myproj/apps/src/fstats.c"

#include <mkl_blas.h>^M
#include <mkl_service.h>^M
#include <mkl_lapack.h>^M
#include <mkl_vml_functions.h>^M
#include <omp.h>^M
#include "fresize.c"

char *module_name = "dailysynframe";

#define DIE(msg) {fflush(stdout); fprintf(stderr, "%s, status = %d \n", msg, status); return(status);}
#define PARAMETER_ERROR(PNAME)
#define     PI     4.0 * atan(1.0)

ModuleArgs_t module_args[] =
  {
    {ARG_STRING, "in", "NOTSPECIFIED", "in"},
    {ARG_STRING, "out", "NOTSPECIFIED", "out"},
    {ARG_STRING, "synoptic", "NOTSPECIFIED", "synoptic"},
    {ARG_STRING, "smallsyn", "NOTSPECIFIED", ""},
    {ARG_STRING, "drmethod", "NOCORRDIFFROT", "drmethod"},
    {ARG_INT, "magresoln", "NOTSPECIFIED", "magresoln"},
    {ARG_INT, "synresoln", "NOTSPECIFIED", "synresoln"},
    {ARG_INT, "xx1", "-1", "xx1"},
    {ARG_INT, "yy1", "-1", "yy1"},
    {ARG_END}
  };

  float zgrid(int jph, int ith, int cmp, int rup, int dbl, 
       float phd[], float thd[], float lad[], float cth[], 
       float sth[], float csc[], float scs[]);
void frebinbox(float *image_in, float *image_out, int nx, int ny, int nbinx, int nbiny);

int DoIt(void)
{
  DRMS_RecordSet_t *inRS, *inRSfinal;
  DRMS_Record_t *inRec, *inRecfinal, *outRec = NULL, *smallRD; 
  DRMS_Segment_t *inSeg, *inSegfinal, *outSeg, *smalloutSeg;
  DRMS_Array_t *inArray, *inArrayfinal;
  TIME t_rec, t_rec0;
  TIME halfw = 7200.0; //half window = 2.0 hours
  char *t_window = "240m"; //hard-coded T-window
  char *inQueryfinal, *trec_str = NULL, *smallRecQuery;
  float crlt, crln;

  float clog0;
  float aa, bb, cc, dd, ee;
  int status = DRMS_SUCCESS, nrecs, irec;
  int i, j, crn;
  int xdim_syn, ydim_syn, xmg, ymg;
  char sdatestime[100], sdatetmp[100], timetmp[100];
  memset(sdatetmp, 0, sizeof(sdatetmp));
  memset(timetmp, 0, sizeof(timetmp));
  memset(sdatestime, 0, sizeof(sdatestime));

  char *sdate, *stime, *inQuery, *outQuery, *synQuery;
  char *drmethod;
  int magresoln = params_get_int(&cmdparams, "magresoln");
  int synresoln = params_get_int(&cmdparams, "synresoln");
  int xx1 = params_get_int(&cmdparams, "xx1");
  int yy1 = params_get_int(&cmdparams, "yy1");
  int nbin = 5;

  inQuery = (char *)params_get_str(&cmdparams, "in");
  outQuery = (char *)params_get_str(&cmdparams, "out");
  synQuery = (char *)params_get_str(&cmdparams, "synoptic");
  smallRecQuery = (char *)cmdparams_get_str(&cmdparams, "smallsyn", &status);
  drmethod = (char *)params_get_str(&cmdparams, "drmethod");

     aa = 13.1988; bb = 0.0; cc = 0.0;
  if (strcmp(drmethod, "Meunier") == 0)  //  Meunier
    {
     aa = 13.562; bb = -2.04; cc = -1.4875;
    }

  if (strcmp(drmethod, "Phil") == 0)     // Phil, ApJ 241:811-819, 1980
   {
    aa = 2.917; bb = -0.40; cc = -0.40;  //in mu rad/s;sidereal
    dd = 0.202006;                       // 1 deg/day=0.202006 mu rad/s
    aa = aa/dd; bb = bb/dd; cc=cc/dd;    // in deg / day
    ee = 0.930505;                       // synodic/sidereal=0.930505
    aa = aa*ee; bb = bb*ee; cc = cc*ee;  // in synodic
   }                                     // end "phil"

  if (strcmp(drmethod, "Pevtsov") == 0)
    {
     aa = 13.51; bb = -1.72;  cc = -2.31;
    }

  if (strcmp(drmethod, "Snodgrass") == 0)
   {
    aa = 2.897; bb = -0.339; cc = -0.485;//in mu rad/s,sidereal
    dd = 0.202006;                       // 1 deg/day=0.202006 mu rad/s
    aa = aa/dd;  bb=bb/dd; cc=cc/dd;     // in deg / day
    ee = 0.930505;                       // synodic/sidereal=0.930505
    aa = aa*ee; bb = bb*ee; cc = cc*ee;  //in synodic
   }

  outRec = drms_create_record(drms_env, outQuery, DRMS_PERMANENT, &status);
  if (status) DIE("Output recordset not created");
  smallRD = drms_create_record(drms_env, smallRecQuery, DRMS_PERMANENT, &status);
  if (status) DIE("Output record not created");

  inRS = drms_open_records(drms_env, inQuery, &status);
  if (status || inRS->n == 0)
     DIE("No input data found");
  inRec = inRS->records[0];
  t_rec = drms_getkey_time(inRec, "T_REC", &status);
  t_rec0 = t_rec - halfw;

  inQueryfinal = (char *)malloc(100 * sizeof(char));
  trec_str = (char *)malloc(30 * sizeof(char));
  sprint_time(trec_str, t_rec0, "TAI", 0);
  sprintf(inQueryfinal, "%s[%s/%s]", inRec->seriesinfo->seriesname, trec_str, t_window);
  printf("%s\n", inQueryfinal);
  drms_close_records(inRS, DRMS_FREE_RECORD);
  free(trec_str);

  inRSfinal = drms_open_records(drms_env, inQueryfinal, &status);
  if (status || inRS->n == 0) DIE("No input data found");
  nrecs = inRSfinal->n;
  int nref = (int)(nrecs/2);
  inRecfinal = inRSfinal->records[nref];
  t_rec = drms_getkey_time(inRecfinal, "T_REC", &status);
  t_rec0 = drms_getkey_time(inRecfinal, "T_OBS", &status);
  crn = drms_getkey_int(inRecfinal, "CAR_ROT", &status);
  crlt = drms_getkey_float(inRecfinal, "CRLT_OBS", &status);
  crln = drms_getkey_float(inRecfinal, "CRLN_OBS", &status);
  clog0 = drms_getkey_float(inRecfinal, "CRVAL1", &status);

//  drms_copykey(outRec, inRecfinal, "CADENCE");
  drms_copykey(outRec, inRecfinal, "DATASIGN");
  drms_copykey(outRec, inRecfinal, "DSUN_OBS");
  drms_copykey(outRec, inRecfinal, "OBS_VR");
  drms_copykey(outRec, inRecfinal, "OBS_VW");
  drms_copykey(outRec, inRecfinal, "OBS_VN");
  drms_copykey(outRec, inRecfinal, "QUALITY");
  drms_copykey(outRec, inRecfinal, "MAPMMAX");
  drms_copykey(outRec, inRecfinal, "MAPLGMAX");
  drms_copykey(outRec, inRecfinal, "MAPLGMIN");
  drms_copykey(outRec, inRecfinal, "SINBDIVS");
  drms_copykey(outRec, inRecfinal, "MAPBMAX");
  drms_copykey(outRec, inRecfinal, "MAPRMAX");
  drms_copykey(outRec, inRecfinal, "LGSHIFT");
  drms_copykey(outRec, inRecfinal, "INTERPO");
  drms_copykey(outRec, inRecfinal, "MCORLEV");
  drms_copykey(outRec, inRecfinal, "MOFFSET");
  drms_copykey(outRec, inRecfinal, "CARSTRCH");
  drms_copykey(outRec, inRecfinal, "DIFROT_A");
  drms_copykey(outRec, inRecfinal, "DIFROT_B");
  drms_copykey(outRec, inRecfinal, "DIFROT_C");
  drms_copykey(outRec, inRecfinal, "INSTRUME");
  drms_copykey(outRec, inRecfinal, "BLD_VERS");

//  smallRD;
int itmp;
float ftmp;
double dtmp;

  drms_copykey(smallRD, inRecfinal, "DATASIGN");
  drms_copykey(smallRD, inRecfinal, "DSUN_OBS");
  drms_copykey(smallRD, inRecfinal, "OBS_VR");
  drms_copykey(smallRD, inRecfinal, "OBS_VW");
  drms_copykey(smallRD, inRecfinal, "OBS_VN");
  drms_copykey(smallRD, inRecfinal, "QUALITY");
  itmp = drms_getkey_int(inRecfinal, "MAPMMAX", &status);
  drms_setkey_int(smallRD, "MAPMMAX", itmp/nbin);
  drms_copykey(smallRD, inRecfinal, "MAPLGMAX");
  drms_copykey(smallRD, inRecfinal, "MAPLGMIN");
  ftmp = drms_getkey_int(inRecfinal, "SINBDIVS", &status);
  drms_setkey_float(smallRD, "SINBDIVS", ftmp/(nbin-1));
  drms_copykey(smallRD, inRecfinal, "MAPBMAX");
  drms_copykey(smallRD, inRecfinal, "MAPRMAX");
  drms_copykey(smallRD, inRecfinal, "LGSHIFT");
  drms_copykey(smallRD, inRecfinal, "INTERPO");
  drms_copykey(smallRD, inRecfinal, "MCORLEV");
  drms_copykey(smallRD, inRecfinal, "MOFFSET");
  drms_copykey(smallRD, inRecfinal, "CARSTRCH");
  drms_copykey(smallRD, inRecfinal, "DIFROT_A");
  drms_copykey(smallRD, inRecfinal, "DIFROT_B");
  drms_copykey(smallRD, inRecfinal, "DIFROT_C");
  drms_copykey(smallRD, inRecfinal, "INSTRUME");
  drms_copykey(smallRD, inRecfinal, "BLD_VERS");

  inSeg = drms_segment_lookupnum(inRecfinal, 0);
  inArray = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
  if (status) DIE("No input data found");

  float *aveData;
  int naxis = inArray->naxis;
  xmg = inArray->axis[0]; ymg = inArray->axis[1];
  xdim_syn = 3600; ydim_syn = 1440;
  int inDims[2] = {xmg, ymg};
  int dxsz = 2 * inDims[0];     // jph in IDL, zgrid.pro
  int ith = inDims[1];
  int ppd = xdim_syn/360;       // pixels per degree
  int count = 0, xbeg = 30;
  if (xx1 == -1) xx1 = 60;              // in degrees
  if (yy1 == -1) yy1 = 0;                  // in pixels
  int hwd = xx1;       // in degree
  xx1 *= ppd;          // in pixels
  xbeg *= ppd;
  aveData = (float *)malloc(xmg * ymg * sizeof(float));
  int ii, jj;
  int smallDims[2], xout, yout;
  xout = xdim_syn/nbin; yout = ydim_syn/(nbin-1);
  smallDims[0] = xout; smallDims[1] = yout;

  for (i = 0; i < nrecs; i++)
    {
        inRecfinal = inRSfinal->records[i];
        inSeg = drms_segment_lookupnum(inRecfinal, 0);
        inArray = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
        if (status) 
           {
              printf(" No data file found, status=%d\n", status);
              drms_free_array(inArray);
              continue;
            }
        count += 1;
        float *inData = (float *)inArray->data;
        int crnn = drms_getkey_int(inRecfinal, "CAR_ROT", &status);
        float clogn = drms_getkey_float(inRecfinal, "CRVAL1", &status);
        int xshift = (int)(ppd * ((clogn - clog0) - 360.0 * (crnn - crn)));
//printf("count=%d, xshift=%d\n", i, xshift);
        for (jj = 0; jj < ymg; jj++)
            {
            for (ii = xbeg; ii < xmg - xbeg; ii++)
                {
                  aveData[jj * xmg + ii] += inData[jj * xmg + ii - xshift];
                }
            }
        drms_free_array(inArray);
    }
   drms_close_records(inRSfinal, DRMS_FREE_RECORD);
   for (jj = 0; jj < ymg; jj++)
       for (ii = 0; ii < xmg; ii++)
           aveData[jj * xmg + ii] /= count;

  float thd[ith], csc[ith], phd[dxsz];
  float lad[ith], sth[ith], cth[ith], scs[dxsz];
  double dtor = PI/180.;
  float constrate = 13.1988; // Carrington rotation rate in degree/day
  float ratelow, sinlat;
  DRMS_RecordSet_t *synRS = NULL;
  DRMS_Record_t *synRec;
  DRMS_Segment_t *synSeg;
  DRMS_Array_t *synArray, *supsynArray, *outArray, *smalloutArray;
  float *synData, *supsynData, *outData, *smalloutData;
  int i1, j1;
  int synleftst = ppd * hwd, supleftst = ppd * clog0;

  zgrid(dxsz, ith, 0, 0, 0, phd, thd, lad, cth, sth, csc, scs);

  snprintf(timetmp, sizeof(timetmp), "%s[%d/3]", synQuery, crn-2);
  synQuery = timetmp;
  printf("inputname= %s\n", synQuery);
  synRS = drms_open_records(drms_env, synQuery, &status);
  if (status || synRS->n == 0) DIE("No input data found");
                     // start combining the synoptic charts
  int nds = synRS->n;
  int supxDim = nds * xdim_syn;
  int supsynDim[2] = {supxDim, ydim_syn};
  int synDim[2] = {xdim_syn, ydim_syn};
  supsynArray = drms_array_create(DRMS_TYPE_FLOAT, 2, supsynDim, NULL, &status);
  supsynData = supsynArray->data;

  outArray = drms_array_create(DRMS_TYPE_FLOAT, 2, synDim, NULL, &status);
  outData = outArray->data;
  smalloutArray = drms_array_create(DRMS_TYPE_FLOAT, 2, smallDims, NULL, &status);
  smalloutData = smalloutArray->data;

  for (i = 0; i < nds; i++)
    {
       synRec = synRS->records[i];
       synSeg = drms_segment_lookupnum(synRec, 0);
       synArray = drms_segment_read(synSeg, DRMS_TYPE_FLOAT, &status);
       synData = synArray->data;
       int ii = (nds - 1 - i) * xdim_syn;
       for (j1 = 0; j1 < ydim_syn; j1++)
          for (i1 = 0; i1 < xdim_syn; i1++)
            {
               supsynData[supxDim * j1 + ii + i1] = synData[xdim_syn * j1 + i1];
            }                  
     }
                   // end of the combination: the super synoptic map--supsynData

                   // start to generate the right portion of synchronic map
                   // the differential rotation is taken into account.
   for (j = 0; j < ydim_syn/2; j++)
     {
        sinlat = sin(lad[j]*dtor);
        ratelow = aa + bb * sinlat * sinlat + cc * pow(sinlat, 4); // in degree/day

        for (i = synleftst; i < xdim_syn; i++)
        {
             float lon = (i - synleftst) * constrate/ratelow;
             float delta_x = lon - (int)lon;
             int lonpixel = (int)(supleftst + lon);
             float x1 = supsynData[j * supxDim + lonpixel];
             float x2 = supsynData[j * supxDim + lonpixel + 1];
             outData[j*xdim_syn + i] = (1.0 - delta_x) * x1 + delta_x * x2;
             jj = ydim_syn - 1 - j;
             x1 = supsynData[jj * supxDim + lonpixel];
             x2 = supsynData[jj * supxDim + lonpixel + 1];
             outData[jj*xdim_syn + i] = (1.0 - delta_x) * x1 + delta_x * x2;
        }
      }
    
    int magleft = ppd * (90 - hwd);
    for (j = 0; j < ydim_syn; j++)
        for (i = 0; i < 2*synleftst; i++)
            {
                outData[j * xdim_syn + i] = aveData[j * inDims[0] + (int)magleft + i];
            }

  double statMin, statMax, statMedn, statMean, statSig, statSkew, statKurt;
  int statNgood;
  if (fstats(xdim_syn*ydim_syn, outData, &statMin, &statMax, &statMedn, &statMean, &statSig,
      &statSkew, &statKurt, &statNgood)) printf("\n Statistics computation failed\n");

  frebinbox(outData, smalloutData, xdim_syn, ydim_syn, nbin, nbin-1);
  double smallstatMin, smallstatMax, smallstatMedn, smallstatMean, smallstatSig, smallstatSkew, smallstatKurt;
  int smallstatNgood;
  if (fstats(xout*yout, smalloutData, &smallstatMin, &smallstatMax, &smallstatMedn, &smallstatMean, &smallstatSig,
      &smallstatSkew, &smallstatKurt, &smallstatNgood)) printf("\n Statistics computation failed\n");


//    outRec = drms_create_record(drms_env, outQuery, DRMS_PERMANENT, &status);
//    if (status) DIE("Output recordset not created");

    drms_setkey_time(outRec, "T_REC", t_rec);
    drms_setkey_time(outRec, "T_OBS", t_rec0);
    drms_setkey_int(outRec, "CAR_ROT", crn);
    drms_setkey_float(outRec, "CRLT_OBS", crlt);
    drms_setkey_float(outRec, "CRLN_OBS", crln);
    float loncen = xdim_syn/2 + 0.5;
    drms_setkey_int(outRec, "CRPIX1", loncen);
    float latcen = ydim_syn/2 + 0.5;
    drms_setkey_double(outRec, "CRPIX2", latcen);
    float lonstep = -360.0/xdim_syn;
    drms_setkey_float(outRec, "CDELT1", lonstep);
    float latstep = 2.0/ydim_syn;
    drms_setkey_float(outRec, "CDELT2", latstep);
    double lonfirst = 360.0 * crn - clog0 + hwd - 360.0 + 360.0/xdim_syn;
    drms_setkey_double(outRec, "LON_FRST", lonfirst);
    double lonlast = 360.0 * crn - clog0 + hwd;
    drms_setkey_double(outRec, "LON_LAST", lonlast);
    drms_setkey_float(outRec, "LON_STEP", lonstep);

// synoptic map info
    float hwnwidth = drms_getkey_float(synRec, "HWNWIDTH", &status);
    drms_setkey_float(outRec, "HWNWIDTH", hwnwidth);
    float eqpoints = drms_getkey_float(synRec, "EQPOINTS", &status);
    drms_setkey_float(outRec, "EQPOINTS", eqpoints);
    float syndro_a = drms_getkey_float(synRec, "DIFROT_A", &status);
    drms_setkey_float(outRec, "OSYNDR_A", syndro_a);
    float syndro_b = drms_getkey_float(synRec, "DIFROT_B", &status);
    drms_setkey_float(outRec, "OSYNDR_B", syndro_b);
    float syndro_c = drms_getkey_float(synRec, "DIFROT_C", &status);
    drms_setkey_float(outRec, "OSYNDR_C", syndro_c);
// image statistics
    i = xdim_syn*ydim_syn;
    drms_setkey_int(outRec, "TOTVALS", i);
    drms_setkey_int(outRec, "DATAVALS", statNgood);
    i = xdim_syn*ydim_syn-statNgood;
    drms_setkey_int(outRec, "MISSVALS", i);
    drms_setkey_double(outRec, "DATAMIN", statMin);
    drms_setkey_double(outRec, "DATAMAX", statMax);
    drms_setkey_double(outRec, "DATAMEDN", statMedn);
    drms_setkey_double(outRec, "DATAMEAN", statMean);
    drms_setkey_double(outRec, "DATARMS", statSig);
    drms_setkey_double(outRec, "DATASKEW", statSkew);
    drms_setkey_double(outRec, "DATAKURT", statKurt);
// frame info.  
    drms_setkey_string(outRec, "FRTIMWDN", t_window); 
    drms_setkey_string(outRec, "SYNDRORA", drmethod);
    drms_setkey_int(outRec, "FRAVEPNT", nrecs);
    drms_setkey_float(outRec, "FRWINDOW", 2.0 * hwd);
    drms_setkey_float(outRec, "SYNDRO_A", aa);
    drms_setkey_float(outRec, "SYNDRO_B", bb);
    drms_setkey_float(outRec, "SYNDRO_C", cc);

    drms_keyword_setdate(outRec);

       printf("       WRITING OUTPUT\n");

    outSeg = drms_segment_lookupnum(outRec, 0);
    outArray->parent_segment = outSeg;
    status = drms_segment_write(outSeg, outArray, 0);
    if (status) DIE("problem writing file");

// writting the small size map 
    drms_setkey_time(smallRD, "T_REC", t_rec);
    drms_setkey_time(smallRD, "T_OBS", t_rec0);
    drms_setkey_int(smallRD, "CAR_ROT", crn);
    drms_setkey_float(smallRD, "CRLT_OBS", crlt);
    drms_setkey_float(smallRD, "CRLN_OBS", crln);
    loncen = xout/2 + 0.5;
    drms_setkey_int(smallRD, "CRPIX1", loncen);
    latcen = yout/2 + 0.5;
    drms_setkey_double(smallRD, "CRPIX2", latcen);
    lonstep = -360.0/xout;
    drms_setkey_float(smallRD, "CDELT1", lonstep);
    latstep = 2.0/yout;
    drms_setkey_float(smallRD, "CDELT2", latstep);
    lonfirst = 360.0 * crn - clog0 + hwd - 360.0 + 360.0/xout;
    drms_setkey_double(smallRD, "LON_FRST", lonfirst);
    lonlast = 360.0 * crn - clog0 + hwd;
    drms_setkey_double(smallRD, "LON_LAST", lonlast);
    drms_setkey_float(smallRD, "LON_STEP", lonstep);
// synoptic map info
    hwnwidth = drms_getkey_float(synRec, "HWNWIDTH", &status);
    drms_setkey_float(smallRD, "HWNWIDTH", hwnwidth);
    eqpoints = drms_getkey_float(synRec, "EQPOINTS", &status);
    drms_setkey_float(smallRD, "EQPOINTS", eqpoints);
    syndro_a = drms_getkey_float(synRec, "DIFROT_A", &status);
    drms_setkey_float(smallRD, "OSYNDR_A", syndro_a);
    syndro_b = drms_getkey_float(synRec, "DIFROT_B", &status);
    drms_setkey_float(smallRD, "OSYNDR_B", syndro_b);
    syndro_c = drms_getkey_float(synRec, "DIFROT_C", &status);
    drms_setkey_float(smallRD, "OSYNDR_C", syndro_c);
// image statistics
    i = xout*yout;
    drms_setkey_int(smallRD, "TOTVALS", i);
    drms_setkey_int(smallRD, "DATAVALS", smallstatNgood);
    i = xout*yout-smallstatNgood;
    drms_setkey_int(smallRD, "MISSVALS", i);
    drms_setkey_double(smallRD, "DATAMIN", smallstatMin);
    drms_setkey_double(smallRD, "DATAMAX", smallstatMax);
    drms_setkey_double(smallRD, "DATAMEDN", smallstatMedn);
    drms_setkey_double(smallRD, "DATAMEAN", smallstatMean);
    drms_setkey_double(smallRD, "DATARMS", smallstatSig);
    drms_setkey_double(smallRD, "DATASKEW", smallstatSkew);
    drms_setkey_double(smallRD, "DATAKURT", smallstatKurt);
// frame info.  
    drms_setkey_string(smallRD, "FRTIMWDN", t_window);
    drms_setkey_string(smallRD, "SYNDRORA", drmethod);
    drms_setkey_int(smallRD, "FRAVEPNT", nrecs);
    drms_setkey_float(smallRD, "FRWINDOW", 2.0 * hwd);
    drms_setkey_float(smallRD, "SYNDRO_A", aa);
    drms_setkey_float(smallRD, "SYNDRO_B", bb);
    drms_setkey_float(outRec, "SYNDRO_C", cc);

    drms_keyword_setdate(smallRD);

       printf("       WRITING OUTPUT--small-size version\n");

    smalloutSeg = drms_segment_lookupnum(smallRD, 0);
    smalloutArray->parent_segment = smalloutSeg;
    status = drms_segment_write(smalloutSeg, smalloutArray, 0);
    if (status) DIE("problem writing file");

    drms_free_array(smalloutArray);
    drms_free_array(outArray);
    drms_free_array(supsynArray);
    drms_free_array(synArray);

    drms_close_record(smallRD, DRMS_INSERT_RECORD);
    drms_close_record(outRec, DRMS_INSERT_RECORD);
    drms_close_records(synRS, DRMS_FREE_RECORD);
  return 0;
}        // end DoIt

/*                
 *            zgrid()
 *
 */
float zgrid(int jph, int ith, int cmp, int rup, int dbl, 
   float phd[], float thd[], float lad[], float cth[], 
   float sth[], float csc[], float scs[])
{

  double dtor = PI/180.;
  int i, j;

  float dcth = 2.0/ith;
  float dph = 360./jph;
  float cpc[jph], thr, dphr;
  float lthd[ith], lphd[jph], lscs[jph];
  float llad[ith], lsth[ith], lcth[ith];

#if 0
  if (dbl == 1)  // checking condition "dbl"
     {
      double thd[ith], phd[jph], csc[ith], cpc[jph];
      double dph = 360./jph;
      double dcth = 2.0/ith;
      double lad[ith], sth[ith], cth[ith], thr;
      double lthd[ith], lphd[jph], lscs[jph];
      double llad[ith], lsth[ith], lcth[ith];
      double dphr;
      float scs[jph];
     }  // end if "dbl" 
#endif

  for (i =0; i < ith; i++)
     {
      cth[i] = (i + 0.5) * dcth - 1.0;   // from south to north
      thr = acos(cth[i]);
      sth[i] = sin(thr);
      thd[i] = thr/dtor;
      lad[i] = 90-thd[i];
     }

  for (j = 0; j < jph; j++)
     phd[j] = (j + 0.5) * dph;          // from left to right

  if (cmp == 1)                   // checking condition "cmp"
      {
       for (i = 0; i < ith; i++)
         csc[i] = 1/sth[i];
       for (j = 0; j < jph; j++)
         {
          dphr = (phd[j] - cmp) * dtor;
          cpc[j] = cos(dphr);
          scs[j] = 1/cpc[j];
         }
      }   //end condition: cmp = 1

  if (rup == 1)   // condition "rup"
      {
       for (i = 0; i < ith; i++)
         {
          thd[i] = lthd[ith - i];
          lad[i] = llad[ith - 1];
          cth[i] = lcth[ith - i];
          sth[i] = lsth[ith - 1];
         }
       for (j = 0; j < jph; j++)
         {
          phd[j] = lphd[jph-j];
          scs[j] = lscs[jph-j];
         }
      } // end condition "rup"
  return 0;
} // end function

/* ################## Wrapper for Jesper's code ################## */

void frebin(float *image_in, float *image_out, int nx, int ny, int nbin)
{
  struct fresize_struct fresizes;
  int nxout, nyout;
  int nlead = nx;

  nxout = nx / nbin; nyout = ny / nbin;
  init_fresize_gaussian(&fresizes, (nbin / 2) * 2, (nbin / 2) * 2, nbin);
  fresize(&fresizes, image_in, image_out, nx, ny, nlead, nxout, nyout, nxout, (nbin / 2) * 2, (nbin / 2) * 2, DRMS_MISSING_FLOAT);
  free_fresize(&fresizes);

}

void frebinbox(float *image_in, float *image_out, int nx, int ny, int nbinx, int nbiny)
{
  int nxout, nyout;
  int ii, jj, i, j;
  nxout = nx / nbinx; nyout = ny / nbiny;
  for (j = 0; j < nyout; j++) {
    int yOff, jy;
    jy = j * nbiny;
    for (i = 0; i < nxout; i++) {
      int ix;
      ix = i * nbinx;
      float aveval = 0.0;
      int number = 0;
      for (jj = 0; jj < nbiny; jj++) {
        yOff = (jy + jj) * nx;
        for (ii = 0; ii < nbinx; ii++) {
          int iData = yOff + ix + ii;
          if (!(isnan(image_in[iData]))) {
            aveval += image_in[iData];
            number += 1;
          }
        }
      }
    image_out[j*nxout+i] = aveval/number;
    }
  }
}

// ************ END ********** END ******** END ********** END ***********
