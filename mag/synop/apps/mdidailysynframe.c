/**

mdidailysynframe in='su_yang.mdi_Blmap[1997.05.01_12:00:00_TAI/24h]' out='su_yang.mdi_Mldailysynframe' synoptic='su_yang.mdi_Mlsynop'

*/

// -----------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "jsoc_main.h"
#include "fstats.h"

char *module_name = "mdidailysynframe";

#define DIE(msg) {fflush(stdout); fprintf(stderr, "%s, status = %d \n", msg, status); return(status);}
#define PARAMETER_ERROR(PNAME)
#define     PI     4.0 * atan(1.0)

ModuleArgs_t module_args[] =
  {
    {ARG_STRING, "in", "NOTSPECIFIED", "in"},
    {ARG_STRING, "out", "NOTSPECIFIED", "out"},
    {ARG_STRING, "synoptic", "NOTSPECIFIED", "synoptic"},
    {ARG_STRING, "drmethod", "NOCORRDIFFROT", "drmethod"},
    {ARG_INT, "magresoln", "NOTSPECIFIED", "magresoln"},
    {ARG_INT, "synresoln", "NOTSPECIFIED", "synresoln"},
    {ARG_STRING, "qualmask", "0x00000200", ""},
    {ARG_INT, "xx1", "-1", "xx1"},
    {ARG_INT, "yy1", "-1", "yy1"},
    {ARG_END}
  };

  float zgrid(int jph, int ith, int cmp, int rup, int dbl, 
       float phd[], float thd[], float lad[], float cth[], 
       float sth[], float csc[], float scs[]);

int DoIt(void)
{
  DRMS_RecordSet_t *inRS, *inRSfinal;
  DRMS_Record_t *inRec, *inRecfinal, *outRec = NULL; 
  DRMS_Segment_t *inSeg, *inSegfinal, *outSeg;
  DRMS_Array_t *inArray, *inArrayfinal;
  TIME t_rec, t_rec0;
  TIME halfw = 43200.0; //half window = 12.0 hours
  char *t_window = "1440m"; //hard-coded T-window -- 24-hour
  char *inQueryfinal, *trec_str = NULL;
  float crlt, crln;

  double clog0;
  float aa, bb, cc, dd, ee;
  int status = DRMS_SUCCESS, nrecs, irec;
  int i, j, crn;
  int xdim_syn, ydim_syn, xmg, ymg;
  char sdatestime[100], sdatetmp[100], timetmp[100];
  unsigned int quality, qualvld, qualmask, qualmask1, qualmask2, qualmask3;

  memset(sdatetmp, 0, sizeof(sdatetmp));
  memset(timetmp, 0, sizeof(timetmp));
  memset(sdatestime, 0, sizeof(sdatestime));

  char *sdate, *stime, *inQuery, *outQuery, *synQuery;
  char *drmethod;
  int magresoln = params_get_int(&cmdparams, "magresoln");
  int synresoln = params_get_int(&cmdparams, "synresoln");
  int xx1 = params_get_int(&cmdparams, "xx1");
  int yy1 = params_get_int(&cmdparams, "yy1");

  inQuery = (char *)params_get_str(&cmdparams, "in");
  outQuery = (char *)params_get_str(&cmdparams, "out");
  synQuery = (char *)params_get_str(&cmdparams, "synoptic");
  drmethod = (char *)params_get_str(&cmdparams, "drmethod");
  qualmask = strtoul(cmdparams_get_str(&cmdparams, "qualmask", &status), (char **) 0, 0);
  qualvld = strtoul("0x00000201", (char **) 0, 0);

  qualmask1 = strtoul("0x00000000", (char **) 0, 0);
  qualmask3 = strtoul("0x00000200", (char **) 0, 0);
  qualmask2 = strtoul("0x00000201", (char **) 0, 0);

// printf("qualvld=%d\n", qualvld);

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

  inRS = drms_open_records(drms_env, inQuery, &status);
  if (status || inRS->n == 0) DIE("No input data found");
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

// find the middle data

  int count = 0;
  int *recp;
  int nref, rec_cen;
  recp = (int *)malloc(nrecs * sizeof(int));
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

        quality = drms_getkey_int(inRecfinal, "QUALITY", &status);
        if (drms_ismissing_int(quality))
        {
          printf("  Can't find QUALITY in log file; quality check disabled\n");
          continue;
        }
          else if ((quality != qualmask) & (quality != qualvld) & (quality != qualmask1))
          {
              printf("  Bad QUALITY = 0x%08x; rejected\n", quality);
              continue;
           }
        recp[count] = i;
        count += 1;
        drms_free_array(inArray);
    }
  
  if (count==0) DIE("No input remapped data found");
  rec_cen = (int)(count/2);
  nref = recp[rec_cen];
  printf("middle data id=%d\n", nref);
  inRecfinal = inRSfinal->records[nref];
  inSeg = drms_segment_lookupnum(inRecfinal, 0);
  inArray = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
  int naxis = inArray->naxis;
  xmg = inArray->axis[0]; ymg = inArray->axis[1];
  t_rec = drms_getkey_time(inRecfinal, "T_REC", &status);
  t_rec0 = drms_getkey_time(inRecfinal, "T_OBS", &status);
  crn = drms_getkey_int(inRecfinal, "CAR_ROT", &status);
  crlt = drms_getkey_float(inRecfinal, "CRLT_OBS", &status);
  crln = drms_getkey_float(inRecfinal, "CRLN_OBS", &status);
  clog0 = drms_getkey_double(inRecfinal, "CRVAL1", &status);

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

  float *aveData;
  xmg = inArray->axis[0]; ymg = inArray->axis[1];
  xdim_syn = 3600; ydim_syn = 1440;
  int inDims[2] = {xmg, ymg};
  int dxsz = 2 * inDims[0];     // jph in IDL, zgrid.pro
  int ith = inDims[1];
  int ppd = xdim_syn/360;       // pixels per degree
  int xbeg = 30;
  if (xx1 == -1) xx1 = 50;              // in degrees
  if (yy1 == -1) yy1 = 0;                  // in pixels
  int hwd = xx1;       // in degree
  TIME tobs_total = 0.0, tobs_ave;
  xx1 *= ppd;          // in pixels
  xbeg *= ppd;
  aveData = (float *)malloc(xmg * ymg * sizeof(float));
  int ii, jj;

  for (i = 0; i < count; i++)
    {
        inRecfinal = inRSfinal->records[recp[i]];
        inSeg = drms_segment_lookupnum(inRecfinal, 0);
        inArray = drms_segment_read(inSeg, DRMS_TYPE_FLOAT, &status);
        float *inData = (float *)inArray->data;
        int crnn = drms_getkey_int(inRecfinal, "CAR_ROT", &status);
        float clogn = drms_getkey_float(inRecfinal, "CRVAL1", &status);
        int xshift = (int)(ppd * ((clogn - clog0) - 360.0 * (crnn - crn)));
        TIME Tobs = drms_getkey_time(inRecfinal, "T_OBS", &status);
        quality = drms_getkey_int(inRecfinal, "QUALITY", &status);
        tobs_total += Tobs;
//printf("COUNT=%d, I=%d, RECP=%d, QUALITY = 0x%08x\n", count, i, recp[i], quality);
        for (jj = 0; jj < ymg; jj++)
            {
            for (ii = xbeg; ii < xmg - xbeg; ii++)
                {
                  aveData[jj * xmg + ii] += inData[jj * xmg + ii - xshift];
                }
            }
        drms_free_array(inArray);
    }

   tobs_ave = tobs_total/count;
   for (jj = 0; jj < ymg; jj++)
       for (ii = 0; ii < xmg; ii++)
           aveData[jj * xmg + ii] /= count;

  drms_close_records(inRSfinal, DRMS_FREE_RECORD);
  float thd[ith], csc[ith], phd[dxsz];
  float lad[ith], sth[ith], cth[ith], scs[dxsz];
  double dtor = PI/180.;
  float constrate = 13.1988; // Carrington rotation rate in degree/day
  float ratelow, sinlat;
  DRMS_RecordSet_t *synRS = NULL;
  DRMS_Record_t *synRec;
  DRMS_Segment_t *synSeg;
  DRMS_Array_t *synArray, *supsynArray, *outArray;
  float *synData, *supsynData, *outData;
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

  for (i = 0; i < nds; i++)
    {
       synRec = synRS->records[i];
       synSeg = drms_segment_lookup(synRec, "synopmag");
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

    drms_setkey_time(outRec, "T_REC", t_rec);
    trec_str = (char *)malloc(30 * sizeof(char));
    sprint_time(trec_str, tobs_ave, "TAI", 0);
    drms_setkey_time(outRec, "T_OBS", tobs_ave);
    drms_setkey_int(outRec, "CAR_ROT", crn);
    drms_setkey_float(outRec, "CRLT_OBS", crlt);
    drms_setkey_float(outRec, "CRLN_OBS", crln);
    double loncen = xdim_syn/2 + 0.0;
    drms_setkey_double(outRec, "CRPIX1", loncen);
       // origin is at the left corner of the first pixel
    double latcen = ydim_syn/2 + 0.0;
    drms_setkey_double(outRec, "CRPIX2", latcen);
       // origin is at the left corner of the first pixel
    double lonstep = -360.0/xdim_syn;
    drms_setkey_double(outRec, "CDELT1", lonstep);
    double latstep = 2.0/ydim_syn;
    drms_setkey_double(outRec, "CDELT2", latstep);
    double lonfirst = 360.0 * (crn - 1) - clog0 + hwd;
    drms_setkey_double(outRec, "LON_FRST", lonfirst);
    double lonlast = 360.0 * crn - clog0 + hwd - 360.0/xdim_syn;
    drms_setkey_double(outRec, "LON_LAST", lonlast);
    double loncenter = lonfirst + 180.0 - 0.5 * 360.0/xdim_syn;
    drms_setkey_double(outRec, "CRVAL1", loncenter);
    drms_setkey_double(outRec, "CARRTIME", loncenter);
    drms_setkey_double(outRec, "LON_STEP", lonstep);

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
    drms_setkey_int(outRec, "FRAVEPNT", count);
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

    free(recp); free(aveData); free(trec_str);
    drms_free_array(outArray);
    drms_free_array(supsynArray);
    drms_free_array(synArray);

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
// ************ END ********** END ******** END ********** END ***********
