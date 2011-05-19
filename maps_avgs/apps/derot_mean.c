/* derot_mean.c */

/* **************************************************************
 *
 * Module Name: derot_mean
 * Purpose:  Produces remapped, weight averages of input data sets
 * Programmer: John Beck    Beck@Sun.Stanford.EDU
 *
 * Notes: This is based on the DSDS module dr_mean (written by 
 * J. Beck & R. Bogart)
 *
 * Usage: derot_mean InData= OutData= [flags]
 * 
 * **************************************************************
 */

/**
\defgroup derot_mean derot_mean - Produce average of multiple images, removing rotational motion
@mainpage

@version

@author J.G. Beck


\par Synopsis:
\code 
derot_mean in=<seriesname> out=<seriesname> t_start=<time string> [t_end=<time string] step_t=<step> (other options described below)
\endcode

\details

\b Derot_mean produces an average of multiple images, removing the
translational motions due to solar rotation.  The averaging function,
rotation rate and output geometry are customizable.

The only flag is "-v" for turning on/off verbose reporting.

\par Options:

\par Parameters The parameters can be broken into three categories:
rotation model, averaging function and output geometry.  The parameters
are discussed below.

\par Parameters controlling the rotation model & meridional motion used.
\li \c a0=<floating number> - the equatorial rotation rate in units of
micro-Radians per second.  For a solid-body rotation model set only this term.
\li \c a2=<floating number> - a differential rotation term in unit of
micro-Radians per second.  This coeficient is multiplied by sin^2 latitude.
\li \c a4=<floating number> - a differential rotation term in unit of
micro-Radians per second.  This coeficient is multiplied by sin^4 latitude.
\li \c merid_v=<floating number> - the meridional circulation term, in unit of
meters per second; positive is northward flow.

\par Parameters controlling the weighting function used for the
temporal averaging.
\li \c wt_func=<string> - the weighting function used (see below for details)
\li \c wt_len=<float>   - the length of the weighting window, in seconds, 
\li \c wt_parm=<float>  - the parameter of the weighting function, typically
in units of seconds, although this parameter is function specific (see below)

\par Parameters controlling the output image geometry
\li \c cols=
\li \c rows=
\li \c proj=<projection name> old, MDI-style
\li \c crot=
\li \c crpix1=
\li \c crpix2=
\li \c crval1=
\li \c crval2=
\li \c cdelt1=
\li \c cdelt2=
\li \c ctype1=
\li \c ctype2=
\li \c ang_rsun=
\li \c dist=

\par Examples:

\b Example 1:
   To do an  average on hmi.V_45s, at 2010.09.09_12:00:00_TAI with default rotation, weighting and geometry, and store the result in the data series foo
\code
  derot_mean in=hmi.V_45s out=foo t_start=2010.09.09_12:00:00_TAI
\endcode

\b Example 2:
To do a series of one hour boxcar averages of hmi.M_45s data for 10 Oct 2010
spaced at 15 minute intervals,  
\code
 derot_mean  in="hmi.M_45s" out="su_foo.bar" \
 # derot_mean  in="su_beck.HGR1_test" out="su_beck.drmRotTest" \
   start_t="2010.10.10_00:00:00_TAI" end_t="2010.10.10_23:59:00_TAI" \
   step_t=15 wt_func="box" wt_len="3600" wt_parm="1800" \
   rows=4096 cols=4096 crpix1=2048 crpix2=2048 crval1=0.0 crval2=0.0 \
   cdelt1=0.03 cdelt2=0.03 v=1 


*/

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <complex.h>
#include <sys/time.h>
#include <jsoc_main.h>
//#include <drms_types.h>
#include "./derot_mean.h"
#include "./cartography.c"
#include "./ccinterp.c"
//  #include "./heliographic_coords.c"
#include <astro.h>
#include <assert.h>

char *module_name = "derot_mean";
char *version_id  = "0.0.1";

ModuleArgs_t module_args[] =
{
  {ARG_STRING,  "in",      "",       "Input data series."},
  {ARG_STRING,  "out",     "",       "Output data series"},
  {ARG_TIME,    "start_t", "",       "First output record time "},     // 1
  {ARG_TIME,    "end_t",   "-1",     "Last output record time (or upper limit)"},     
  {ARG_INT,     "step_t",  "900",    "Steps between records in seconds"},
  {ARG_INT,     "cols",    "-1",     "Number of columns in output image"},  
  {ARG_INT,     "rows",    "-1",     "Number of rows"},  
  {ARG_FLOAT,   "a0",      " 2.851", "Solid Body Rotation rate (uRad/s)"},
  {ARG_FLOAT,   "a2",      "0.0",    "A2 rotation term (uRad/S)"},  // set to zero for test
  {ARG_FLOAT,   "a4",      "0.0",    "A4 rotation term (uRad/s)"},  // set to zero for test
//  {ARG_FLOAT,   "a2",      "-0.343", "A2 rotation term (uRad/S)"},
//  {ARG_FLOAT,   "a4",      "-0.474", "A4 rotation term (uRad/s)"},
  {ARG_FLOAT,   "merid_v", "0.0",    "Meridional Flow in m/s"},         // 9
  {ARG_STRING,  "wt_func", "hath",   "Weighting Filter Function Name"}, // 10
  {ARG_INT,     "wt_len",  "1860.0", "Width of Weighting Window in Seconds"},
  {ARG_FLOAT,   "wt_parm", "900.0",  "Weighting Filter Parameter"},
  {ARG_STRING,  "proj",    "ortho",  "Projection Name (MDI style)"},
  {ARG_FLOAT,   "crot",    "0.0",    "WCS Coord Rotation Angle (-pangle)"},   
  {ARG_FLOAT,   "crpix1",  "None",    "Reference Pixel (Columnwise)"},     
  {ARG_FLOAT,   "crpix2",  "None",    "Reference Pixel (Rowwise)"},   
  {ARG_FLOAT,   "crval1",  "",    "Coord at Ref Pixel (deg/arcsec)"},
  {ARG_FLOAT,   "crval2",  "",    "Coord at Ref Pixel (deg/arcsec)"},
  {ARG_FLOAT,   "cdelt1",  "None",    "Pixel size at Ref Pix (in appropriate units)"},
  {ARG_FLOAT,   "cdelt2",  "None",    "Pixel size at Ref Pix (in appropriate units)"},
  {ARG_STRING,  "ctype1",  "CRLN-SIN", "Coordinate Name for col axis"},
  {ARG_STRING,  "ctype2",  "CRLT-SIN", "Coordinate Name for row axis"}, // 22
  {ARG_FLOAT,   "ang_rsun", "960.",  "Solar Radius in arcsecs (if AERIAL)"},
  {ARG_FLOAT,   "dist",    "1.0",    "Obs-Solar Distance in AU (if AERIAL)"}, 
//  {ARG_STRING,  "qualfile", "", "", ""},     Not used
  {ARG_FLAG,    "v",        "",      "Run in Verbose Mode"},
  {ARG_END}
};

ParmList parms;

/* These functions involve reading, checking cmdline params & subsequent calcs */
int    getcheck_cmd_parms(void);
int    get_cparms(void);
int    check_cparms(void);
int    cparm_getcheck_float(char kword[], double *out);
int    cparm_getcheck_int(char kword[], int *out);
int    cparm_getcheck_time(char kword[], TIME *out);
int    cparm_getcheck_str(char kword[], char out[]);
int    wcs_map_proj(void);
int    calc_weight_function(void);

// These functions involve setting up output arrays & db
int    setup_outdb(OutImgs **outimg, DRMS_RecordSet_t **outDB);
int    check_outdb(DRMS_Record_t *rec);
int    init_outimgs(OutImgs **outimg);
int    calc_out_geom(OutImgs outimg);

// These functions involve opening & testing input db for conformity
int    open_check_indb(DRMS_RecordSet_t **inDB);
int    get_query_str(char *qry);
int    check_indb(DRMS_Record_t *rec);

// These involve the actual derotation summation
int    addin(DRMS_RecordSet_t *inDB, OutImgs *out);
int    readin(DRMS_Record_t *irec, TIME *timg, double *Bo, double *Lo, 
	      double *Rsun, double *phi,
              double *xc, double *yc, char ctype1[], char ctype2[]);
int    rec_getcheck_double(DRMS_Record_t *rec, char kword[], double *out);
int    rec_getcheck_int(DRMS_Record_t *rec, char kword[], int *out);
int    rec_getcheck_str(DRMS_Record_t *rec, char kword[], char out[]);
int    rec_getcheck_time(DRMS_Record_t *rec, char kword[], TIME *out);

// These involve  writing out records
int    write_out_dbs(OutImgs *outimgs, DRMS_RecordSet_t *outDB, DRMS_RecordSet_t *inDB);
int    set_bscale_bzero(OutImgs *outimg);

// Some Utility & cleanup functions
double wcs_parm_unit(double parm, char *unit);
int    free_outimgs(OutImgs *outimg);
int    If_Err_Print(void);
void   V_print(char *str);

/* 
**********************************************************************
*
*    THIS IS THE DOIT FUNCTION WHERE THE TOP CONTROL TAKES PLACE
*
*    Note, this is where all error writing occurs.  Errors are passed
*    using the parms structure
*
**********************************************************************
*/

int DoIt(void) 
{
  int status, done;
  OutImgs *out_imgs;
  DRMS_RecordSet_t *inDB  = NULL;
  DRMS_RecordSet_t *outDB = NULL;
  
  V_print("\n===== Checking Commandline Arguments =====\n");
  if(getcheck_cmd_parms() !=  kMyMod_Success) // read command parms
    return (status = If_Err_Print());

  V_print("\n===== Checking Output Data Series =====\n");
  if(setup_outdb(&out_imgs, &outDB) != kMyMod_Success) // open & check output DB
    return (status = If_Err_Print());

  V_print("\n===== Checking Input Data Series =====\n");
  if(open_check_indb(&inDB) !=  kMyMod_Success ) // open & check  input DB
    return (status = If_Err_Print());

  V_print("\n===== Summing Derotated Mean =====\n");
  if(addin(inDB, out_imgs) !=  kMyMod_Success) // Add input images to output
    return (status = If_Err_Print());

  V_print("\n===== Writing Records to Output =====\n\n\n");
  if(write_out_dbs(out_imgs, outDB, inDB) !=  kMyMod_Success)  // Write out all output images
    return (status = If_Err_Print());

  V_print("Done.  derot_mean exiting normally\n");
  return kMyMod_Success;
}



/* 
 ******************************************************************* 
 * 
 *  Functions involving reading & calculating cmd-line arguments
 * 
 ******************************************************************** 
*/
//GETCHECK_CMD_PARMS
int  getcheck_cmd_parms(void)
{
  int status;

  strcpy(parms.SigMsg," ");
  strcpy(parms.Msg," ");
  if((status = get_cparms()))
    return status;

  if((status = wcs_map_proj()))
    return status;

  if((status = check_cparms()))
    return status;

  if((status = calc_weight_function()))
    return status;
  
  V_print("OK\n");
  return kMyMod_Success;    // Exit getcheck_cmd_params normally
}

// GET_CPARMS
int  get_cparms(void)
{
  int status;
  /* Have to allocate memory for the strings */
  
  /* 
     Note: this only returns an error if the reading is an error
     it does not check that the values are valid.  That is done
     later as it requires more complicated logic 
  */
  if((status = cparm_getcheck_str("in",(parms.iser))))
    return status;
  if((status = cparm_getcheck_str("out",(parms.oser))))
    return status;

  if((status = cparm_getcheck_time("start_t",&(parms.Tstart))))
    return status;
  if((status = cparm_getcheck_time("end_t",&(parms.Tend))))
    return status;
  if((status = cparm_getcheck_int("step_t",&(parms.Tstep))))
    return status;
  if((status = cparm_getcheck_int("rows",&(parms.rows))))
    return status;
  if((status = cparm_getcheck_int("cols",&(parms.cols))))
    return status;
  if((status = cparm_getcheck_float("a0",&(parms.A0))))
    return status;
  if((status = cparm_getcheck_float("a2",&(parms.A2))))
    return status;
  if((status = cparm_getcheck_float("a4",&(parms.A4))))
    return status;
  if((status = cparm_getcheck_float("merid_v",&(parms.Meri_V))))
    return status;

  if((status = cparm_getcheck_str("wt_func",(parms.WtFunc))))
    return status;
  if((status = cparm_getcheck_int("wt_len",&(parms.WtLen))))
    return status;
  if((status = cparm_getcheck_float("wt_parm",&(parms.WtParm))))
    return status;

  if((status = cparm_getcheck_str("proj",(parms.projname))))
    return status;
  if((status = cparm_getcheck_str("ctype1",(parms.CTYPE1))))
    return status;
  if((status = cparm_getcheck_str("ctype2",(parms.CTYPE2))))
    return status;
  if((status = cparm_getcheck_float("crot",&(parms.CROTA2))))
    return status;
  if((status = cparm_getcheck_float("crpix1",&(parms.CRPIX1))))
    return status;
  if((status = cparm_getcheck_float("crpix2",&(parms.CRPIX2))))
    return status;
  if((status = cparm_getcheck_float("crval1",&(parms.CRVAL1))))
    return status;
  if((status = cparm_getcheck_float("crval2",&(parms.CRVAL2))))
    return status;
  if((status = cparm_getcheck_float("cdelt1",&(parms.CDELT1))))
    return status;
  if((status = cparm_getcheck_float("cdelt2",&(parms.CDELT2))))
    return status;
  if((status = cparm_getcheck_float("ang_rsun",&(parms.ang_rad))))
    return status;
  if((status = cparm_getcheck_float("dist",&(parms.dist))))
    return status;

  parms.verbose = cmdparams_isflagset(&cmdparams,"v");

  return kMyMod_Success;     // exit get_cparms normally
}

// CHECK_CPARMS
int check_cparms(void)
{
  char temp[100];

  if(parms.Tend < 0)                   //  if Tend is not set, set it to Tstart
    parms.Tend = parms.Tstart;
  else if(parms.Tend < parms.Tstart)   // First check that Tend >= Tstart
  {
    fprintf(stderr,"Error: Tend must be greater  or equal to Tstart\n");
    return (parms.status = kMyMod_ValErr);
  }
  if (parms.Tstep == 0)
  {
    fprintf(stderr,"Error: Tstep must be greater or equal to zero\n");
    return (parms.status = kMyMod_ValErr);
  }
  else if (parms.Tstep < 0)
  {
    V_print("Warning: Tstep < 0.  Absolute value taken\n");
    parms.Tstep = abs(parms.Tstep);
  }

  // Then verify that required geometry values are valid
  // (Note: CRPIXs are checked in check_outdb() where they can be set
  // to 0.5*cols or rows

  if( isnan(parms.CROTA2))
    parms.CROTA2 = 0.0;

  if( isnan(parms.CRVAL1) || isnan(parms.CRVAL2) )
  {
    fprintf(stderr,"Error: Missing or Invalid CRVALs\n");
    return (parms.status = kMyMod_ValErr);
  }
  if( isnan(parms.CDELT1) || isnan(parms.CDELT2) )
  {
    fprintf(stderr,"Error: Missing or Invalid CDELTs\n");
    return (parms.status = kMyMod_ValErr);
  }
  if( isnan(parms.A0) || isnan(parms.A2) || 
      isnan(parms.A4) || isnan(parms.Meri_V) )
  {
    fprintf(stderr,"Error: Invalid Rotation Coefs or V_meridional\n");
    return (parms.status = kMyMod_ValErr);
  }

  if(  isnan(parms.WtParm) || parms.WtLen < 1)
  {
    fprintf(stderr,"Error: Invalid Weighting Function Parameters\n");
    return (parms.status = kMyMod_ValErr);
  }
 
  // Convert CRVALs & CDELTs to appropriate (Projection specific) units
  if(parms.projcode == AERIAL)
  {
    if(parms.CRVAL1 != 0 || parms.CRVAL2 != 0)
    {
      sprintf(temp,"Warning: AERIAL map requires CRVALn to be Zero - setting to Zero\n");
      V_print(temp);
      parms.CRVAL1 = parms.CRVAL2 = 0.0;
    }
    if(parms.ang_rad < 0 && parms.dist < 0)
    {
      sprintf(parms.Msg,"%sWarning: AERIAL map requires ang_rsun or dist0\n",parms.Msg);
      return (parms.status = kMyMod_ValErr);
    }
    if(parms.ang_rad <= 0)
      parms.ang_rad = asin( 1.0/(parms.dist*214.93950));
    if(parms.dist <= 0)
      parms.dist = 1.0 / (sin(parms.ang_rad)*214.93950);    

    parms.CDELT1 *= (M_PI / 180.0 / 3600.0);   // Convert from arcsec to rad
    parms.CDELT2 *= (M_PI / 180.0 / 3600.0);   // Convert from arcsec to rad

    if(parms.rsun <= 0  && (parms.CDELT1 == parms.CDELT2) )
      parms.rsun = parms.ang_rad / parms.CDELT1;
  }
  else if (parms.projcode == LAMBERT || parms.projcode == ORTHOGRAPHIC ||  
	   parms.projcode == STEREOGRAPHIC || parms.projcode == POSTEL || 
	   parms.projcode == GNOMONIC || parms.projcode == RECTANGULAR)
  {
    parms.CRVAL1 *= (M_PI / 180.0);  // Convert from deg to rad
    parms.CRVAL2 *= (M_PI / 180.0);  // Convert from deg to rad
    parms.CDELT1 *= (M_PI / 180.0);  // Convert from deg to rad
    parms.CDELT2 *= (M_PI / 180.0);  // Convert from deg to rad
  }
  else if (parms.projcode == CYLEQA || parms.projcode == SINEQA ||  
	   parms.projcode == MERCATOR || parms.projcode == CASSINI)
  {
    parms.CRVAL1 *= (M_PI / 180.0);  // Convert from deg to rad
    parms.CDELT1 *= (M_PI / 180.0);  // Convert from deg to rad
  }

  //  If projcode is CYLEQA, SINEQA, MERCATOR, CASSINI do not convert!!
  
  // Do some angle converstions & fill in some values
  parms.CROTA2  *= (M_PI / 180.0);  // Convert from deg to rad
  parms.pangle -1.0 * parms.CROTA2; // Change sign 
  parms.A0  =   parms.A0 / 1.0e6;   // convert from uRad/s to Rad/s
  parms.A2  =   parms.A2 / 1.0e6;   // convert from uRad/s to Rad/s
  parms.A4  =   parms.A4 / 1.0e6;   // convert from uRad/s to Rad/s
  parms.Meri_V  =   parms.Meri_V / 6.955e8;   // convert from m/s to Rad/s
  parms.size = parms.rows * parms.cols; // Calc number of pixels
  parms.Nmaps = ceil(((1 + parms.Tend - parms.Tstart)/parms.Tstep));
  if(parms.Nmaps  < 1) 
    parms.Nmaps = 1;

  return kMyMod_Success;     // exit check_cparms normally
}

// CPARM_GETCHECK_FLOAT
int cparm_getcheck_float(char kword[], double *out)
{
  int stat;
  double temp;
  temp = cmdparams_get_float(&cmdparams,kword,&stat);
  if(stat)
  {
    sprintf(parms.Msg,"%sError %d: could not read cmd param %s\n",parms.Msg,stat,kword);
    return (parms.status = stat);
  }
  (*out) = temp;
  return kMyMod_Success;     // exit cparm_getcheck_float normally

}

// CPARM_GETCHECK_INT
int cparm_getcheck_int(char kword[], int *out)
{
  int temp, stat;
  
  temp = cmdparams_get_int(&cmdparams,kword,&stat);
  if(stat)
  {
    sprintf(parms.Msg,"%sError %d: could not read cmd param %s\n",parms.Msg,stat,kword);
    return (parms.status = stat);
  }
  (*out) = temp;
  return kMyMod_Success;     // exit cparm_getcheck_int normally
  
}

// CPARM_GETCHECK_TIME
int cparm_getcheck_time(char kword[], TIME *out)
{
  int stat;
  TIME temp;
  temp = cmdparams_get_time(&cmdparams,kword,&stat);
  if(stat)
  {
    sprintf(parms.Msg,"%sError %d: could not read cmd param %s\n",parms.Msg,stat,kword);
    return (parms.status = stat);
  }
  (*out) = temp;
  return kMyMod_Success;     // exit cparm_getcheck_time normally
  
}

// CPARM_GETCHECK_STR
int cparm_getcheck_str(char kword[], char out[])
{
  int stat;
  const char *temp;
  temp = cmdparams_get_str(&cmdparams,kword,&stat);
  if(stat)
  {
    sprintf(parms.Msg,"%sError %d: could not read cmd param %s\n",parms.Msg,stat,kword);
    return (parms.status = stat);
  }
  strcpy(out,temp);
  return  kMyMod_Success;     // exit cparm_getcheck_str normally
  
}

// WCS_MAP_PROJ
int wcs_map_proj(void)
{
  int i, done;
  char comp1[12],comp2[12];
  static struct wcs_projnames 
  {
    int code;
    char CTYPE1[9];
    char CTYPE2[9];
    char projname[20];
  } wcs[]= {
    {RECTANGULAR,  "CRLN-CAR",  "CRLT-CAR", "Plate-Caree"},
    {RECTANGULAR,  "HGLN-CAR",  "HGLT-CAR", "Plate-Caree"},
    {CASSINI,      "CRLN-CAS",  "CRLT-CAS", "Cassini"},
    {CASSINI,      "HGLN-CAS",  "HGLT-CAS", "Cassini"},
    {MERCATOR,     "CRLN-MER",  "CRLT-MER", "Mercator"},
    {MERCATOR,     "HGLN-MER",  "HGLT-MER", "Mercator"},
    {CYLEQA,       "CRLN-CEA",  "CRLT-CEA", "CylEqA"},
    {CYLEQA,       "HGLN-CEA",  "HGLT-CEA", "CylEqA"},
    {SINEQA,       "CRLN-SFL",  "CRLT-SFL", "Sanson-Flamsteed"},
    {SINEQA,       "HGLN-SFL",  "HGLT-SFL", "Sanson-Flamsteed"},
    {GNOMONIC,     "CRLN-TAN",  "CRLT-TAN", "Gnomonic"},
    {GNOMONIC,     "HGLN-TAN",  "HGLT-TAN", "Gnomonic"},
    {POSTEL,       "CRLN-ARC",  "CRLT-ARC", "Postel"},
    {POSTEL,       "HGLN-ARC",  "HGLT-ARC", "Postel"},
    {STEREOGRAPHIC,"CRLN-STG",  "CRLT-STG", "Stereographic"},
    {STEREOGRAPHIC,"HGLN-STG",  "HGLT-STG", "Stereographic"},
    {ORTHOGRAPHIC, "CRLN-SIN",  "CRLT-SIN", "Orthographic"},
    {ORTHOGRAPHIC, "HGLN-SIN",  "HGLT-SIN", "Orthographic"},
    {LAMBERT,      "CRLN-ZEA",  "CRLT-ZEA", "Lambert"},
    {LAMBERT,      "HGLN-ZEA",  "HGLT-ZEA", "Lambert"},
    {AERIAL,       "HPLN-TAN",  "HPLT-TAN", "Aerial"},
    {-1, "END ", ""}
  };
  
  done = i = 0;
  while(strncasecmp(wcs[i].CTYPE1,"END ",4) && !done)  // First check by CTYPE
  {
    if(!strncasecmp(parms.CTYPE1,wcs[i].CTYPE1,8) &&
       !strncasecmp(parms.CTYPE2,wcs[i].CTYPE2,8))
    {
      parms.projcode = wcs[i].code;
      strcpy(parms.projname,wcs[i].projname);
//      printf("PROJ_CODE: %s %s - %s\n",parms.CTYPE1,parms.CTYPE2,parms.projname);
      done = 1;
    }
    i++;
  }
  i = 0;                              // If not done, check by projname
  while(strncasecmp(wcs[i].CTYPE1,"END ",4) && !done) 
  {
    if(!strncasecmp(parms.projname,wcs[i].projname,3))
    {
      parms.projcode = wcs[i].code;
      strcpy(parms.CTYPE1,wcs[i].CTYPE1);
      strcpy(parms.CTYPE2,wcs[i].CTYPE2);
      printf("PROJ_NAME: %s - %s %s\n",parms.projname,parms.CTYPE1,parms.CTYPE2);
      done = 1;
    }
    i++;
  }
  if(done)
    return kMyMod_Success;     // exit wcs_map_proj normally

  sprintf(parms.Msg,"%sError: Undetermined Map Projection- %s, %s, %s,\n",
	  parms.Msg, parms.CTYPE1, parms.CTYPE2, parms.projname);
  return (parms.status = kMyMod_ValErr);
}


// CALC_WEIGHT_FUNCTION
int calc_weight_function(void)
{
  int j, mid;
  double arg, Dmid, dt, Sum = 0.0;
  int npts = parms.WtLen;  
  char *filter = parms.WtFunc;
  double fwhm = parms.WtParm;
  double *Wt;
  
  if(( parms.Wt = (double *) malloc(sizeof(double)*npts))==NULL)
  {
    sprintf(parms.Msg,"%sError: Allocating Memory for weighting function",parms.Msg);
    return  (parms.status = kMyMod_MallocErr);  
  }
  Wt = parms.Wt;
  Dmid = 0.5*(npts - 1);
  mid = Dmid;
  
  if(!strncasecmp(filter,"box",3))
  {
    for(j = 0; j < npts; j++) 
      Wt[j] = 1.0;
  }
  else  if(!strncasecmp(filter,"tri",3))
  {
    dt = 1.0 / fwhm;
    for(j = 0; j <= mid; j++) 
    {
      arg = 1.0 - (Dmid - j) * dt;
      Wt[j] = Wt[npts-j-1] = (arg >= 0.0) ? arg: 0.0;
    }
  }
  else  if(!strncasecmp(filter,"sinc2",5))
  {
    dt = 2 * 0.442946 / fwhm;
    for(j = 0; j <= mid; j++) 
    {
      arg = (Dmid - j) * dt;
      if(arg != 0.0)
      {
	Wt[j] = Wt[npts - j-1] = sin(arg)/arg;
	Wt[npts-j-1] *= Wt[j];
	Wt[j] *= Wt[j];
      }
      else
	Wt[j] = Wt[npts - j-1] = 1.0;	    
    }
  }
  else  if(!strncasecmp(filter,"sinc",4))
  {
    dt = 2 * 0.603355 / fwhm;
    for(j = 0; j <= mid; j++) 
    {
      arg = (mid - j) * dt;
      if(arg != 0.0)
	Wt[j] = Wt[npts - j-1] = sin(arg)/arg;
      else
	Wt[j] = Wt[npts - j-1] = 1.0;	    
    }
  }
  else if(!strncasecmp(filter,"hath",4))
  {
    dt = exp(-2.0);
    for(j = 0; j < npts; j++) 
    {
      arg = (j - mid)/(fwhm);
      arg *= -0.5 * arg;
      Wt[j] = exp(arg) - dt*(3.0 + arg);
    }
  }
  else if(!strncasecmp(filter,"gauss",5))
  {
    dt = 2 * sqrt (log (2.0)) / fwhm;
    for(j = 0; j <= mid; j++) 
    {
      arg = (mid - j) * dt;
      arg *= -arg;
      Wt[j] = Wt[npts-j-1] = exp (arg);
    }
  }
  else
  {
    sprintf(parms.Msg,"%sUnknown Weight Function: %s\n",parms.Msg,filter);
    return (parms.status = kMyMod_WrongType);
  }
/*
  Sum = 0.0;
  for(j = 0; j < npts; j++)
    Sum += Wt[j];

  if (Sum == 0.0)
  {
    sprintf(parms.Msg,"%sError Calculating Weight Function\n",parms.Msg);
    return (parms.status = kMyMod_InitErr);
  }
  for(j = 0; j < npts; j++)
    Wt[j] /= Sum;
*/
  return kMyMod_Success;     // exit calc_weight_func normally
}


/* 
******************************************************************** 
* 
*  Functions involving setting up output arrays and checking outDB
* 
******************************************************************** 
*/
// SETUP_OUTDB
int  setup_outdb(OutImgs **outimg, DRMS_RecordSet_t **outdb)
{
  int status;

  strcpy(parms.SigMsg," ");
  strcpy(parms.Msg," ");
  (*outdb) = drms_create_records(drms_env, parms.Nmaps,
				 parms.oser, DRMS_PERMANENT,&status);

  if(!(*outdb) || (*outdb)->n != parms.Nmaps)
  {
    sprintf(parms.Msg,"%sError creating records in series %s\n",
	    parms.Msg,parms.oser);
    return (parms.status = kMyMod_Missing);
  }
  
  if((status = check_outdb((*outdb)->records[0])))
    return status;

  if((status = init_outimgs(outimg)))
    return status;
//   drms_close_records((*outdb), DRMS_FREE_RECORD);
  V_print("OK\n");
  return  kMyMod_Success;     // exit setup_outdb normally
}

// CHECK_OUTDB

int    check_outdb(DRMS_Record_t *rec)
{
  int i, done, segcnt;
  int status = 0;
  char buff[200];
  TIME ttime;
  DRMS_Segment_t *record_segment;
  DRMS_Keyword_t *keywd;
  static struct drmean_keys 
  {
    int type;
    char name[20];
  } keys[]= {
    {DRMS_TYPE_TIME,    "TSTART"},
    {DRMS_TYPE_TIME,    "TEND"},
    {DRMS_TYPE_TIME,    "T_REC_STEP"},
    {DRMS_TYPE_TIME,    "T_REC"},
    {DRMS_TYPE_DOUBLE,  "A0"},
    {DRMS_TYPE_DOUBLE,  "A2"},         
    {DRMS_TYPE_DOUBLE,  "A4"},         
    {DRMS_TYPE_DOUBLE,  "MERI_V"},     
    {DRMS_TYPE_STRING,  "WTFUNC"}, 
    {DRMS_TYPE_DOUBLE,  "WTPARM"},     
    {DRMS_TYPE_INT,     "WTLEN"},      
    {DRMS_TYPE_STRING,  "PROJNAME"},  
    {DRMS_TYPE_DOUBLE,  "RSUN"},       
    {DRMS_TYPE_DOUBLE,  "DIST"},       
    {DRMS_TYPE_DOUBLE,  "ANG_RAD"},    
    {DRMS_TYPE_DOUBLE,  "PANGLE"},     
    {DRMS_TYPE_DOUBLE,  "CROTA2"},      
    {DRMS_TYPE_DOUBLE,  "CRLN"},       
    {DRMS_TYPE_DOUBLE,  "CRLT"},       
    {DRMS_TYPE_STRING,  "CTYPE1"},     
    {DRMS_TYPE_STRING,  "CTYPE2"},     
    {DRMS_TYPE_DOUBLE,  "CRPIX1"},     
    {DRMS_TYPE_DOUBLE,  "CRPIX2"},     
    {DRMS_TYPE_DOUBLE,  "CRVAL1"},     
    {DRMS_TYPE_DOUBLE,  "CRVAL2"},     
    {DRMS_TYPE_DOUBLE,  "CDELT1"},     
    {DRMS_TYPE_DOUBLE,  "CDELT2"},   
    {DRMS_TYPE_RAW,     "END "}
  };
  
  // First check to ensure that there are two segments (one for data
  // and the other for weights and that the dimensions are consistent
  
  if((segcnt = drms_record_numsegments (rec)) != 2)
    {
      sprintf(parms.Msg,"%sError: Bad Output DB - need 2 segments\n",parms.Msg);
      return (parms.status = kMyMod_Missing);
    }
  for(i = 0; i < segcnt; i++) 
  {
    record_segment = drms_segment_lookupnum (rec, i);
    if(record_segment->info->naxis != 2) 
    {
      sprintf(parms.Msg,"%sError: Bad Output DB Seg %d - not 2D\n",parms.Msg,i);
      return (parms.status = kMyMod_Missing);
    }
    if(record_segment->axis[0] != parms.cols ||
	record_segment->axis[1] != parms.rows ) 
    {
      parms.cols = record_segment->axis[0];
      parms.rows = record_segment->axis[1];
      parms.size = parms.rows * parms.cols;
      sprintf(buff,"Specified dimensions do not match segment definition. Set to %d,%d\n",
	      parms.cols, parms.rows);
      V_print(buff);
    }
    if( isnan(parms.CRPIX1))
      parms.CRPIX1 = 0.5*parms.cols;
    if( isnan(parms.CRPIX2))
      parms.CRPIX2 = 0.5*parms.rows;

  }
  
  done = i = 0;
  while(strcasecmp(keys[i].name,"END ") && !done)
  {
    if(!(keywd = drms_keyword_lookup (rec, keys[i].name, 0))) 
    {
      sprintf(parms.Msg,"%sError: Bad Output DB - missing keyword %s\n",
	      parms.Msg, keys[i].name);
      return (parms.status = kMyMod_Missing);
    }
    if(keywd->info->type != keys[i].type)
    {
      sprintf(parms.Msg,"%sError: Bad Output DB - keyword %s is wrong type\n",
	      parms.Msg,keys[i].name);
      return (parms.status = kMyMod_WrongType);
    }
    i++;
  }
  
  // Check that T_REC_step matches step_t variable
  ttime = drms_getkey_time(rec, "T_REC_step", &status);
  if(status || ttime > parms.Tstep)
  {
    sprintf(parms.Msg,"%sError: STEP_T is less than T_REC_step\n",parms.Msg);
    return (parms.status = kMyMod_WrongType);
  }
  
  return  kMyMod_Success;   // exit  check_outdb normally
}
  
// INIT_OUTIMGS
int    init_outimgs(OutImgs **outimg)
{
  int i, status;
  char buff[200], tbuff[200];
  int size = parms.size;
  double Bz;
  OutImgs *out;
  int naxis = 2;
  int naxes[2] = {parms.cols, parms.rows}; 

  
  out = (OutImgs *)malloc(sizeof(OutImgs)*parms.Nmaps);
  status = 0;
  sprintf(buff,"%4dx%4d %s,%s CDELT=[%8.5f,%8.5f] CRPIX=[%8.2f,%8.2f]\n", 
	  parms.cols,parms.rows, 
	  parms.CTYPE1, parms.CTYPE2,
	  parms.CDELT1, parms.CDELT2,
	  parms.CRPIX1, parms.CRPIX2
    );
  V_print(buff);
  V_print("-----------------------\n");

  for(i =0; i < parms.Nmaps && !status; i++)
  {
    out[i].DatArray = drms_array_create(DRMS_TYPE_FLOAT, naxis, naxes, NULL, &status);
    out[i].DatArray->israw = 0;
    if(status)
    {
      sprintf(parms.Msg,"%sError: Creating output arrays",parms.Msg);
      return (parms.status = status);
    }
    out[i].dat = (float *)out[i].DatArray->data;

    out[i].WgtArray = drms_array_create(DRMS_TYPE_FLOAT, naxis, naxes, NULL, &status);
//    out[i].WgtArray->israw = 0;
    if(status)
    {
      sprintf(parms.Msg,"%sError: Creating output arrays",parms.Msg);
      return (parms.status = status);
    }
    out[i].wgt = (float *)out[i].WgtArray->data;
    
    if((out[i].lon = (float *)malloc(sizeof(float)*size))==NULL)
      status = kMyMod_MallocErr;
    if((out[i].lat = (float *)malloc(sizeof(float)*size))==NULL)
      status = kMyMod_MallocErr;
    if((out[i].osun = (short *)malloc(sizeof(short)*size))==NULL)
      status = kMyMod_MallocErr;

    if(status)
    {
      sprintf(parms.Msg,"%sError: Allocating Memory for output arrays",
	      parms.Msg);
      return (parms.status = status);
    }
    out[i].timctr = parms.Tstart + i * parms.Tstep;
    out[i].LN = parms.CRVAL1 - (i * parms.Tstep * parms.A0);
    out[i].LT = parms.CRVAL2 + i * parms.Tstep * parms.Meri_V;
    HeliographicLocation(out[i].timctr, &(out[i].crot), &(out[i].crln), &Bz);
    out[i].irecdt = parms.WtLen;
    out[i].irecno = -1;
    sprint_time(tbuff,out[i].timctr,"TAI",0);
    sprintf(buff,"[%d] %27s[%.19s]; LN,LT=[%8.4f,%8.4f]\n",i, 
	    parms.oser, tbuff, out[i].LN,out[i].LT);
    V_print(buff);
    
    if((status = calc_out_geom(out[i])))
    {
      sprintf(parms.Msg,"%sError: Calculating Geometry for output image",parms.Msg);
      return (parms.status = status);
    }
  }
  (*outimg) = out; 
  return kMyMod_Success;     // exit init_outimgs normally
}


// CALC_OUT_GEOM
int calc_out_geom(OutImgs outimg)
{
  double x, y, xrot, yrot, lon, lat;
  int row, col, i, status;
  static int first = 1;
  static double sin_phi, cos_phi, x0, y0;

  if(first)
  {
    first = 0;
    cos_phi = cos (parms.pangle);
    sin_phi = sin (parms.pangle);
    y0 = parms.CRPIX2;
    x0 = parms.CRPIX1;
  }
  
  status = i = 0;
  for(row=0; row < parms.rows; row++)
  {
    y = (row-y0)*parms.CDELT2;
    for(col = 0; col < parms.cols; col++ ) 
    {
      i = col + row * parms.cols;
      x = (col-x0)*parms.CDELT1;
      xrot = (x * cos_phi - y * sin_phi);
      yrot = (y * cos_phi + x * sin_phi);
      if(parms.projcode == AERIAL)
      {  	// AERIAL: call img2sphere instead of plane2sphere.  
	double RHO,SLT,CLT,SIG,MU,CHI;  // Dummy Vars for passing 
	outimg.osun[i] = img2sphere(xrot,yrot,parms.ang_rad,outimg.LT,
				    outimg.LN,0.0, &RHO,&lat,&lon,&SLT,&CLT,
				    &SIG,&MU,&CHI);
      } 
      else 
      {
	outimg.osun[i] = plane2sphere(xrot,yrot,outimg.LT,outimg.LN,&lat,&lon,
				      parms.projcode);
      }
      if(isnan(lat) || isnan(lon))
	outimg.osun[i] = kMyMod_ValErr;
      outimg.dat[i] = 0.0;
      outimg.wgt[i] = 0.0;
      outimg.lat[i] = lat;
      outimg.lon[i++] = lon;
    } 
  }
  return kMyMod_Success;     // exit geom_outimgs normally
}



/* 
 ******************************************************************* 
 * 
 *  Functions involving opening and  checking input DB
 * 
 ******************************************************************** 
*/
// OPEN_CHECK_INDB
int  open_check_indb(DRMS_RecordSet_t **indb)
{
  int status;
  char temp[2000];
  char qry_str[DRMS_MAXQUERYLEN];
  DRMS_Record_t *iRec;
  
  strcpy(parms.SigMsg," ");
  strcpy(parms.Msg," ");
  /*  Need to check that query returned actual records */
  if(get_query_str(qry_str) != kMyMod_Success)
    return (parms.status = kMyMod_Missing);

  sprintf(temp,"%s\n",qry_str);
  V_print(temp);

  (*indb) = drms_open_records(drms_env, qry_str, &status);

  if(status)
  {
    sprintf(parms.Msg,"%sError %d: in opening input data records\n",
	    parms.Msg, status);
    return (parms.status = status);
  }
  if((*indb)->n == 0)
  {
    sprintf(parms.Msg,"%sNo input data records -Check date/time\n",parms.Msg);
    return (parms.status = kMyMod_Missing);
  }
  iRec = (*indb)->records[0];
  if(iRec == NULL) 
    printf("iRec is Null\n");
  //  Now, run through and check the existence of required keywords
  if((status = check_indb(iRec)) != kMyMod_Success)
     return status;
     
  V_print("OK\n");
  return kMyMod_Success;     // exit open_check_indb normally
}

// GET_QUERY_STR
int    get_query_str(char *qry)
{
  const char *test;
  char *p;
  char tbuf1[100],tbuf2[100];

  test = cmdparams_get_str(&cmdparams, "in", NULL);
  if (strchr(test,'[') != NULL)
    strcpy(qry,test);
  else
  {
    sprint_time(tbuf1, parms.Tstart - (parms.WtLen)/2 - 60, "TAI", 0);
    sprint_time(tbuf2, parms.Tend + (parms.WtLen)/2 + 60, "TAI", 0);
    sprintf(qry,"%s[%s - %s]",test,tbuf1,tbuf2);
  }
  return kMyMod_Success;     // exit open_check_indb normally
}

// CHECK_INDB    checks the structure of input data series
int    check_indb(DRMS_Record_t *rec)
{
  int i, done, segcnt;
  int status = 0;
  DRMS_Keyword_t *keywd;
  static struct drmean_keys
  {
    int type;
    char name[20];
  } keys[]= {
//    {DRMS_TYPE_DOUBLE,  "R_SUN"},       
    {DRMS_TYPE_DOUBLE,  "DSUN_OBS"},       
//    {DRMS_TYPE_DOUBLE,  "ANG_RAD"},    
//    {DRMS_TYPE_DOUBLE,  "PANGLE"},     
    {DRMS_TYPE_DOUBLE,  "CROTA2"},      
    {DRMS_TYPE_DOUBLE,  "CRLN_OBS"},       
    {DRMS_TYPE_DOUBLE,  "CRLT_OBS"},       
    {DRMS_TYPE_STRING,  "CTYPE1"},     
    {DRMS_TYPE_STRING,  "CTYPE2"},     
    {DRMS_TYPE_DOUBLE,  "CRPIX1"},     
    {DRMS_TYPE_DOUBLE,  "CRPIX2"},     
    {DRMS_TYPE_DOUBLE,  "CRVAL1"},     
    {DRMS_TYPE_DOUBLE,  "CRVAL2"},     
    {DRMS_TYPE_DOUBLE,  "CDELT1"},     
    {DRMS_TYPE_DOUBLE,  "CDELT2"},   
    {DRMS_TYPE_RAW,     "END "}
  };
  
  done = i = 0;
  while(strcasecmp(keys[i].name,"END "))
  {
    if(!(keywd = drms_keyword_lookup (rec, keys[i].name, 0))) 
    {
      sprintf(parms.Msg,"%sError: Bad Input DB - missing keyword: %s\n",
	      parms.Msg, keys[i].name);
      return (parms.status = kMyMod_Missing);
    }
//    if(keywd->info->type != keys[i].type)
//   {
//      sprintf(parms.Msg,"%sError: Bad Input DB - keyword %s is wrong type\n",
//	      parms.Msg, keys[i].name);
//      return (parms.status = kMyMod_WrongType);
//    }
    i++;
  }

  if(!(keywd = drms_keyword_lookup (rec, "R_SUN", 0))) 
  {
    if(!(keywd = drms_keyword_lookup (rec, "RSUN_OBS", 0))) 
    {
      sprintf(parms.Msg,"%sError: Input DB - missing R_SUN && RSUN_OBS one is required\n",
	      parms.Msg);
      return (parms.status = kMyMod_Missing);
    }
  }
  
  return kMyMod_Success;     // exit check_indb normally
}



/* 
 ******************************************************************* 
 * 
 *  Derotation & summation  Functions 
 * 
 ******************************************************************** 
*/
// ADDIN

int    addin(DRMS_RecordSet_t *inDB, OutImgs *outimgs)
{
  int i, nInRecs, status;
  int dt, cols, rows, m, n;
  int WtLen2 = parms.WtLen / 2;
  char repstr[DRMS_MAXQUERYLEN], temp[1000], ctype1[20], ctype2[20];
  float *in_data;
  double Bo, DT, Lo, SB2, lat, lon, phi, Rsun, xc, xx, yc, yy, TEST; 
  TIME timg;
  DRMS_Array_t *data_array = NULL;
  DRMS_Segment_t *dseg;
  DRMS_Keyword_t *keywd;
  DRMS_Record_t    *inRec = NULL;
  
  strcpy(parms.SigMsg," ");
  strcpy(parms.Msg," ");

  nInRecs = inDB->n;
  for(i =0; i < nInRecs; i++)
  {
    inRec = inDB->records[i];
    drms_sprint_rec_query(repstr,inRec);
    if((readin(inRec, &timg, &Bo, &Lo, &Rsun, &phi, 
	       &xc, &yc, ctype1, ctype2) == kMyMod_Success))
    {
      sprintf(temp,"%5.1f %5.1f %5.1f %5.1f %5.1f ",Lo,Bo,xc,yc,Rsun);
      strcat(repstr,temp);
      dseg = drms_segment_lookupnum (inRec, 0);
      if(data_array != NULL)
      {
	drms_free_array(data_array);
	data_array = NULL;
      }
      data_array = drms_segment_read (dseg, DRMS_TYPE_FLOAT, &status);
      in_data = (float *)data_array->data;
      cols = data_array->axis[0];
      rows = data_array->axis[1];
      for(m = 0; m < parms.Nmaps; m++) 
      {
	dt = DT = timg - outimgs[m].timctr;
	if(abs(DT) <= WtLen2)
	{
	  sprintf(temp,"+%d,",m);
	  strcat(repstr,temp);
	  if(abs(DT) < abs(outimgs[m].irecdt))
	  {
	    outimgs[m].irecdt = DT;
	    outimgs[m].irecno = i;
	  }
	  for(n = 0; n < parms.size; n++) 
	  {
	    if(!outimgs[m].osun[n]) 		
	    {
	      lat = outimgs[m].lat[n] + parms.Meri_V * DT;
	      lon = outimgs[m].lon[n] + parms.A0 * DT;
	      if(parms.A2 != 0 || parms.A4 != 0)
	      {
		SB2 = sin(lat) * sin(lat);
		lon -= (SB2 * parms.A2 + parms.A4 * SB2 * SB2) * DT;
	      }
	      sphere2img (lat, lon, Bo, Lo, &xx, &yy, xc, yc, 
			  Rsun, phi, 0.0, 1.0, 0.0, 0.0); 
	      if(!isnan(TEST = ccint2(in_data,cols,rows,xx,yy)))
	      {
		outimgs[m].dat[n] += TEST * parms.Wt[dt + WtLen2];
		outimgs[m].wgt[n] += parms.Wt[dt + WtLen2];
	      }
	    }
	  }
	}
      }
    }
    strcat(repstr,"\n");
    V_print(repstr);
  }

  for(m = 0; m < parms.Nmaps; m++) 
  {
    for(n = 0; n < parms.size; n++) 
    {
      if(outimgs[m].wgt[n] > 0.0)
	outimgs[m].dat[n] = outimgs[m].dat[n] / outimgs[m].wgt[n];
      else
	outimgs[m].dat[n] = 0.0/0.0;
    }
  }    
  return kMyMod_Success;     // exit addin normally
}


// READIN reads input image parameters

int    readin(DRMS_Record_t *irec, TIME *timg, double *Bo, double *Lo, 
	      double *Rsun, double *phi,
              double *xc, double *yc, char ctype1[], char ctype2[])
{
  int status;
  double dsun, cdelt, rang;
  double R2D = 180.0 / M_PI;
//  char temp1[100], tbuff[100];
  
  if(irec == NULL)
    return kMyMod_Missing;

  if((status = rec_getcheck_time(irec, "T_REC",timg)))
    if((status = rec_getcheck_time(irec, "T_OBS",timg)))
      return status;
  
  if(status = rec_getcheck_str(irec,"CTYPE1",ctype1))
    return status;
  if(status = rec_getcheck_str(irec,"CTYPE2",ctype2))
    return status;
  if(status = rec_getcheck_double(irec,"CRVAL1",Lo))
    return status;
  if(status = rec_getcheck_double(irec,"CRVAL2",Bo))
    return status;
  if(status = rec_getcheck_double(irec,"CRPIX1",xc))
    return status;
  if(status = rec_getcheck_double(irec,"CRPIX2",yc))
    return status;
  if(status = rec_getcheck_double(irec,"CROTA2",phi) ||
     drms_ismissing_double(*phi))
    (*phi) = 0.0;
  (*phi) = -1*(*phi) / R2D;    // Convert from degrees to Radians

  if(
     drms_ismissing_double(*Lo) || drms_ismissing_double(*Bo) ||
     drms_ismissing_double(*xc) || drms_ismissing_double(*yc) ||
     drms_ismissing_string(ctype1) || drms_ismissing_string(ctype2) ||
     drms_ismissing_time(*timg) )
    return kMyMod_Missing;

/* To get the Solar Radius - in pixels instead read DSUN_OBS and
    convert.
  if(status = rec_getcheck_double(irec,"RSUN_OBS",Rsun))
    if(status = rec_getcheck_double(irec,"R_SUN",Rsun))
*/
  if(status = rec_getcheck_double(irec,"DSUN_OBS",&dsun))
    return status;
  rang = 3600.*57.295778 * asin(6.96e8/dsun);
  if(status = rec_getcheck_double(irec,"CDELT1",&cdelt))
    return status;
  (*Rsun) = rang / cdelt;
  
  return kMyMod_Success;     // exit readin normally
}


// REC_GETCHECK_DOUBLE
int rec_getcheck_double(DRMS_Record_t *rec, char kword[], double *out)
{
  int stat;
  char ctmp[1000];
  double temp;

  temp = drms_getkey_double(rec, kword,&stat);
  if(stat)
  {  
    drms_sprint_rec_query(ctmp,rec);
    sprintf(ctmp,"%s:  MISSING KEYWORD %s ",ctmp, kword);
    V_print(ctmp);
    return stat;
  }
  (*out) = temp;
  return kMyMod_Success;     // exit rec_getcheck_double normally
  
}

// REC_GETCHECK_INT
int rec_getcheck_int(DRMS_Record_t *rec, char kword[], int *out)
{
  int stat;
  char ctmp[1000];
  int temp;

  temp = drms_getkey_int(rec, kword,&stat);
  if(stat)
  {
    drms_sprint_rec_query(ctmp,rec);
    sprintf(ctmp,"%s:  MISSING KEYWORD %s ",ctmp, kword);
    V_print(ctmp);
    return stat;
  }
  (*out) = temp;
  return kMyMod_Success;     // exit rec_getcheck_int normally
  
}
// REC_GETCHECK_STR
int rec_getcheck_str(DRMS_Record_t *rec, char kword[], char out[])
{
  int stat;
  char ctmp[1000];
  const char *temp;
  temp = drms_getkey_string(rec,kword,&stat);
  if(stat)
  {
    drms_sprint_rec_query(ctmp,rec);
    sprintf(ctmp,"%s:  MISSING KEYWORD %s ",ctmp, kword);
    V_print(ctmp);
    return stat;
  }
  strcpy(out,temp);
  return  kMyMod_Success;     // exit rec_getcheck_str normally
  
}

// REC_GETCHECK_TIME
int rec_getcheck_time(DRMS_Record_t *rec, char kword[], TIME *out)
{
  int stat;
  char ctmp[1000];
  TIME temp;

  temp = drms_getkey_time(rec, kword,&stat);
  if(stat)
  {
    drms_sprint_rec_query(ctmp,rec);
    sprintf(ctmp,"%s:  MISSING KEYWORD %s ",ctmp, kword);
    V_print(ctmp);
    return stat;
  }
  (*out) = temp;
  return kMyMod_Success;     // exit rec_getcheck_time normally
  
}

/* 
 ******************************************************************* 
 * 
 *  Checking & Writing Out Img Functions
 * 
 ******************************************************************** 
*/
// WRITE_OUT_DBS
int    write_out_dbs(OutImgs *outimgs, DRMS_RecordSet_t *outDB, DRMS_RecordSet_t *inDB)
{
  int i, nrecs, status;
  char cunit[20];
  char buff[200];
  DRMS_Record_t *orec, *irec;
  DRMS_Segment_t *oseg;
  double SCL, R2D, R2AS;
//  DRMS_RecordSet_t *outdb;

  strcpy(parms.SigMsg," ");
  strcpy(parms.Msg," ");

//  outdb = drms_create_records(drms_env, parms.Nmaps,
//				 parms.oser, DRMS_PERMANENT,&status);
  R2D = 180.0 / M_PI;
  R2AS = 3600.0 *180.0 / M_PI;

// Have to re-do SCL factor to be consistent with input,
// comp with check_cparms... the key CDELT2 is off!!

  nrecs = parms.Nmaps;
  if(parms.projcode == AERIAL)
  {
    strcpy(cunit,"arcsec");
    SCL = R2D * 3600.00;
  }
  else
  {
     strcpy(cunit,"deg");
     SCL = R2D;
  }

  for(i = 0; i < nrecs; i++)
  {
    orec = outDB->records[i];
//    orec = drms_create_record(drms_env,parms.oser,DRMS_PERMANENT, &status);
    irec = inDB->records[outimgs[i].irecno];
    drms_copykeys(orec,irec,0, kDRMS_KeyClass_Explicit ); 
    drms_setkey_time(orec,"T_REC",outimgs[i].timctr);
    drms_setkey_double(orec,"A0",parms.A0*1.0e6);
    drms_setkey_double(orec,"A2",parms.A2*1.0e6);
    drms_setkey_double(orec,"A4",parms.A4*1.0e6);
    drms_setkey_double(orec,"MERI_V",parms.Meri_V);
    drms_setkey_string(orec,"WTFUNC",parms.WtFunc);
    drms_setkey_double(orec,"WTPARM",parms.WtParm);
    drms_setkey_int(orec,"WTLEN",parms.WtLen);
    drms_setkey_string(orec,"PROJNAME",parms.projname);
    drms_setkey_double(orec,"RSUN",parms.rsun);
    drms_setkey_double(orec,"DIST",parms.dist);
    drms_setkey_double(orec,"ANG_RAD",parms.ang_rad*R2AS);
    drms_setkey_double(orec,"PANGLE_RAD",parms.pangle*R2D);
    drms_setkey_double(orec,"CROTA2",parms.CROTA2*R2D);
    drms_setkey_double(orec,"CRLN",outimgs[i].LN*SCL);
    drms_setkey_double(orec,"CRLT",outimgs[i].LT*SCL);
    drms_setkey_string(orec,"CTYPE1",parms.CTYPE1);
    drms_setkey_string(orec,"CTYPE2",parms.CTYPE2);
    drms_setkey_double(orec,"CRPIX1",parms.CRPIX1);
    drms_setkey_double(orec,"CRPIX2",parms.CRPIX2);
    drms_setkey_double(orec,"CRVAL1",outimgs[i].LN*SCL);
    drms_setkey_double(orec,"CRVAL2",outimgs[i].LT*SCL);
    drms_setkey_double(orec,"CDELT1",parms.CDELT1*SCL);
    drms_setkey_double(orec,"CDELT2",parms.CDELT2*SCL);
    drms_setkey_string(orec,"CUNIT1",cunit);
    drms_setkey_string(orec,"CUNIT2",cunit);
    drms_setkey_double(orec,"CADENCE",parms.Tstep);
    drms_setkey_double(orec,"CRLN_OBS",outimgs[i].crln);
    drms_setkey_int(orec,"CAR_ROT",outimgs[i].crot);

    if((parms.status = set_bscale_bzero(&outimgs[i])) != kMyMod_Success)
    {
      sprintf(parms.Msg,"%sError: setting BSCALE & BZERO\n",parms.Msg);
      return (parms.status = kMyMod_InitErr);
    } 
    sprintf(buff,"[%02d] Data:   Max=%9.4f, Min=%9.4f, Bscale=%9.4f, Bzero=%9.4f\n",
	    i, outimgs[i].Dmax, outimgs[i].Dmin, 
	    outimgs[i].DatArray->bscale, outimgs[i].DatArray->bzero);
    V_print(buff);
    sprintf(buff,"[%02d] Weight: Max=%9f, Min=%9f, Bscale=%9f, Bzero=%9f\n",
	    i, outimgs[i].Wmax, outimgs[i].Wmin, 
	    outimgs[i].WgtArray->bscale, outimgs[i].WgtArray->bzero);
    V_print(buff);

    oseg = drms_segment_lookupnum (orec, 0);
    oseg->bzero = outimgs[i].DatArray->bzero;
    oseg->bscale = outimgs[i].DatArray->bscale;
    outimgs[i].DatArray->parent_segment = oseg;
    outimgs[i].DatArray->israw = 0;
    drms_segment_write(oseg,outimgs[i].DatArray,0);
    drms_free_array(outimgs[i].DatArray);

    oseg = drms_segment_lookupnum (orec, 1);
    oseg->bzero = outimgs[i].WgtArray->bzero;
    oseg->bscale = outimgs[i].WgtArray->bscale;
    outimgs[i].WgtArray->parent_segment = oseg;
    outimgs[i].WgtArray->israw = 0;
    drms_segment_write(oseg,outimgs[i].WgtArray,0);
    drms_free_array(outimgs[i].WgtArray);

  }

  drms_close_records(outDB, DRMS_INSERT_RECORD);
  status = free_outimgs(outimgs);

  V_print("OK\n");
  return kMyMod_Success;     // exit write_out_dbs normally
}

// SET_BSCALE_BZERO
int    set_bscale_bzero(OutImgs *outimg)
{
  int i;
  double TAPERNG = 3.0e4;
  char buff[200];

  outimg->DatArray->bscale = 1.0;
  outimg->DatArray->bzero = 0.0;
  outimg->WgtArray->bscale = 0.001;
  outimg->WgtArray->bzero = 0.0;
  outimg->Dmax = SHRT_MIN;
  outimg->Dmin = SHRT_MAX;
  outimg->Wmax = SHRT_MIN;
  outimg->Wmin = SHRT_MAX;

  for (i = 0; i < parms.size; i++)
  {
    if(!isnan(outimg->dat[i]))
    {
      if(outimg->Dmax < outimg->dat[i])
	outimg->Dmax = outimg->dat[i];
      if(outimg->Dmin > outimg->dat[i])
	outimg->Dmin = outimg->dat[i];
    }
    if(!isnan(outimg->wgt[i]))
    {
      if(outimg->Wmax < outimg->wgt[i])
	outimg->Wmax = outimg->wgt[i];
      if(outimg->Wmin > outimg->wgt[i])
	outimg->Wmin = outimg->wgt[i];
    }
  }

  if (outimg->Wmin < 0.0)
  { 
    outimg->WgtArray->bscale = (outimg->Wmax - outimg->Wmin)/TAPERNG;
    outimg->WgtArray->bzero = outimg->Wmin;
  }
  else
  {
    outimg->WgtArray->bscale = (outimg->Wmax)/TAPERNG;
    outimg->WgtArray->bzero = 0.0;
  }
  
  if(outimg->Dmax >= TAPERNG || outimg->Dmin <= -1*TAPERNG || 
     abs(outimg->Dmax - outimg->Dmin) < 0.01*TAPERNG)
  {
    outimg->DatArray->bscale = (outimg->Dmax - outimg->Dmin)/TAPERNG/2.0;
    outimg->DatArray->bzero = 0.5*(outimg->Dmax + outimg->Dmin);
    if(abs(outimg->DatArray->bzero) < 0.08*(outimg->DatArray->bscale)*TAPERNG)
      outimg->DatArray->bzero = 0.0;
  }

  return kMyMod_Success;   // Exit set_bscale_bzero normally - setvals
}



/* 
 ******************************************************************* 
 * 
 *  Misc Utility Functions 
 * 
 ******************************************************************** 
*/


// WCS_PARM_UNIT
double wcs_parm_unit(double X, char *unit)
{
  static double d2r =  M_PI / 180.0;
  static double m2r =  M_PI / 180.0 / 60.0;
  static double a2r =  M_PI / 180.0 / 3600.0;
  
  if(!strcmp(unit,"arcsec")) 
    return X * a2r;
  else if(!strcmp(unit,"arcmin"))
    return X * m2r;
  else if(!strcmp(unit,"deg")) 
    return X * d2r;
  else if(!strcmp(unit,"rad")) 
    return X;
  
  return 0.0/0.0;
}

// FREE_OUTIMGS
int    free_outimgs(OutImgs *outimg)
{
  int i;

  for(i =0; i < parms.Nmaps; i++)
  {
    free(outimg[i].lon);
    free(outimg[i].lat);
    free(outimg[i].osun);
  }
  free(outimg);
  return kMyMod_Success;     // exit free_outimgs normally
}


// IF_ERR_PRINT
int If_Err_Print(void)
{
  if(parms.status)
  {
    fprintf(stderr, "Error %d: %s\n",parms.status,parms.Msg);
    fflush(stderr);
  }
  return parms.status;
}

//   V_PRINT
void   V_print(char *str)
{
  if (!parms.verbose && !(parms.verbose = cmdparams_isflagset(&cmdparams,"v")))
    return;
  printf("%s",str);
  fflush(stdout);
  
}



/*
 *  Revision History
 *  v 0.0  09.12.07     John Beck      created this file
 *
 */


