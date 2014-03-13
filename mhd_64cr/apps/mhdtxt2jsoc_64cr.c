/* from ingest_from_fits.c */

#include "jsoc_main.h"
#include "drms_types.h"

#if 0
#include "cfitsio.h"
#include "fitsio.h"
#endif

/* 1 = moderate resolution with definitive CR map, others=NRT low-resol.*/
#define NONDAILY 1
/* 1 = new map series having two segment, or not. */
#define NEWMAP 1

int set_statistics(DRMS_Segment_t *, DRMS_Array_t *, int);

char *module_name = "mhdtxt2jsoc";
#define DIE(msg) {fflush(stdout);fprintf(stderr,"%s, status=%d\n",msg,status); return(status);}
ModuleArgs_t module_args[] =
{
#if NONDAILY == 1
     {ARG_STRING, "in", "d3equidist144x72.txt", "Input text file."},
     {ARG_STRING, "out", "hmi.MHDcorona", "Target DRMS data series."},
#if NEWMAP == 1
     {ARG_STRING, "synomap", "hmi.BrBlossynop_720s",  "synoptic map series as string"},
#else
     {ARG_STRING, "synomap", "hmi.Synoptic_Mr_720s",  "synoptic map series as string"},
#endif
     {ARG_STRING, "crdate", "-1",  "date or CR as string"},
     {ARG_INT,    "sphindx", "-1",  "terms of pfss spherical harmonic poly. (negative for raw map)"},
     {ARG_INT,    "inx", "144",  "long. grid size of input MHD data, in the text file"},
     {ARG_INT,    "iny", "72",  "lat.  grid size of input MHD data, in the text file"},
     {ARG_INT,    "inz", "80",  "rad.  grid size of input MHD data, in the text file"},
#else
     {ARG_STRING, "in", "d3equidist.txt", "Input text file."},
     {ARG_STRING, "out", "hmi.MHDcorona_daily_nrt", "Target DRMS data series."},
#if NEWMAP == 1
     {ARG_STRING, "synomap", "hmi.BrBlosdailysynframe_720s_nrt", "synoptic map series as string"},
#else
     {ARG_STRING, "synomap", "hmi.Mrdailysynframe_720s_nrt",  "synoptic map series as string"},
#endif
     {ARG_STRING, "crdate", "$",  "date as string"},
     {ARG_INT,    "sphindx", "5",  "terms of pfss spherical harmonic poly. (negative for raw map)"},
     {ARG_INT,    "inx", "72",  "long. grid size of input MHD data, in the text file"},
     {ARG_INT,    "iny", "36",  "lat.  grid size of input MHD data, in the text file"},
     {ARG_INT,    "inz", "80",  "rad.  grid size of input MHD data, in the text file"},
#endif
     {ARG_INT,    "ichrbnd", "0",  "index for choice of sub-alfvenic boundary treatement"},
     {ARG_END}
};

int DoIt(void)
  {
  CmdParams_t *params    = &cmdparams;
  const char *instrc     = params_get_str(params, "in");
  const char *outstrc    = params_get_str(params, "out");
  const char *crdatestrc = params_get_str(params, "crdate");
  const char *synostrc   = params_get_str(params, "synomap");
  const int   sphindxc = params_get_int(params, "sphindx");
  const int   inxc     = params_get_int(params, "inx");
  const int   inyc     = params_get_int(params, "iny");
  const int   inzc     = params_get_int(params, "inz");
  const int   ichrbndc = params_get_int(params, "ichrbnd");
/* copy from runtime-option constant to variable */
  char   *instr, *outstr, *crdatestr, *synostr;
  instr     = strdup(instrc);
  outstr    = strdup(outstrc);
  crdatestr = strdup(crdatestrc);
  synostr   = strdup(synostrc);
  int     isphindx = sphindxc;
  int     nx = inxc;
  int     ny = inyc;
  int     nz = inzc;
  int     ichrbnd = ichrbndc;

/* JSOC-DRMS variables */
  DRMS_RecordSet_t *inRS, *outRS;
  DRMS_Record_t    *outRec, *inRec;
  DRMS_Segment_t   *outSeg, *inSeg;
  DRMS_Array_t     *inArray, *outArray;

/* some working variables */
  int  status;
  int  i, j, k, n;
  double ddummy;
  char *strdummy;
/* version info. in CVS tree */
  char *mhdcorona_version();

#if NONDAILY == 1
  int ncr = atoi(crdatestr);
  if (ncr < 0)
  {
    printf("CR is not given, so quit. %d %s\n",ncr,crdatestr);
    return (-1);
  }
#endif

/* MHD variable arrays */
  printf("Nx Ny Nz = %d %d %d\n",nx,ny,nz);
  float *mhd8all;
  mhd8all=(float *)malloc(sizeof(float)*(nx*ny*nz*8));

/* read plain text file storing equi-radius,lat,lon 3D-MHD variables */
  int icount=0;
{
  FILE *fp=NULL;  
  fp=fopen(instr,"r");
  char line[256];
  if(fp == NULL)
  {
    printf("Fail to open file %s\n",instr);
    return 1;
  }
  while(fgets(line,256,fp) != NULL)
  {
    float rdummy;
    sscanf(line,"%e",&rdummy);
    mhd8all[icount]=rdummy;
    icount=icount+1;
  }
  fclose(fp);
}
  printf("Num of data expected in text file %s : %d \n",instr,nx*ny*nz*8);
  printf("Num of data loaded                   : %d \n",icount);
/* adjust unit and scale */
  for(k = 0; k < nz; k++){for(j = 0; j < ny; j++){for(i = 0; i < nx; i++){
    int idrs = i + j * nx + k * nx * ny;
    mhd8all[idrs+0*nx*ny*nz] =   mhd8all[idrs+0*nx*ny*nz] / 1.0e8; // for density. 1e8 must be appearing in a string keyword BUNIT_000
    mhd8all[idrs+3*nx*ny*nz] = - mhd8all[idrs+3*nx*ny*nz];         // for V_theta, positive for northward, heliographic lat.
    mhd8all[idrs+6*nx*ny*nz] = - mhd8all[idrs+6*nx*ny*nz];         // for B_theta, positive for northward, heliographic lat.
  }}}

/* get some keyword values from the input Synchronic Mag series */
  double  lonfirst,lonlast,carrtime;
  char *trecstr, *tobsstr;
  int carrot, cols, rows;
  char *mapcver, *mapbldv;

  char strtmp[100];
  sprintf(strtmp,"%s[%s]",synostr,crdatestr);
  char *inQuery=strdup(strtmp);
  printf("now looking for synoptic map %s\n",inQuery);
  inRS = drms_open_records(drms_env, inQuery, &status);
  if (status || inRS->n == 0) {DIE("No input synoptic data found, so quit. Bye.");}
  int nrecs = inRS->n;
  if (nrecs > 1) {printf(" num. of data available is more than 1, %3d !!! Will use the first one anyway.\n",nrecs);}
  inRec = inRS->records[0]; // take first one, only. Hmmm. Maybe OK.
#if NEWMAP == 1
  inSeg = drms_segment_lookup(inRec,"Br"); // new series will have two or more segments
#else
  inSeg = drms_segment_lookupnum (inRec, 0); // take first segment data 
#endif
  cols = inSeg->axis[0];
  rows = inSeg->axis[1];
  lonfirst= drms_getkey_double(inRec,"LON_FRST",&status);
  lonlast = drms_getkey_double(inRec,"LON_LAST",&status);
  carrtime= drms_getkey_double(inRec,"CARRTIME",&status);
  printf(" LON_FRST of the syn.mag : %f\n",lonfirst);
  printf(" LON_LAST                : %f\n",lonlast);
  printf(" CARRTIME                : %f\n",carrtime);
{  
  TIME t_rec = drms_getkey_time(inRec,"T_REC",&status);
  TIME t_obs = drms_getkey_time(inRec,"T_OBS",&status);
  char timestr[26];
  sprint_time(timestr,t_rec,"TAI",0);
  trecstr=strdup(timestr);
  sprint_time(timestr,t_obs,"TAI",0);
  tobsstr=strdup(timestr);
  printf(" T_REC                   : %s\n",trecstr);
  printf(" T_OBS                   : %s\n",tobsstr);
}
  carrot = drms_getkey_int(inRec,"CAR_ROT",&status);
  printf(" CAR_ROT                 : %d\n",carrot);

#if NEWMAP == 1
  mapcver=drms_getkey_string(inRec,"CODEVER",&status);
  mapbldv=drms_getkey_string(inRec,"BLD_VERS",&status);
#endif

/* now making output record at designated series */
  printf(" Now open and create output at series %s\n",outstr);
  outRec = drms_create_record (drms_env, outstr, DRMS_PERMANENT, &status);
  if (!outRec) {fprintf (stderr, "Error creating record in series %s; abandoned\n",outstr);return 1;}

/* first of all, copy keywords from syn map. */
  drms_copykeys(outRec, inRec, 0, kDRMS_KeyClass_Explicit); // copy as much as keywords defined in both in and out.

/* we need to modify adjust some keyword, to adjust differences due to the differences in grid size etc. */  
  drms_setkey_string(outRec,"MAP_BLDV" ,mapbldv); 
  drms_setkey_string(outRec,"MAP_CVER" ,mapcver); 
#if NONDAILY == 1
  lonfirst = (double)(carrot-1)*360.0 + 360.0/((double)(nx))*0.5;
  lonlast  = (double)(carrot  )*360.0 - 360.0/((double)(nx))*0.5;
#else
  lonfirst = lonfirst;  // + (360.0/((double)(nx)) - 360.0/((double)(cols))); // need check!!!!!
  lonlast  = lonlast ;  // - (360.0/((double)(nx)) - 360.0/((double)(cols)));
#endif
  printf(" LON_FRST of the MHD     : %f\n",lonfirst);
  printf(" LON_LAST                : %f\n",lonlast);
  drms_setkey_double(outRec,"LON_FRST", lonfirst);
  drms_setkey_double(outRec,"LON_LAST", lonlast);
/* time of process */
{
  char timestr[26];
  TIME t_rec = drms_getkey_time(inRec,"DATE",&status);
  sprint_time(timestr,t_rec,"TAI",0);
  drms_setkey_string(outRec,"MAP_DATE" ,timestr);   // keep the date of syn.map creation
  sprint_time(timestr,CURRENT_SYSTEM_TIME,"UTC",1); // current time
  drms_setkey_string(outRec,"DATE" ,timestr);
}
/* some keywords whose value depend on grid number etc. */
  ddummy=((double)(nx))*0.5+0.5 ; drms_setkey_double(outRec,"CRPIX1",ddummy); // MIND check!!!
  ddummy=((double)(ny))*0.5+0.5 ; drms_setkey_double(outRec,"CRPIX2",ddummy);
  ddummy=((double)(nz))*0.5+0.5 ; drms_setkey_double(outRec,"CRPIX3",ddummy);
  ddummy=180.0+(double)(carrot-1)*360.0 ; drms_setkey_double(outRec,"CRVAL1",ddummy); // Carrington "Time"
  ddummy=0.0                            ; drms_setkey_double(outRec,"CRVAL2",ddummy); // equator be zero, in this context
  ddummy=(1.0+5.0)/2.0                  ; drms_setkey_double(outRec,"CRVAL3",ddummy); // 3.0 Rs
  ddummy=-360.0/((double)(nx))      ; drms_setkey_double(outRec,"CDELT1",ddummy);
  ddummy= 180.0/((double)(ny))      ; drms_setkey_double(outRec,"CDELT2",ddummy);
  ddummy= 4.000/((double)(nz))      ; drms_setkey_double(outRec,"CDELT3",ddummy);
  strdummy="degree"       ; drms_setkey_string(outRec,"CUNIT1", strdummy);
  strdummy="degree"       ; drms_setkey_string(outRec,"CUNIT2", strdummy);
  strdummy="solRad"       ; drms_setkey_string(outRec,"CUNIT3", strdummy);
  strdummy="CRLN-CAR"     ; drms_setkey_string(outRec,"CTYPE1", strdummy);
  strdummy="CRLT-CAR"     ; drms_setkey_string(outRec,"CTYPE2", strdummy);
  strdummy="HECR"         ; drms_setkey_string(outRec,"CTYPE3", strdummy);
  strdummy="3D-SPHERICAL" ; drms_setkey_string(outRec,"WCSNAME",strdummy);
  strdummy = mhdcorona_version() ; drms_setkey_string(outRec,"MHD_VER1",strdummy); // CVS version info. of this C-wrapper
#if NONDAILY == 1
  strdummy="v1.0, Jan. 10, 2014" ; drms_setkey_string(outRec,"MHD_VER2",strdummy); // version (and date) of MHD code (fortran)
#else
  strdummy="v1.0, Dec. 10, 2013" ; drms_setkey_string(outRec,"MHD_VER2",strdummy);
#endif

/* new keywords added to the MHD sereies */
{
  char *inputmapstr;
  char strtmp[100];
#if NONDAILY == 1
#if NEWMAP == 1
  sprintf(strtmp,"%s[%d]{Br}",synostr,carrot);
#else
  sprintf(strtmp,"%s[%d]",synostr,carrot);
#endif
#else
#if NEWMAP == 1
  sprintf(strtmp,"%s[%s]{Br}",synostr,trecstr);
#else
  sprintf(strtmp,"%s[%s]",synostr,trecstr);
#endif
#endif
  inputmapstr=strdup(strtmp);
  drms_setkey_string(outRec,"INPUTMAP", inputmapstr);
}
/* given and fixed words about MHD etc.*/
  strdummy = "3D version of the code in ApJS (2005) vol.161, 480";
  drms_setkey_string(outRec,"MHDMODEL",strdummy);
  strdummy = "Polytropic, with specific heat ratio 1.05.";
  drms_setkey_string(outRec,"MHD_SET1", strdummy);
  strdummy = " ";
  if (ichrbnd == 0) {strdummy=" sub-Alfvenic boundary treatment: case O";}
  if (ichrbnd == 1) {strdummy=" sub-Alfvenic boundary treatment: case AB'";}
  if (ichrbnd == 6) {strdummy=" sub-Alfvenic boundary treatment: case A";}
  drms_setkey_string(outRec,"MHD_SET2", strdummy);
  strdummy = " ";
  drms_setkey_string(outRec,"MHD_SET3", strdummy);
  if (isphindx < 0)
  {
    strdummy = "PFSS with iterative Laplace solver";
  }
  else
  {
    char strtmp[100];
    sprintf(strtmp,"PFSS with %d-order spherical harmonics polynomial",isphindx);
    strdummy = strdup(strtmp);
  }
  drms_setkey_string(outRec,"MHDIBMAG", strdummy);
  drms_setkey_int(outRec,"MHDMGIDX", isphindx);
/* comment */
  strdummy  = " ";
  drms_setkey_string(outRec,"COMMENT", strdummy);

/* per-segment keywords */
  drms_setkey_string(outRec,"BUNIT_000","1e8/cc"); // N .... UNIT or BUNIT ??
  drms_setkey_string(outRec,"BUNIT_001","1e6K");   // T
  drms_setkey_string(outRec,"BUNIT_002","km/s");   // V
  drms_setkey_string(outRec,"BUNIT_003","km/s"); 
  drms_setkey_string(outRec,"BUNIT_004","km/s"); 
  drms_setkey_string(outRec,"BUNIT_005","Gauss");  // B
  drms_setkey_string(outRec,"BUNIT_006","Gauss"); 
  drms_setkey_string(outRec,"BUNIT_007","Gauss"); 
  drms_setkey_string(outRec,"DSCRPT_000","number density");
  drms_setkey_string(outRec,"DSCRPT_001","temperature");
  drms_setkey_string(outRec,"DSCRPT_002","radial component of plasma flow");
  drms_setkey_string(outRec,"DSCRPT_003","latitudinal component of V (positive for northward)");
  drms_setkey_string(outRec,"DSCRPT_004","longitudinal component of V");
  drms_setkey_string(outRec,"DSCRPT_005","radial component of magnetic field");
  drms_setkey_string(outRec,"DSCRPT_006","latitudinal component of B (positive for northward)");
  drms_setkey_string(outRec,"DSCRPT_007","longitudinal component of B");

/* writing 8 segments */
  char *Resname[] = {"N","T","Vr","Vt","Vp","Br","Bt","Bp"};
  int ivar;
  for(ivar = 0; ivar < 8; ivar++)
  {
    float *dat1;
    int axes[3];
    axes[0] = nx;
    axes[1] = ny;
    axes[2] = nz;
    dat1 = (float *)calloc(nx*ny*nz, sizeof(float));
    int ndatsize=nx*ny*nz;
    for(n = 0; n < ndatsize; n++){dat1[n] = mhd8all[(n+ndatsize*ivar)];}
    outArray = drms_array_create (DRMS_TYPE_FLOAT, 3, axes, dat1, &status);
#if 0
    outArray->israw  = 0;
    outArray->bzero  = bzero[j]; // currently not used... saved as float
    outArray->bscale = bscal[j];
#endif
    outSeg = drms_segment_lookup (outRec, Resname[ivar]);

    if (!outSeg)
    {
      printf("Error in setting data segment %s\n", Resname[ivar]);
      DIE(" bye \n")
    }
    else
    {
      set_statistics(outSeg, outArray, 1); // 1 at the third argument is for "quick" run
      printf("Now writing segment : %s\n", Resname[ivar]);
      if (drms_segment_write (outSeg, outArray, 0))
      {
        printf("Error in writing to %dth segment %s\n", ivar,Resname[ivar]);
        DIE(" bye \n")
      }
    }
  } // end of ivar-loop
  printf("Done \n");

/* trailer and cleanup */
  drms_free_array(outArray);
  drms_close_record (outRec, DRMS_INSERT_RECORD);
  drms_close_records(inRS,  DRMS_FREE_RECORD); // here close the in Recs... OK?
  return (DRMS_SUCCESS);
} /* end of DoIt */

/* ----- */

char *mhdcorona_version(){return strdup("$Id: mhdtxt2jsoc_64cr.c,v 1.2 2014/03/13 18:59:30 keiji Exp $");}

/* ----------- end of this file ----------- */
