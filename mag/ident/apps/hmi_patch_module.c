/*
 * This jsoc module is in development to generate patches
 * from active region bitmaps
 *
 * turmon: created modeled on hmi_segment, additions
 * incorporated from xudong sun march 2010
 *
 */

#include "jsoc_main.h"

// mex and jsoc inter-operability
#include "mex.h"
#include "mexhead.h"
#include "mxt_jsoc.c"
// argument placements for the hmi_patch worker routine
#include "hmi_patch.h"
// indexing of per-patch statistics
#include "roi_stats_mag_defs.h"


// type for the mask that indicates region membership
typedef unsigned char HARP_region_t;

// number of pixels to pad AR with, all around.  nonnegative.
#define HARP_PIXEL_PAD 1

// standard trick to swallow the semicolon
// (ensure nonzero exit even if zero status)
#define DIE(msg) do { \
	fflush(stdout); \
        fprintf(stderr, "%s: %s, status=%d\n", module_name, msg, status); \
        return(status ? status : 1); \
        } while (0)
// report nonfatal error
//  (standard trick to swallow the semicolon)
//  TBD: should this affect the error status?
#define WARN(msg) do { \
	fflush(stdout); \
        fprintf(stderr, "%s: ERROR: %s.  Continuing.\n", module_name, msg); \
        } while (0)

// some standard JSOC stuff
#define	DTOR	(M_PI / 180.)
#define RAD2ARCSEC	(648000. / M_PI)
// Macros for WCS transformations.  
// ASSUMES: 
//   crpix1, crpix2 = CRPIX1, CRPIX2
//   sina,cosa = sin and cos of CROTA2 resp.
//   crvalx and crvaly are CRVAL1 and CRVAL2, 
//   cdelt = CDELT1 == CDELT2
// PIX_X and PIX_Y are CCD pixel addresses, 
// WX and WY are arc-sec W and N on the Sun from disk center.
#define PIX_X(wx,wy) ((((wx-crvalx)*cosa + (wy-crvaly)*sina)/cdelt)+crpix1)
#define PIX_Y(wx,wy) ((((wy-crvaly)*cosa - (wx-crvalx)*sina)/cdelt)+crpix2)
#define WX(pix_x,pix_y) (((pix_x-crpix1)*cosa - (pix_y-crpix2)*sina)*cdelt+crvalx)
#define WY(pix_x,pix_y) (((pix_y-crpix2)*cosa + (pix_x-crpix1)*sina)*cdelt+crvaly)


/* ################## Timing from T. Larson ################## */

static double getwalltime(void)
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec * 1000.0 + tv.tv_usec/1000.0;
}

static double getcputime(double *utime, double *stime)
{
  struct rusage ru;
  getrusage(RUSAGE_SELF, &ru);
  *utime = ru.ru_utime.tv_sec * 1000.0 + ru.ru_utime.tv_usec / 1000.0;
  *stime = ru.ru_stime.tv_sec * 1000.0 + ru.ru_stime.tv_usec / 1000.0;
  return *utime + *stime;
}

/********************************************************
 *
 * Boilerplate for input argument processing
 *
 *********************************************************/


/* 
 * command-line parameter names
 * (following the order in the driver arg list)
 */
// actually goes by "y" there...
#define kRecSetX        "x"
// ctr is determined from KW's
#define kParActive      "active"
#define kParKer         "ker"
#define kParKwt         "kwt"
#define kParTau         "tau"
// outputs
#define kSeriesOut      "bb"

char *module_name = "active region patches";
char *mex2c_progname = "hmi_patch_module";

// convolution kernel (gaussian on sphere)
static char patch_ker[] = 
  "[149.303191943,145.203125261,141.215651931,137.337679988,133.566202376,129.898294617,"
  "126.331112544,122.861890093,119.487937162,116.206637522,113.015446788,109.911890448,"
  "106.893561945,103.958120806,101.103290833,98.3268583348,95.6266704117,93.0006332846,"
  "90.4467106729,87.9629222148,85.5473419322,83.1980967365,80.9133649772,78.6913750282,"
  "76.5304039151,74.4287759784,72.3848615746,70.3970758123,68.4638773236,66.5837670684,"
  "64.7552871722,62.977019796,61.2475860363,59.5656448563,57.9298920457,56.3390592098,"
  "54.7919127855,53.2872530853,51.8239133663,50.4007589264,49.0166862237,47.6706220209,"
  "46.3615225536,45.0883727203,43.8501852956,42.6460001649,41.4748835792,40.3359274318,"
  "39.2282485537,38.1509880291,37.103310529,36.0844036638,35.0934773529,34.1297632128,"
  "33.1925139606,32.2810028348,31.3945230318,30.5323871579,29.6939266959,28.8784914871,"
  "28.0854492271,27.3141849754,26.5641006784,25.8346147063,25.1251614012,24.435190639,"
  "23.7641674031,23.1115713686,22.4768965,21.8596506579,21.2593552176,20.6755446985,"
  "20.1077664022,19.5555800625,19.0185575032,18.4962823064,17.9883494897,17.4943651921,"
  "17.0139463685,16.5467204926,16.0923252683,15.6504083488,15.2206270629,14.8026481497,"
  "14.3961475004,14.0008099062,13.6163288146,13.2424060915,12.8787517898,12.5250839249,"
  "12.1811282558,11.8466180725,11.5212939892,11.2049037433,10.8972019996,10.59795016,"
  "10.3069161789,10.0238743827,9.74860529543,9.4808954679,9.2205373127,8.967328943,"
  "8.72107401606,8.481581581,8.24866593066,8.0221464577,7.80184751447,7.58759827684,"
  "7.37923261175,7.17658894838,6.97951015286,6.78784340643,6.60144008692,6.42015565355,"
  "6.24384953482,6.07238501949,5.90562915063,5.74345262245,5.5857296801,5.43233802212,"
  "5.2831587056,5.13807605398,4.99697756732,4.85975383508,4.7262984513,4.59650793205,"
  "4.47028163522,4.34752168246,4.22813288329,4.1120226613,3.99910098234,3.88928028473,"
  "3.78247541135,3.6786035436,3.5775841372,3.47933885971,3.38379152983,3.29086805828,"
  "3.20049639037,3.11260645015,3.02713008602,2.94400101795,2.86315478601,2.78452870046,"
  "2.70806179309,2.63369476995,2.5613699654,2.49103129737,2.42262422387,2.35609570071,"
  "2.29139414038,2.22846937201,2.16727260251,2.10775637869,2.04987455052,1.99358223529,"
  "1.93883578283,1.88559274167,1.8338118261,1.7834528842,1.73447686665,1.68684579649,"
  "1.64052273964,1.59547177631,1.5516579731,1.50904735591,1.46760688364,1.42730442254,"
  "1.38810872129,1.34998938677,1.3129168605,1.27686239572,1.24179803509,1.20769658901,"
  "1.17453161456,1.14227739496,1.11090891967,1.08040186494,1.05073257501,1.02187804373,"
  "0.993815896726,0.966524374068,0.939982313368,0.914169133394,0.889064818098,0.864649901102,"
  "0.8409054506,0.81781305468,0.795354807045,0.77351329313,0.752271576595,0.731613186197,"
  "0.711522103016,0.69198274803,0.672979970042,0.654499033923,0.636525609194,0.619045758909,"
  "0.602045928848,0.585512937011,0.56943396339,0.553796540035,0.53858854138,0.523798174845,"
  "0.50941397169,0.495424778121,0.481819746644,0.46858832765,0.455720261237,0.443205569253,"
  "0.431034547563,0.419197758516,0.407686023633,0.396490416491,0.385602255795,0.375013098653,"
  "0.364714734023,0.354699176353,0.344958659382,0.335485630122,0.326272743,0.317312854164,"
  "0.308599015939,0.300124471445,0.291882649353,0.283867158793,0.276071784397,0.268490481477,"
  "0.261117371344,0.253946736743,0.246973017422,0.240190805824,0.233594842889,0.227180013977,"
  "0.220941344904,0.214873998084,0.208973268777,0.203234581442,0.197653486188,0.192225655325,"
  "0.186946880005,0.181813066964,0.176820235341,0.171964513597,0.16724213651,0.162649442256,"
  "0.158182869571,0.153838954985,0.149614330143,0.145505719188,0.141509936221,0.137623882835,"
  "0.133844545707,0.130168994265,0.126594378415,0.12311792633]";

ModuleArgs_t module_args[] =
{
   {ARG_STRING,  kRecSetX,   "",         "Input AR bitmap data records."},
   {ARG_FLOAT,   kParActive, "2",        "Active pixel label (RS)"},
   {ARG_FLOATS,  kParKer,    patch_ker,  "Smoothing kernel (RV)"},
   {ARG_FLOATS,  kParKwt,    "[4,1,1]",  "Kernel x,y,z weighting (RV(3))"},
   {ARG_FLOAT,   kParTau,    "0.018",    "Grouping threshold"},
   {ARG_STRING,  kSeriesOut, "",         "Output data series."},
   {ARG_INT,     "VERB",     "1",         "Verbosity: 0=errs/warns; 1=all"},
   {ARG_END}
};


/********************************************************
 *
 * Utility functions
 *
 *********************************************************/

/*
 * patch_fixup_patches: make the halfwidths of all patches
 * integers by bumping the patch size up by one in whatever
 * direction is possible.  If the patch is edge-to-edge, this
 * will be impossible, so in this case it is shrunken by one.
 * Also, pads AR's on all sides by HARP_PIXEL_PAD pixels.
 *
 * The patches are embedded in an MxN image, and all coordinates
 * are zero-based.
 */

static
void
patch_fixup_patches(double *bbs, int np, int M, int N)
{
  int p;
  int xsiz, ysiz;
  int pad = HARP_PIXEL_PAD; // a #define above

  // loop over patches
  for (p = 0; p < np; p++) {
    // take care of padding
    // push downward
    bbs[np*0+p] = bbs[np*0+p] - pad;
    if (bbs[np*0+p] < 0) bbs[np*0+p] = 0;
    bbs[np*1+p] = bbs[np*1+p] - pad;
    if (bbs[np*1+p] < 0) bbs[np*1+p] = 0;
    // push upward
    bbs[np*2+p] = bbs[np*2+p] + pad;
    if (bbs[np*2+p] >= M) bbs[np*2+p] = M-1;
    bbs[np*3+p] = bbs[np*3+p] + pad;
    if (bbs[np*3+p] >= N) bbs[np*3+p] = N-1;

    // number of pixels in this patch
    // (it is already an integer)
    xsiz = rint(bbs[np*2+p] - bbs[np*0+p] + 1);
    ysiz = rint(bbs[np*3+p] - bbs[np*1+p] + 1);
    // x: if even, make it odd
    if (xsiz % 2 == 0) {
      // attempt to expand by one; contract by one if that fails
      if (bbs[np*0+p] > 0)
	bbs[np*0+p] -= 1.0; // push start back by one
      else if (bbs[np*2+p] < M-1)
	bbs[np*2+p] += 1.0; // push end out by one
      else
	bbs[np*0+p] += 1.0; // pull start up by one
    }
    // y: if even, make it odd
    if (ysiz % 2 == 0) {
      // attempt to expand by one; contract by one if that fails
      if (bbs[np*1+p] > 0)
	bbs[np*1+p] -= 1.0; // push start back by one
      else if (bbs[np*3+p] < N-1)
	bbs[np*3+p] += 1.0; // push end out by one
      else
	bbs[np*1+p] += 1.0; // pull start up by one
    }
  }
  return;
}

/*
 * patch_free_bitmaps: free a list of bitmaps allocated 
 * by patch_extract_bitmaps
 */
static
void
patch_free_bitmaps(HARP_region_t **yps, int np)
{
  int p;

  for (p = 0; p < np; p++)
    free(yps[p]);
}

/*
 * patch_extract_bitmaps: extract a list of bitmaps around a given
 * list of bounding boxes.  bitmaps are extracted from MxN image Img;
 * there are np bounding boxes.
 * Rgn is encoded as:
 *    NaN/0  --> offdisk
 *    1..np  --> patch number
 * Mask is encoded as:
 *    NaN/0  --> offdisk
 *    1      --> not active
 *    2      --> active
 */
static
int
patch_extract_bitmaps(HARP_region_t **yps, double *Mask, double *Rgn, 
		      int M, int N, double *bbs, int np)

{
  double R1, M1;     // image values
  HARP_region_t L1;  // patch label value
  HARP_region_t *p1; // head of one patch
  int p;             // patch counter
  int x0, y0;        // patch origin
  int xsiz, ysiz;    // number of pixels in patch
  int x, y;          // coordinate counters
  
  // all yps = NULL, signaling not allocated yet
  bzero(yps, np*sizeof(*yps));
  // loop over patches
  for (p = 0; p < np; p++) {
    // start of this patch
    x0 = bbs[np*0 + p];
    y0 = bbs[np*1 + p];
    // number of pixels in this patch
    xsiz = rint(bbs[np*2+p] - bbs[np*0+p] + 1);
    ysiz = rint(bbs[np*3+p] - bbs[np*1+p] + 1);
    // copy the bitmap over the patch
    p1 = yps[p] = calloc((size_t) (xsiz * ysiz), sizeof(**yps));
    if (p1 == NULL) 
      break; // failure
    for (x = 0; x < xsiz; x++)
      for (y = 0; y < ysiz; y++) {
	// FIXME: this region encoding can be improved
	R1 = Rgn [(y0+y)*M + (x0+x)];
	M1 = Mask[(y0+y)*M + (x0+x)];
	// This is a bit unwieldy
	if (R1 == (p+1))
	  // note, > 128 does not fit in a char
	  L1 = (HARP_region_t) (64 + p + 1); // map region #p to 64+p
	else if (R1 > 0)
	  L1 = 3; // other regions (HARPs) -> 3
	else if (M1 == 2)
	  L1 = 2; // other ARs pixels (not HARP) -> 2
	else if (M1 == 1)
	  L1 = 1; // on-disk, but not AR/HARP at all
	else
	  L1 = 0; // off-disk
        // plug it in
	p1[y*xsiz+x] = L1;
      }
  }
  // error iff p1 == NULL
  if (p1 == NULL) {
    // tidy up the partial allocations
    for (p = 0; p < np; p++)
      if (yps[p])
	free(yps[p]);
    return 1;
  } else {
    // return status: zero for ok
    return 0;
  }
}    


/********************************************************
 *
 * Driver function
 *
 *********************************************************/

int DoIt(void)
{
  // status and input arguments
  int status = DRMS_SUCCESS;
  CmdParams_t *params = &cmdparams;
  const char *xRecQuery = params_get_str(params, kRecSetX);
  const char *yRecQuery = params_get_str(params, (char *) kSeriesOut);
  // DRMS data
  DRMS_RecordSet_t *xRD, *yRD;
  DRMS_Record_t  *xRec,  *yRec;
  DRMS_Record_t *link_mg, *link_pg;
  DRMS_Link_t *link_mask;
  DRMS_Segment_t *xSeg,   *ySeg,   *magSeg;
  DRMS_Array_t   *xArray, *yArray, *magArray;
  
  // mexFunction arguments
  // (the constants are defined in hmi_patch.c)
  const int nrhs = MXT_Hpat_NARGIN_MAX;     // MIN == MAX
  const int nlhs = MXT_Hpat_NARGOUT_MAX;    // bounding boxes + labeled HARP
  mxArray *prhs[MXT_Hpat_NARGIN_MAX];       // input args (too long is OK)
  mxArray *plhs[MXT_Hpat_NARGOUT_MAX];      // output args
  double *x_tmp, *mag_tmp;                  // holding place for array data

  // other variables
  int ds, Nds;                     // iteration
  int M, N, pDims[2];              // image sizes
  double rSun, x0, y0, p0, beta;   // KWs
  int numPatch;                    // number of HARPs found
  int p;                           // patch counter
  int sn;                          // keyword counter
  double *bbs;                     // bounding box array (np-by-4)
  double *stats;                   // stats array (np-by-RS_num_stats)
  HARP_region_t **bmps;            // points to each HARP mask
  int verbflag;
  TIME t_rec;
  char time_buf[64]; // time string for t_rec
  char *msg;         // error message
  // wcs info
  float crvalx, crvaly, cdelt, crpix1, crpix2, crota2, crlt, sina, cosa;
  double rsun_ref, dsun_obs;
  // Time measuring -- ugly, but functional
  double wt0, wt1, wt;
  double ut0, ut1, ut;
  double st0, st1, st;
  double ct0, ct1, ct;
  // for Runtime1 keyword
  time_t now;

  // timings
  wt0 = getwalltime();
  ct0 = getcputime(&ut0, &st0);

  // verbosity level
  verbflag = params_get_int(&cmdparams, "VERB");

  // set up fixed arguments via mxt_param... calls
  // other parameters are allocated here, but set up from KWs
  prhs[MXT_Hpat_ARG_ctr]   = mxCreateDoubleMatrix(1, 3, mxREAL); // from KWs
  prhs[MXT_Hpat_ARG_p0]    = mxCreateDoubleMatrix(1, 1, mxREAL); // from KWs
  prhs[MXT_Hpat_ARG_beta]  = mxCreateDoubleMatrix(1, 1, mxREAL); // from KWs
  prhs[MXT_Hpat_ARG_active]= mxt_param_to_scalar(params, kParActive);
  prhs[MXT_Hpat_ARG_ker]   = mxt_param_to_matrix(params, kParKer, 1, -1);
  prhs[MXT_Hpat_ARG_kwt]   = mxt_param_to_matrix(params, kParKwt, 1, -1);
  prhs[MXT_Hpat_ARG_tau]   = mxt_param_to_scalar(params, kParTau);

  // check args
  if (!prhs[MXT_Hpat_ARG_ctr]    || 
      !prhs[MXT_Hpat_ARG_p0]     || 
      !prhs[MXT_Hpat_ARG_beta]   || 
      !prhs[MXT_Hpat_ARG_active] || 
      !prhs[MXT_Hpat_ARG_ker]    || 
      !prhs[MXT_Hpat_ARG_kwt]    || 
      !prhs[MXT_Hpat_ARG_tau])
    DIE("failed to set up a miscellaneous parameter");

  // open activity mask series
  xRD = drms_open_records(drms_env, (char *) xRecQuery, &status);
  if (status || xRD->n == 0)
    DIE("No AR bitmap input data found");
  Nds = xRD->n;

  // open patches

  // Xudong Mar 24: One bitmap corresponds to multiple patch records, 
  // and the number is undetermined. So we can't create
  // Nds records at a time. We'll create them one by one in the loop.

  // loop over full-disk activity masks
  for (ds = 0; ds < Nds; ds++) {
    printf("begin record %d of %d\n", ds+1, Nds);
    fflush(stdout);
    // mark the time
    wt1 = getwalltime();
    ct1 = getcputime(&ut1, &st1);

    // input record and time
    xRec = xRD->records[ds];
    t_rec = drms_getkey_time(xRec, "T_REC", &status);
    // summary
    sprint_time(time_buf, t_rec, "TAI", 0);
    printf("\tprocessing image taken at %s\n", time_buf);
    fflush(stdout);

    // Ephemeris courtesy of Phil's macros above
    // No error checking for now
    crvalx = drms_getkey_float(xRec, "CRVAL1",   &status); // center in arcsec
    crvaly = drms_getkey_float(xRec, "CRVAL2",   &status);
    cdelt  = drms_getkey_float(xRec, "CDELT1",   &status); // arcsec, assumimg dx=dy
    crpix1 = drms_getkey_float(xRec, "CRPIX1",   &status); // disk center on ccd
    crpix2 = drms_getkey_float(xRec, "CRPIX2",   &status);
    crota2 = drms_getkey_float(xRec, "CROTA2",   &status); // rotation
    crlt   = drms_getkey_float(xRec, "CRLT_OBS", &status); // rotation
    if (status) {
      WARN("could not get WCS from labeling");
      continue;
    }
    sina = sin(crota2 * DTOR); 
    cosa = cos(crota2 * DTOR); 
    p0 = -crota2;
    beta = crlt;
    x0 = PIX_X(0.0,0.0) - 1.0; // zero-based coordinates
    y0 = PIX_Y(0.0,0.0) - 1.0;
    rsun_ref = drms_getkey_double(xRec, "RSUN_REF", &status);
    if (status) {
      rsun_ref = 6.96e8; // this is OK
      status = 0;
    }
    dsun_obs = drms_getkey_double(xRec, "DSUN_OBS", &status);
    if (status) {
      WARN("could not get DSUN_OBS from mask");
      continue;
    }
    rSun = asin(rsun_ref / dsun_obs) * RAD2ARCSEC / cdelt;
    if (verbflag)
      printf("\tX0=%.3f, Y0=%.3f, RSUN=%.3f, P0=%.2f, beta=%.2f\n", 
             x0, y0, rSun, p0, beta);

    // find links to source data; used later on
    if ((link_mg = drms_link_follow(xRec, "MDATA", &status)) == NULL)
      DIE("require the MDATA link to get patch statistics");
    if ((link_pg = drms_link_follow(xRec, "PDATA", &status)) == NULL)
      WARN("will be unable to set PDATA link on patches");

    // look up AR mask
    xSeg = drms_segment_lookupnum(xRec, 0);
    if (!xSeg)
      DIE("problem getting AR bitmap");
    // store as a small integer, but read in as double
    xArray = drms_segment_read(xSeg, DRMS_TYPE_DOUBLE, &status);
    if (status)
      DIE("problem reading AR bitmap");
    M = xArray->axis[1]; N = xArray->axis[0];

    // set up varying arguments: AR mask
    if ((prhs[MXT_Hpat_ARG_y] = mxCreateDoubleMatrix(M, N, mxREAL)) == NULL)
      DIE("failed calloc for input mask");
    x_tmp = mxGetPr(prhs[MXT_Hpat_ARG_y]); // hold data segment
    mxSetPr(prhs[MXT_Hpat_ARG_y], (double *) xArray->data);

    // look up magnetogram
    magSeg = drms_segment_lookupnum(link_mg, 0);
    if (!magSeg)
      DIE("problem getting mgram corresponding to AR mask");
    // read in as double
    magArray = drms_segment_read(magSeg, DRMS_TYPE_DOUBLE, &status);
    if (status)
      DIE("problem reading magnetogram");
    if ((M != magArray->axis[1]) || (N != magArray->axis[0]))
      DIE("dimension mismatch: mask != magnetogram");

    // set up varying arguments: magnetogram
    if ((prhs[MXT_Hpat_ARG_mag] = mxCreateDoubleMatrix(M, N, mxREAL)) == NULL)
      DIE("failed calloc for input magnetogram");
    mag_tmp = mxGetPr(prhs[MXT_Hpat_ARG_mag]); // hold data segment
    mxSetPr(prhs[MXT_Hpat_ARG_mag], (double *) magArray->data);

    // set up varying arguments: center
    //   center was allocated outside the loop
    // mexfunction coordinates are one-based (origin at (1,1)), so add 1.
    // mexfunction ctr[1] direction varies fastest, ctr[0] direction
    // varies slowest, so flip around.
    mxGetPr(prhs[MXT_Hpat_ARG_ctr])[0] = y0+1;
    mxGetPr(prhs[MXT_Hpat_ARG_ctr])[1] = x0+1;
    mxGetPr(prhs[MXT_Hpat_ARG_ctr])[2] = rSun;
 
    // p-angle and beta, also allocated outside the loop
    mxGetPr(prhs[MXT_Hpat_ARG_p0  ])[0] = p0;
    mxGetPr(prhs[MXT_Hpat_ARG_beta])[0] = beta;

    // print a sample value
    if (verbflag)
      printf("\tMask[%d,%d] = %f\n", M/2, N/2, 
	     ((double *) xArray->data)[M/2+N/2*M]);

    // call the function
    printf("patch_module: hmi_patch calling mexfunction.\n");fflush(stdout);
    msg = main_hmi_patch(nlhs, plhs, nrhs, prhs);
    if (msg == NULL) {
      printf("patch_module: hmi_patch returned OK from mexfunction.\n");
      fflush(stdout);
    } else {
      /* error within called function */
      printf("main loop: returned via error from mexfunction.\n");
      /* fail */
      status = 1;
      DIE("mexFunction failed (mexErrMsgTxt)");
    }
    // the head of the bounding box list, for clarity
    // in contrast with the center, these bounding boxes have 
    // zero-based coordinates: origin at (0,0)
    // the origin is an argument to the bounding box generator, and
    // it's easiest to set it up that way.
    bbs = mxGetPr(plhs[MXT_Hpat_ARG_bb]);
    // resize bounding boxes so they have integer halfwidths
    numPatch = mxGetM(plhs[MXT_Hpat_ARG_bb]);
    patch_fixup_patches(bbs, numPatch, M, N);

    // the head of the patch statistics
    stats = mxGetPr(plhs[MXT_Hpat_ARG_stats]);

    // extract patches from bitmap
    bmps = calloc(numPatch, sizeof(*bmps)); // pointers to per-patch bitmaps
    if (0 && verbflag) {
      // can load directly into ds9 as a check
      printf("# Region file format: DS9 version 3.0\n");
      for (p = 0; p < numPatch; p++)
        printf("image; box %6.1f %6.1f %6.1f %6.1f\n", 
	       1+(bbs[0*numPatch+p] + bbs[2*numPatch+p])/2,  // ycen, ds9 coords
	       1+(bbs[1*numPatch+p] + bbs[3*numPatch+p])/2,  // xcen, ds9
	       bbs[2*numPatch+p] - bbs[0*numPatch+p] + 1,
	       bbs[3*numPatch+p] - bbs[1*numPatch+p] + 1);
    }
    // args: bitmap-list, mask, mask sizes, bounding boxes, number of boxes
    printf("patch_module: extracting %d maps.\n", numPatch);
    fflush(stdout);
    status = patch_extract_bitmaps(bmps, 
				   mxGetPr(prhs[MXT_Hpat_ARG_y]),   // mask
				   mxGetPr(plhs[MXT_Hpat_ARG_yrgn]),// region tag
				   M, N, bbs, numPatch);
    if (status) {
      DIE("failed (bad calloc during patch extraction)");
    }
    printf("patch_module: maps extracted.\n");
    fflush(stdout);

    // write patches 
    // (create multiple output records for one input record)

    // abbreviations for each corner coordinate
    double *yLO = bbs;                  // lower left y
    double *xLO = bbs + 1*numPatch;	// lower left x
    double *yHI = bbs + 2*numPatch;    	// upper right y
    double *xHI = bbs + 3*numPatch;     // upper right x

    for (p = 0; p < numPatch; p++) {
      // One at a time: "create_record" instead of "create_records"
      yRec = drms_create_record(drms_env, (char *) yRecQuery, DRMS_PERMANENT, &status);
      if (status)
	DIE("Output recordset for patches not created");
      printf("\tpatch: looking up segment\n");
      ySeg = drms_segment_lookupnum(yRec, 0);
        
      // get sizes (as integers)
      pDims[0] = rint(yHI[p] - yLO[p] + 1.0);
      pDims[1] = rint(xHI[p] - xLO[p] + 1.0);

      // NB: a subsequent drms_free of yArray would free bmps[p] also!
      // (this is not done now; we free bmps[...] separately)
      printf("\tpatch: array_create (patch(0,0) = %d)\n", (int) bmps[p][0]);
      yArray = drms_array_create(DRMS_TYPE_CHAR, 2, pDims,
				 bmps[p], &status);
      if (status)
	DIE("problem creating drms array for patches");                       
        
      ySeg->axis[0] = yArray->axis[0];
      ySeg->axis[1] = yArray->axis[1];

      yArray->parent_segment = ySeg;
    
      // get mask to insert link
      // FIXME: there must be a better way
      if ((link_mask = hcon_lookup_lower(&yRec->links, "MASK")) == NULL)
	WARN("will be unable to set MASK link on patches");
    
      // insert links to source data
      if (link_mg)
	drms_setlink_static(yRec, "MDATA", link_mg->recnum);
      if (link_pg)
	drms_setlink_static(yRec, "PDATA", link_pg->recnum);
      if (link_mask)
        drms_setlink_static(yRec, "MASK", xRec->recnum);
            
      // Essential prime keys
      drms_copykey(yRec, xRec, "T_REC");
      drms_setkey_int(yRec, "PNUM", p+1); // number from 1
		
	  // date and build version
	  // updated by xudong oct 3 2010
	  drms_setkey_string(yRec, "BLD_VERS", jsoc_version);
	  drms_setkey_time(yRec, "DATE", CURRENT_SYSTEM_TIME);
        
      // Geometry
      drms_setkey_int(  yRec, "HWIDTH1", (pDims[0]-1)/2);
      drms_setkey_int(  yRec, "HWIDTH2", (pDims[1]-1)/2);
      drms_setkey_float(yRec,  "CRPIX1", (yHI[p]+yLO[p])/2);
      drms_setkey_float(yRec,  "CRPIX2", (xHI[p]+xLO[p])/2);
      // probably should be set as a link
      drms_setkey_float(yRec,  "CDELT1", (float) cdelt);
      drms_setkey_float(yRec,  "CDELT2", (float) cdelt);

      // Flags
      drms_setkey_int(yRec, "PATCHNUM", numPatch);
      drms_setkey_int(yRec, "MASK", 64+p+1);
            
      // Insert all patch statistics
      for (sn = 0; sn < RS_num_stats; sn++) {
	const char *name1 = RS_index2keyname[sn];
	double stat1 = stats[sn*numPatch+p];
	// insert if HMI keyword name is non-null
	if (name1) {
	  switch (*name1) {
	  case 'i':
	    drms_setkey_int(yRec, name1+1, (int) stat1);
	    break;
	  case 'f':
	    drms_setkey_float(yRec, name1+1, (float) stat1);
	    // printf("\tP %d:\t%s -> %.3g\n", p, name1+1, stat1);
	    break;
	  default:
	    WARN("unknown key type (first char of HMI keyname)");
	    break;
	  }
	}
      } // end for sn

      // Write the complete patch
      status = drms_segment_write(ySeg, yArray, 0);
      if (status)
	DIE("problem writing drms segment for patches");
      drms_close_record(yRec, DRMS_INSERT_RECORD);
    }

    // free the bitmaps
    patch_free_bitmaps(bmps, numPatch); // component bitmaps
    free(bmps); // pointer array, calloc'd above

    // free the AR mask
    // (put back the original, unused, data segment; then free)
    mxSetPr(prhs[MXT_Hpat_ARG_y], (double *) x_tmp);
    mxDestroyArray(prhs[MXT_Hpat_ARG_y]);

    // free the mgram
    // (put back the original, unused, data segment; then free)
    mxSetPr(prhs[MXT_Hpat_ARG_mag], (double *) mag_tmp);
    mxDestroyArray(prhs[MXT_Hpat_ARG_mag]);

    // region codes
    mxDestroyArray(plhs[MXT_Hpat_ARG_yrgn]); 
    // bounding boxes
    mxDestroyArray(plhs[MXT_Hpat_ARG_bb]); 
    // region statistics
    mxDestroyArray(plhs[MXT_Hpat_ARG_stats]); 

    // free inArray (with its data segment)
    if (verbflag) {
      printf("patch_module: within-loop drms_free\n");
      drms_free_array(xArray); 
    }

    // TODO: copy over more keywords?
    // e.g., total number of patches?
    // drms_copykey(yRec, xRec, "T_REC");

    // Timing info
    if (verbflag) {
      wt = getwalltime();
      ct = getcputime(&ut, &st);
      printf("record %d time used: %.3f s wall, %.3f s cpu\n", 
	     ds+1, (wt - wt1)*1e-3, (ct - ct1)*1e-3);fflush(stdout);
    }
    printf("record %d done\n", ds+1); fflush(stdout);
  } 
  printf("segment_module: close record (input)\n");
  drms_close_records(xRD, DRMS_FREE_RECORD);
  //drms_close_records(yRD, DRMS_INSERT_RECORD);

  if (verbflag) {
    wt = getwalltime();
    ct = getcputime(&ut, &st);
    printf("total time used: %.3f s wall, %.3f s cpu\n", 
	   (wt - wt0)*1e-3, (ct - ct0)*1e-3);
  }

  // free the unchanging arguments
  mxDestroyArray(prhs[MXT_Hpat_ARG_ctr]);
  mxDestroyArray(prhs[MXT_Hpat_ARG_p0]);
  mxDestroyArray(prhs[MXT_Hpat_ARG_beta]);
  mxDestroyArray(prhs[MXT_Hpat_ARG_active]);
  mxDestroyArray(prhs[MXT_Hpat_ARG_ker]);
  mxDestroyArray(prhs[MXT_Hpat_ARG_kwt]);
  mxDestroyArray(prhs[MXT_Hpat_ARG_tau]);

  // exit OK
  return DRMS_SUCCESS;
}

