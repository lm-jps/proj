/* standalone code for reading SUVI fits files and fitting EUV limb */
/* just simple C all here in one file, you can add classes and methods using the
functions here as you desire */

/* you can build this on most systems with just 
  cc -o suvilimbs suvilimbs.c aia_despike.c -lm
*/

/* major components are functions to read simple fits files, the limbcompute
function with its attendants, and some boiler plate to read the SUVI images and
output the results in a list, main is at the end */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <ctype.h>
#define	HEAD_LIMIT	100

#include "../../lev0/apps/aia_despike.c"

/*------------------------------------------------------------------------------------------------*/
int aia_despike(
     int *array, unsigned char *mask, int nx, int ny, float frac, int level, int niter, /* void **outarray, */
     int *badblobs, int sizeofbads,
     int *nspikes, int **oldvalues, int **spikelocs, int **newvalues);
/*------------------------------------------------------------------------------------------------*/
int limbcompute(float *x, int nx, int ny, float xcguess, float ycguess, float rguess, float rrange,
  int limbmode, int useprevflag, float fwhm);
/* fits reader and friends */
/*------------------------------------------------------------------------------------------------*/
int wrfits2d(void *x, char *name, int nx, int ny, int ana_type);
int rdfits2d(void **x, char *name, int *nx, int *ny, int *ana_type, float *exp, float *wave, int *gt_y, int *gt_z);
int fits_fatality(FILE *fin);
int fits_problems(int i);
int verbosity=1;
void swapd(char x[],int n);
void swapl(char x[],int n);
void swapb(char x[],int n);
char *bigger_header(int n, char *head);
int	fits_head_malloc_flag=0;
char	*fitshead;
int	ana_type_size[] = {1,2,4,4,8};
FILE	*fopen(), *fin, *fout;
 /*--------------------------------------------------------------------------*/
/*------------------------------------------------------------------------- */
/* limbcompute function and support functions here */
/*------------------------------------------------------------------------------------------------*/
extern float sdisk_xc, sdisk_yc, sdisk_r;
int scounts ;
typedef unsigned char byte;
#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define ABS(x) ((x)>=0?(x):-(x))
#define ALN2I 1.442695022
#define TINY 1.0e-5
#define MAXNITER 100
/* arrays showing convergence, length is MAXNITER */
float xcenters[MAXNITER], ycenters[MAXNITER], radiai[MAXNITER];
long	pcounts[MAXNITER], nextcenter;
extern int max4histogram;
void mmq(float *x, int n, float *max, float *min);
void indexf(int n, float ra[],int indx[]);
long getindexofmax(int *x, int n);
int pit(float *qxbase, float *qybase, int nxq, int npow, double *a, double *cfbase, double *fbase );
float parabolicmax5(int *x, int nx);
void sobel_index(float *x, int nx, int ny, long *iq, int nindex, float **grad, float **ang);
void mm_f(float *x, int n, float *max, float *min);
void mm_int(int *x, int n, int *max, int *min);
void indexf(int n, float ra[],int indx[]);
int circle_fit_lsq(float *x, float *y, float *w, int n, float *xc_out, float *yc_out, float *r_out, float *chi_out);
int key_locations(int *array, int n, int *nv, int **counts, int **values, int ***indexptrs);
long hist_max_xvalue(float *x, int n);
int *hist_array(float *x, int n, int *nh, int *hmin);
int hist_max_freq(float *x, int n);
void getmin9(int *p, int ix, int iy, int nx, float *x0, float *y0);
void zeroinonlimb(float *xlimb, float *ylimb, int nin, float dr, float *g);
void  minihough(int nc, float *xlimb, float *ylimb, int niter);
void decomp(double *x, int n, int nd);
void solve(double *a, double *b, int n, int nd);
int pit(float *qxbase, float *qybase, int nxq, int npow, double *a, double *cfbase, double *fbase );
float parabolicmax5(int *x, int nx);
void printarray(float *x, int n);
void printintarray(int *x, int n);
#define EPSILON 1.e-10
double systime();				/* internal systime */
/*------------------------------------------------------------------------- */
/* main here */
/*------------------------------------------------------------------------------------------------*/
int main(int argc,char *argv[])
  {
  char *progname = NULL;
  char *listname = NULL;
  char *basePath = NULL;
  char *outname = NULL;
  char listnameDefault[] = "/home/suvitc/LimbFit/inp/test/filelist";
  //char basePathDefault[] = "/titan/goes/suvi/GOES-16/fits/2017/02/22/H0000/";
  char basePathDefault[] = "";
  char outnameDefault[] = "/home/suvitc/LimbFit/out/test.limb_fits";
  char samplename[] = "suvi_fm1_l0_171_light_20170222_005004_00207151.fits";
  float *imp, *impf;
  float rguess, xcguess, ycguess, rrange, fwhm, exp, min_exp = 0.5, wave;
  float gt_y_px, gt_z_px;
  int gt_y, gt_z;
  int nx, ny, nl, ana_type, nxfits, nyfits, i, ns;
  double tstart, tend;
  int listflag, singleflag, status, limbmode=0, useprevflag=0, ifit, ireject;
  int auto_guess = 1, despike = 1;
  char  *ffname, *sq;
  FILE *listfile, *imfile, *outfile;
  tstart = systime();
  /* some initial defaults */
  rrange = 50.;  rguess = 390.;  xcguess = ycguess = 641.;   fwhm = 3.0;
   /* process command line arguments */
   progname = argv[0];
   if (argc>=2) {
    for (i=1;i<argc;i++) {
     if (argv[i][0] == '-') {
       switch (argv[i][1])
       {
        case 'G': auto_guess = 0; break;

	case 'S': despike = 0; break;

	case 'r': /* radius guess */
	  ns = sscanf(argv[i+1],"%f", &rguess); i++;
	  if (ns != 1) {
	     fprintf(stderr,"problem parsing r\n");
	     exit(1); }
	  break;

	case 'x': /* xc guess */
	  ns = sscanf(argv[i+1],"%f", &xcguess); i++;
	  if (ns != 1) {
	     fprintf(stderr,"problem parsing x\n");
	     exit(1); }
	  break;

	case 'y': /* yc guess */
	  ns = sscanf(argv[i+1],"%f", &ycguess); i++;
	  if (ns != 1) {
	     fprintf(stderr,"problem parsing y\n");
	     exit(1); }
	  break;

	case 'd': /* radius delta or range */
	  ns = sscanf(argv[i+1],"%f", &rrange); i++;
	  if (ns != 1) {
	     fprintf(stderr,"problem parsing delta r\n");
	     exit(1); }
	  break;

	case 'm': /* mode (for limbmode) */
	  ns = sscanf(argv[i+1],"%d", &limbmode); i++;
	  if (ns != 1) {
	     fprintf(stderr,"problem parsing mode r\n");
	     exit(1); }
	  break;

	case 'e': /* minimum exposure time */
	  ns = sscanf(argv[i+1],"%f", &min_exp); i++;
	  if (ns != 1) {
	     fprintf(stderr,"problem parsing mode r\n");
	     exit(1); }
	  break;

	case 'i':  /* input single file, make sure we have a name here */
	  if ( (i+1) >= argc || argv[i+1][0] == '-') {
	    fprintf(stderr,"no input file name after -i key (or it begins with a -)\n");
	    fprintf(stderr,"%s fatal error\n", progname);
	    exit(1); }
	  ffname = argv[i+1];
	  singleflag = 1;  listflag = 0;
          break;

	case 's': /* use sample image */
	  ffname = samplename;
	  singleflag = 1;  listflag = 0;
	  break;

	case 'l':  /* input list name, make sure we have a name here */
	  if ( (i+1) >= argc || argv[i+1][0] == '-') {
	    fprintf(stderr,"no input file name after -i key (or it begins with a -)\n");
	    fprintf(stderr,"%s fatal error\n", progname);
	    exit(1); }
	  listname = argv[i+1];
 	  singleflag = 0;  listflag = 1;
         break;

	case 'o':  /* output file name, make sure we have a name here */
	  if ( (i+1) >= argc || argv[i+1][0] == '-') {
	    fprintf(stderr,"no output name after -o key (or it begins with a -)\n");
	    fprintf(stderr,"%s fatal error\n", progname);
	    exit(1); }
	  outname = argv[i+1];
          break;

	case 'p':  /* optional path for image files, make sure we have a name here */
	  if ( (i+1) >= argc || argv[i+1][0] == '-') {
	    fprintf(stderr,"no image path after -p key (or it begins with a -)\n");
	    fprintf(stderr,"%s fatal error\n", progname);
	    exit(1); }
	  basePath = argv[i+1];
          break;

        case 'q': verbosity--; break;

        case 'v': verbosity++; break;

       }
     }
    }
  }
  /* some defaults if we didn't have them on the command line */
  if (listname == 0) { listflag = 0; } else { listflag = 1; }
  if (basePath == 0) basePath = basePathDefault;
  if (outname == 0) outname = outnameDefault;
  /* check if anything to do */
  if (singleflag == 0 && listflag == 0) {
    fprintf(stderr,"suvilimbs usage\n\nreads one or more SUVI images and attempts to fit the limb(s)\n");
    fprintf(stderr,"results output one line per file, images with short exposure times\n");
    fprintf(stderr,"can be rejected\n\ncommand line parameters\n\n");
    fprintf(stderr,"-x guess for x center (default is 641)\n-y guess for y center (default is 641)\n");
    fprintf(stderr,"-r guess for radius (default is 390)\n-d annulus half width (default is 50)\n");
    fprintf(stderr,"-m limbmode (default is 0)\n-e minimum exposure time (default is 0.5)\n");
    fprintf(stderr,"-i name of single file (overrides -l)\n-l name for list of file names (overrides -i)\n");
    fprintf(stderr,"-o output file name\n-p optional path for file names (added to input name)\n");
    fprintf(stderr,"\nexample:  suvilimbs -l suvi.list -o limbfits2.txt -x 640 -y 640 -r 389\n");
    exit(1); }
  /* setup an output file */
  outfile = fopen(outname, "w");
  ifit = ireject = 0;
  if (singleflag) {
    sq = (char *) malloc(strlen(basePath) + strlen(ffname)  + 10);
    if (strlen(basePath) >0) {
      strcpy(sq, basePath);  strncat(sq, "/", 1); strncat(sq,ffname, 511);
    } else { strcpy(sq, ffname); }
    /* add a fits to the end if not there already */
    if (strstr(sq, ".fits") == NULL) strncat(sq, ".fits", 5);
    if (verbosity>1) fprintf(stdout,"final name %s\n", sq);
//    status = rdfits2d((void **) &imp, sq, &nxfits, &nyfits, &ana_type, &exp, &wave, &gt_y, &gt_z);
    if (status)
      { fprintf(stderr, "error reading fits file %s\n", ffname); exit(1); }
    if (verbosity>2) fprintf(stdout,"wave, exp = %f, %f\n", wave, exp);
    if (auto_guess) {
      /* SUVI pix = 180.0e-9 * nano_radians / pi * 3600 / 2.5  */
      gt_y_px = gt_y * 8.25059225e-5;
      gt_z_px = gt_z * 8.25059225e-5;
      xcguess = gt_y_px + 641.0;
      ycguess = gt_z_px + 641.0;
    }
    /* check if we meet minimum exposure time */
    if (exp >= min_exp || exp < 0) {
      if (despike) {
        int *impi, *oldvals, *spikelocs, *newvals, nspikes;
        /* aia_despike requires I*4 input for now */
        if (ana_type != 2) {
          int	npix = nxfits * nyfits;
          impi = (int *) malloc(npix * sizeof(int));
          if (npix <= 0) {
          fprintf(stderr, "empty fits file? # pixels = %d\n", npix); exit(1); }
          if (ana_type == 1) {
            short *ps = (short *) imp;
            int *pi = impi;
            while (npix--) { *pi++ = (int) *ps++; }
            free(imp);
          } else if (ana_type == 3) {
            float *impf = (float *) imp;
            int *pi = impi;
            while (npix--) { *pi++ = (int) *impf++; }
            free(imp);
          } else {
            fprintf(stderr,"fits file is not F*4, I*2, or I*4\n");  exit(1);
          }
          ana_type = 2;
        } else impi = (int *) imp;
        status = aia_despike(impi, NULL, nxfits, nyfits, 0.8, 7, 3,
                 NULL, 0, &nspikes, &oldvals, &spikelocs, &newvals);
        imp = (float *) impi; /* ana_type preserves intness */
        if (verbosity > 0) printf("nspikes: %d\n", nspikes);
        // status = wrfits2d(imp, "despiked.fits", nxfits, nyfits, ana_type);
      }
      /* limbcompute requires F*4 input so convert if necessary */
      if (ana_type != 3) {
	/* convert to F*4 */
	int	npix = nxfits * nyfits;
	float *impf = (float *) malloc(npix * sizeof(float));
	if (npix <= 0) {
          fprintf(stderr, "empty fits file? # pixels = %d\n", npix); exit(1); }
	if (ana_type == 1) {
	  short *ps = (short *) imp;
	  float *pf = impf;
	  while (npix--) { *pf++ = (float) *ps++; }
	  free(imp);
	  imp = impf;
	} else if (ana_type == 2) {
	  int *pi = (int *) imp;
	  float *pf = impf;
	  while (npix--) { *pf++ = (float) *pi++; }
	  free(imp);
	  imp = impf;
	} else {
	  fprintf(stderr,"fits file is not F*4, I*2, or I*4\n");  exit(1);
	}
      }
      limbcompute(imp, nxfits, nyfits, xcguess, ycguess, rguess, rrange, limbmode, useprevflag, fwhm);
      /* note that result is saved in globals sdisk_xc, sdisk_yc, sdisk_r */
      fprintf(outfile, "%s, %9.2f, %9.2f, %9.2f\n", sq, sdisk_xc, sdisk_yc, sdisk_r);
      ifit++;
    } else ireject++;
    free(sq);
    free(imp);
  } else {
    if (verbosity>2) fprintf(stdout,"list case\n");
    sq = (char *) malloc(strlen(basePath) + strlen(listname)  + 10);
    if (strlen(basePath) >0) {
      strcpy(sq, basePath);  strcat(sq, "/"); strcat(sq,listname);
    } else { strcpy(sq, listname); }
    listfile = fopen(sq, "r");  
    if (listfile == NULL) {
      fprintf(stderr,"could not open list file %s\n", sq);
      exit(1);
    } else {
      char buf[256];
      int ns, status;
      free(sq);
      while( fgets(buf, 255, listfile) != NULL ) {
	char *buf2 = buf;
	/* skip any white space */
	while (isspace( (int) *buf2)) buf2++;
	if (strlen(buf2) <= 0) continue;
	if( strcmp(buf2,"\n") == 0 ) continue;
	if( buf2[0] != '#' ) {	/* if not a comment line */
          /* we expect a file name and, perhaps, guesses for xc, yc, r */
          ffname = (char *) malloc(256);
	  ns = sscanf(buf2, "%250s", ffname);
	  if (ns < 1) { fprintf(stderr,"error decoding a line of list file, line: %s\n", buf);  exit(1); }
	  //if (ns != 3) { fprintf(stderr, "error decoding a line of list file, line: %s\n", buf);  exit(1); }
	  sq = (char *) malloc(strlen(basePath) + strlen(ffname)  + 10);
	  if (strlen(basePath) >0) {
	    strcpy(sq, basePath);  strcat(sq, "/"); strcat(sq,ffname);
	  } else { strcpy(sq, ffname); }
	  /* add a fits to the end if not there already */
	  if (strstr(sq, ".fits") == NULL) strcat(sq, ".fits");
	  if (verbosity>1) fprintf(stdout,"final name %s\n", sq);
//	  status = rdfits2d((void **) &imp, sq, &nxfits, &nyfits, &ana_type, &exp, &wave, &gt_y, &gt_z);
	  if (status) {
            fprintf(stderr, "error reading fits file %s\n", ffname); exit(1); }
          if (verbosity>2) fprintf(stdout,"wave, exp = %f, %f\n", wave, exp);
          if (auto_guess) {
            /* SUVI pix = 180.0e-9 * nano_radians / pi * 3600 / 2.5  */
            gt_y_px = gt_y * 8.25059225e-5;
            gt_z_px = gt_z * 8.25059225e-5;
            xcguess = gt_y_px + 641.0;
            ycguess = gt_z_px + 641.0;
          }
    /* check if we meet minimum exposure time */
          if (exp >= min_exp || exp < 0) {
            if (despike) {
              int *impi, *oldvals, *spikelocs, *newvals, nspikes;
              /* aia_despike requires I*4 input for now */
              if (ana_type != 2) {
                int   npix = nxfits * nyfits;
                impi = (int *) malloc(npix * sizeof(int));
                if (npix <= 0) {
                  fprintf(stderr, "empty fits file? # pixels = %d\n", npix);
                  exit(1);
                }
                if (ana_type == 1) {
                  short *ps = (short *) imp;
                  int *pi = impi;
                  while (npix--) { *pi++ = (int) *ps++; }
                  free(imp);
                } else if (ana_type == 3) {
                  float *impf = (float *) imp;
                  int *pi = impi;
                  while (npix--) { *pi++ = (int) *impf++; }
                  free(imp);
                } else {
                  fprintf(stderr,"fits file is not F*4, I*2, or I*4\n");
                  exit(1);
                }
                ana_type = 2;
              } else impi = (int *) imp;
              status = aia_despike(impi, NULL, nxfits, nyfits, 0.8, 7, 3,
                       NULL, 0, &nspikes, &oldvals, &spikelocs, &newvals);
              imp = (float *) impi; /* ana_type preserves intness */
              if (verbosity > 0) printf("nspikes: %d\n", nspikes);
           // status = wrfits2d(imp, "despiked.fits", nxfits, nyfits, ana_type);
            }       
	    if (ana_type != 3) {
	      /* convert to F*4 */
	      int	npix = nxfits * nyfits;
	      float *impf = (float *) malloc(npix * sizeof(float));
	      if (npix <= 0) { fprintf(stderr, "empty fits file? # pixels = %d\n",  npix);  exit(1); }
	      if (ana_type == 1) {
		short *ps = (short *) imp;
		float *pf = impf;
		while (npix--) { *pf++ = (float) *ps++; }
		free(imp);
		imp = impf;
	      } else if (ana_type == 2) {
		int *pi = (int *) imp;
		float *pf = impf;
		while (npix--) { *pf++ = (float) *pi++; }
		free(imp);
		imp = impf;
	      } else {
		fprintf(stderr,"fits file is not F*4, I*2, or I*4\n"); exit(1);
	      }
	    }
	    limbcompute(imp, nxfits, nyfits, xcguess, ycguess, rguess, rrange, limbmode, useprevflag, fwhm);
	    /* note that result is saved in globals sdisk_xc, sdisk_yc, sdisk_r */
	    fprintf(outfile, "%s, %9.2f, %9.2f, %9.2f\n", sq, sdisk_xc, sdisk_yc, sdisk_r);
	    ifit++;
	    } else ireject++;
	  free(sq);
	  free(imp);
        }
    
      }
    }
  }
  tend = systime();
  if (verbosity>0) fprintf(stdout,"total time for %d fits was %7.2f, %d were rejected (exp < %f)\n", ifit, tend-tstart, ireject, min_exp);
  exit(0);
  }
 
