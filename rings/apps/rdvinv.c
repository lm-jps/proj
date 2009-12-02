/*****************************************************************
* Module - rdvinv.c
* Responsible - Charles Baldner, Rick Bogart, Sarbani Basu
*
* Velocity inversion code for ring diagram power spectra fits.

* Usage:
*   Input is specified by argument 'ds' (a drms record query) or 
*   by argument 'in', which is a file name, and overrides ds argument 
*   if present.  Output is specified by argment 'out'.  If 'out' is 
*   a drms record set that can be opened, the inversions will be 
*   saved to the drms record set; otherwise, 'out' will be used 
*   as the base filename for the inversion outputs.
*
* Undocumented features (bugs):
*   If the input data is not in a drms series, an output series 
*   can't be created - rdvinv will exit with an error message.
*
*   Arguments 'car', 'lat', and 'lon' are not currently used.
*
*   Argument 'out', when a file name, cannot contain any path 
*   information.
*
*   When writing output to a file, rather than a drms record set, 
*   only one set of files can be specified - if the input contains 
*   more than one data set, each inversion will overwrite the last.
*
******************************************************************/


#include <stdio.h>
#include <string.h>
#include <math.h>
#include <jsoc_main.h>
#include "keystuff.c"

/* prototypes */
int read_fit(FILE *fpt, int *npts, int **n, double **l, double **f,
    double **ef, double **ux, double **eux, double **uy, double **euy);
int interp(int *n, double *l, double *f, double *ef, double *ux,
    double *eux, double *uy, double *euy, int npts);
int autoweed_vel(int* n, double* l, double *ux, double *uy,
    int *mask, int npts);
int nearst(double xb, double x[], int ntab);
int divdif(double xb, double x[], double f[], int *nuse, int ntab,
    double fb[], double aeps, double *dfb, double *ddfb);
extern void ola_(double *, int *, double *, double *, double *,
      double *, double *, int *, char *, char *, char *, char *, char *,
      int *, int *, int *, double *, double *,
      int *, double *, double *, double *,
      int, int, int, int, int);

char *module_name = "rdvinv";
char *version_id = "0.2";

ModuleArgs_t module_args[] = {
  {ARG_STRING, "ds", "su_rsb.rdfits", "data set"},
  {ARG_STRING, "seg", "fit.out", "data segment name"},
  {ARG_INT, "car", "Not Specified", "Carrington rotation"},
  {ARG_DOUBLE, "lat", "Not Specified", "AR latitude"},
  {ARG_DOUBLE, "lon", "Not Specified", "AR longitude"},
  {ARG_DOUBLE, "amu", "0.005", "Error trade-off parameter"},
  {ARG_DOUBLE, "ob", "1.0", "Lower frequency limit (mHz)"},
  {ARG_DOUBLE, "oe", "5.2", "Upper frequency limit (mHz)"},
  {ARG_DOUBLE, "rb", "0.97", "Lower radius limit"},
  {ARG_DOUBLE, "re", "1.00", "Upper radius limit"},
  {ARG_INT, "num", "40", "Number of target inversion points"},
  {ARG_STRING, "in", "Not Specified", "Input filename - if set, overides drms lookup"},
  {ARG_STRING, "out", "fort.10.hmixy", "Output filename"},
  {ARG_STRING, "kernel", "/tmp21/baldner/kernels/ringker.kur_cm_opal_78_mod.ascii", ""},
  {ARG_STRING, "ave", "Not Specified", "output file for averaging kernels (if not set, kernels are not written)"},
  {ARG_STRING, "coef", "Not Specified", ""},
  {ARG_FLAG,    "v",    "", "run in verbose mode"},
  {ARG_END}
};

int read_fit(FILE *fpt, int *npts, int **n, double **l, double **f,
    double **ef, double **ux, double **eux, double **uy, double **euy)       {
  /*
     Function - read_fit
     Takes an open file pointer, and reads n, l, frequencies, and velocities
  */
  int i, nlines;
  char buffer[8192];
  if( ferror(fpt) )     {
    //perror(fpt);
    return -1;
  } 
  if( feof(fpt)) rewind(fpt);
  
  nlines = 0;
  while(!feof(fpt))     {
    fgets(buffer, 8192, fpt);
    if(buffer[0] != '#' && !feof(fpt)) nlines++;
  }
  if(nlines == 0) return 1;
  
  *n = (int*) malloc(nlines*sizeof(int));
  *l = (double*) malloc(nlines*sizeof(double));
  *f = (double*) malloc(nlines*sizeof(double));
  *ef = (double*) malloc(nlines*sizeof(double));
  *ux = (double*) malloc(nlines*sizeof(double));
  *eux = (double*) malloc(nlines*sizeof(double));
  *uy = (double*) malloc(nlines*sizeof(double));
  *euy = (double*) malloc(nlines*sizeof(double));
  
  rewind(fpt);
  
  for(i=0; i<nlines; i++)       {
    fgets(buffer, 8192, fpt);
    while(buffer[0] == '#') fgets(buffer, 8192, fpt);
    sscanf(buffer, "%i %lf %*f %lf %lf %lf %lf %lf %lf", &(*n)[i], &(*l)[i], &(*f)[i], 
        &(*ef)[i], &(*ux)[i], &(*eux)[i], &(*uy)[i], &(*euy)[i]);
  } 
  *npts = nlines;
  
  return 0;
}

int interp(int *n, double *l, double *f, double *ef, double *ux, 
      double *eux, double *uy, double *euy, int npts)	{
  /*
     Function - interp
     Takes power spectra fit parameters, and interpolates them to integer l.
     Data are returned in the input arrays - **data is overwritten!**
  */
  int i, j, nuse, ierr, offset=0;
  int n_num[13];
  double fb[10], ll;
  double dfb, ddfb, aeps=1e-7;
  double *inp_l, *inp_f, *inp_ef, *inp_ux, *inp_uy, *inp_eux, *inp_euy;
  FILE * outBug;
  
  // count number of frequencies for each n
  // note that it is assumed that all frequencies for each 
  //	n are grouped together in the arrays
  for(i=0; i<13; i++) n_num[i]=0;
  for(i=0; i < npts; i++)	{
	  n_num[n[i]] ++;
  }
  
  inp_l=(double*) malloc(npts*sizeof(double));
  inp_f=(double*) malloc(npts*sizeof(double));
  inp_ef=(double*) malloc(npts*sizeof(double));
  inp_ux=(double*) malloc(npts*sizeof(double));
  inp_uy=(double*) malloc(npts*sizeof(double));
  inp_eux=(double*) malloc(npts*sizeof(double));
  inp_euy=(double*) malloc(npts*sizeof(double));
  // loop over n
//  outBug=fopen("tesbug","w");
  for(i=0; i < 10; i++)	{
    for(j=0;j<npts;j++)	{
      inp_l[j]=0.0;
      inp_f[j]=0.0;
      inp_ef[j]=0.0;
      inp_ux[j]=0.0;
      inp_uy[j]=0.0;
      inp_eux[j]=0.0;
      inp_euy[j]=0.0;
    }
    for(j=0;j<n_num[i];j++)	{
      inp_l[j]=l[j+offset];
      inp_f[j]=f[j+offset];
      inp_ef[j]=ef[j+offset];
      inp_ux[j]=ux[j+offset];
      inp_uy[j]=uy[j+offset];
      inp_eux[j]=eux[j+offset];
      inp_euy[j]=euy[j+offset];
    }
  for(j=0;j<n_num[i];j++)	{
    ll=(int)(inp_l[j]+0.5);
    nuse=3;
    ierr=divdif(ll, inp_l, inp_f, &nuse, 
	n_num[i], fb, aeps, &dfb, &ddfb);
    f[j+offset]=fb[nuse];
    nuse=3;
    ierr=divdif(ll, inp_l, inp_ef, &nuse,
	n_num[i], fb, aeps, &dfb, &ddfb);
    ef[j+offset]=fb[nuse];
    nuse=3;
    ierr=divdif(ll, inp_l, inp_ux, &nuse, 
	n_num[i], fb, aeps, &dfb, &ddfb);
    ux[j+offset]=fb[nuse];
    nuse=3;
    ierr=divdif(ll, inp_l, inp_uy, &nuse, 
	n_num[i], fb, aeps, &dfb, &ddfb);
    uy[j+offset]=fb[nuse];
    nuse=3;
    ierr=divdif(ll, inp_l, inp_eux, &nuse, 
	n_num[i], fb, aeps, &dfb, &ddfb);
    eux[j+offset]=fb[nuse];
    nuse=3;
    ierr=divdif(ll, inp_l, inp_euy, &nuse, 
	n_num[i], fb, aeps, &dfb, &ddfb);
    euy[j+offset]=fb[nuse];
    l[j+offset]=ll;
      }
    offset+=n_num[i];
  }
  
  return 0;
}

int divdif(double xb, double x[], double f[], int *nuse, int ntab,
    double fb[], double aeps, double *dfb, double *ddfb)	{
  int i,j,k,next,in,ip,nit,ier, nmax=10;
  double err,px,dpx,ddpx,xn[11],xd[11];

/*	Find the nearest point */

  next=nearst(xb,x,ntab);
  fb[1]=f[next];
  xd[1]=f[next];
  xn[1]=x[next];
  ier=0;
  px=1.0;

/*	Initialisation for the derivatives */

  *dfb=0.0; *ddfb=0.0;
  dpx=0.0; ddpx=0.0;

/*	Points between IN and IP are used for interpolation */

  ip=next; in=next;

/*	Maximum number of points to be used for interpolation */
  nit=*nuse; if(nmax<nit) nit=nmax; if(ntab<nit) nit=ntab;
  if(*nuse>nmax || *nuse>ntab) ier=22;
  if(*nuse<1) {
    ier=21;
    nit=6; if(nmax<nit) nit=nmax; if(ntab<nit) nit=ntab;
  }
  *nuse=1;
	  
/*	Calculate successive interpolation polynomial */
  for(j=2; j<=nit; ++j) {
/*	Choose the next nearest point to XB */
    if(in<=0 ) {
      ip=ip+1; next=ip;
    }
    else if(ip >= ntab-1) {
      in=in-1; next=in;
    }
    else if(fabs(xb-x[ip+1]) < fabs(xb-x[in-1]) ) {
      ip=ip+1; next=ip;
    }
    else {
      in=in-1; next=in;
    }

/*	Calculating the divided differences */
    xd[j]=f[next];
    xn[j]=x[next];
    for(k=j-1; k>=1; --k) xd[k]=(xd[k+1]-xd[k])/(xn[j]-xn[k]);

/*	Calculating the derivatives */
    ddpx=ddpx*(xb-xn[j-1])+2.*dpx;
    dpx=dpx*(xb-xn[j-1])+px;
    *dfb = *dfb+dpx*xd[1];
    *ddfb = *ddfb+ddpx*xd[1];

    px=px*(xb-xn[j-1]);
    err=xd[1]*px;
    fb[j]=fb[j-1]+err;
    *nuse=j;

    if(fabs(err) < aeps) return ier;
  }
  return 23;
}

int nearst(double xb, double x[], int ntab)	{
  int low, igh, mid;

  low=0; igh=ntab-1;
  if((xb < x[low]) != (xb < x[igh]) ) {

/*	If the point is within the range of table, then locate it by bisection */

    while(igh-low > 1) {
      mid=(low+igh)/2;
      if((xb < x[mid]) == (xb < x[low])) low=mid;
      else igh=mid;
    }
  }

  if(fabs(xb-x[low]) < fabs(xb-x[igh])) return low;
  else return igh;
}

int autoweed_vel(int* n, double* l, double *ux, double *uy, 
    int *mask, int npts)	{
  /*
     Function - autoweed_vel
     Function to automatically weed velocities.  An int array 
     'mask' is returned of length npts - mask[i] is False for 
     rejected modes.  Ux and Uy are weeded together.
  */
  const int maxn=8;
  const double tol=5.0;

  int i, j, offset;
  int n_num[maxn];
  double num;
  double sumx, sumxx, sumy, sumyy, meanx, meany, stdx, stdy;
  double range, uxmin, uxmax, uymin, uymax;
  double lmin,lmax;
	  
  for(i=0; i<maxn; i++)	{
    n_num[i]=0;
  }
  for(i=0; i<npts; i++)	{
    if(n[i] < maxn) n_num[n[i]]++;
  }

  offset=0;
  for(i=0; i<maxn; i++)	{
    num=0.0;
    sumx=0.0;
    sumy=0.0;
    sumxx=0.0;
    sumyy=0.0;
    // get l limits
    lmin=lmax=l[offset];
    for(j=offset; j<n_num[i]+offset; j++)	{
      if(l[j] < lmin) lmin=l[j];
      if(l[j] > lmax) lmax= l[j];
    }
    range=(lmax-lmin)/5.0;
    //printf("%f %f %f",lmin,lmax,range);
    lmin+=2.0*range;
    lmax-=2.0*range;
    for(j=offset; j<n_num[i]+offset; j++)	{
      if(l[j] <= lmax && l[j] >= lmin)	{
	sumx+=ux[j];
	sumxx+=ux[j]*ux[j];
	sumy+=uy[j];
	sumyy+=uy[j]*uy[j];
	num++;
      }
    }
    meanx=sumx/num;
    meany=sumy/num;
    stdx=num*sumxx-sumx*sumx;
    stdy=num*sumyy-sumy*sumy;
    stdx/=(num*(num-1));
    stdy/=(num*(num-1));
    stdx=sqrt(stdx);
    stdy=sqrt(stdy);
    //printf("%f  %f  %f  %f\n",stdx,stdy,lmin,lmax);
    uxmin=meanx-5.0*stdx;
    uxmax=meanx+5.0*stdx;
    uymin=meany-5.0*stdy;
    uymax=meany+5.0*stdy;
    //printf("%f  %f  %f  %f\n",uxmin,uxmax,uymin,uymax);
    
    for(j=offset; j<n_num[i]+offset; j++)	{
      if(ux[j] <= uxmax && ux[j] >= uxmin &&
	  uy[j] <= uymax && uy[j] >= uymin) mask[j]=1;
      else mask[j]=0;
    }
    offset+=n_num[i];
  }
  return 0;
}

/* main loop */
int DoIt(void)	{
  CmdParams_t *params = &cmdparams;
  int status = 0;
  int verbose = 0;
  int drms_output;
  int rec_ct, rec_i, seg_ct;
  DRMS_RecordSet_t *recordSet = NULL;
  DRMS_Record_t *record, *orecord;
  DRMS_Segment_t *segment, *osegment_ux, *osegment_uy;
  DRMS_Array_t *keytable_coord;
  char *ds = params_get_str (params, "ds");
  char *seg = params_get_str (&cmdparams, "seg");
  char *kernel = params_get_str(&cmdparams, "kernel");
  char filename[DRMS_MAXPATHLEN+5], *recnum_querry;
  char buffer[1024], outfilex[DRMS_MAXPATHLEN], outfiley[DRMS_MAXPATHLEN], *suffix;
  FILE *infile;
  char *keylist_coord = "recnum"; //"LatHG,LonHG,CarrTime,recnum";
  int car = params_get_int(&cmdparams, "car");
  double lat = params_get_double(&cmdparams, "lat");
  double lon = params_get_double(&cmdparams, "lon");
  int *open_record;
  double *data_ptr;
  int *n, *mask, npts, i, j;
  double *l, *f, *ef, *ux, *eux, *uy, *euy;
  int qave, qcoef, lenkern, lenoutx, lenouty, lenave, lencoef, lends, filename_set;
  double ob = params_get_double(&cmdparams, "ob");
  double oe = params_get_double(&cmdparams, "oe");
  double rb = params_get_double(&cmdparams, "rb");
  double re = params_get_double(&cmdparams, "re");
  double amu = params_get_double(&cmdparams, "amu");
  int num = params_get_int(&cmdparams, "num");
  char *out = params_get_str(&cmdparams, "out");
  char *ave = params_get_str(&cmdparams, "ave");
  char *coef = params_get_str(&cmdparams, "coef");


  if(params_isflagset (&cmdparams, "v")) verbose = 1;
  if(strcmp(ave, "Not Specified")==0)	{
    qave = 0;
    lenave = 0;
  } else	{
    qave = 1;
    lenave = strlen(ave);
  }
  if(strcmp(coef, "Not Specified")==0)	{
    qcoef = 0;
    lencoef = 0;
  } else	{
    qcoef = 1;
    lencoef = strlen(coef);
  }
  lenkern = strlen(kernel);

  //Attempt to create output drms record
  orecord = drms_create_record (drms_env, out, DRMS_PERMANENT, &status);
  if(orecord)	{	// if record creation was successful
    seg_ct = drms_record_numsegments(orecord);
    if(seg_ct < 1)	{
      fprintf(stderr, "Error: no data segment in output series %s\n", out);
      drms_close_record(orecord, DRMS_FREE_RECORD);
      return 1;
    }
    osegment_ux = drms_segment_lookup(orecord, "Ux_inv.out");
    osegment_uy = drms_segment_lookup(orecord, "Uy_inv.out");
    // check that the appropriate data segments exist and are of type generic
    if(!osegment_ux || !osegment_uy)	{
      fprintf(stderr, "Could not find 2 output segments in output series %s\n", out);
      drms_close_record(orecord, DRMS_FREE_RECORD);
      return 1;
    }
    if(osegment_ux->info->protocol != DRMS_GENERIC)	{
      fprintf(stderr, "Ux segment is not of type GENERIC in output series %s,\n", out);
      drms_close_record(orecord, DRMS_FREE_RECORD);
      return 1;
    }
    if(osegment_uy->info->protocol != DRMS_GENERIC)	{
      fprintf(stderr, "Uy segment is not of type GENERIC in output series %s,\n", out);
      drms_close_record(orecord, DRMS_FREE_RECORD);
      return 1;
    }
    drms_output = 1;
    fprintf(stderr, "Using drms series for output\n");
  }
  else	{	// if record creation was not successful, assume `out' is a filename
    sprintf(outfilex, "%s.ux", out);
    sprintf(outfiley, "%s.uy", out);
    lenoutx = strlen(outfilex);
    lenouty = strlen(outfiley);
    drms_output = 0;
    fprintf(stderr, "Not using drms series for output\n");
  }
  strcpy(filename, params_get_str(&cmdparams, "in"));

  // If argument 'in' is specified, attempt to open that filename, ignore drms
  if(strcmp(filename, "Not Specified")!=0)	{
    if(drms_output)	{
      drms_close_record(orecord, DRMS_FREE_RECORD);
      fprintf(stderr, "Error:  rdvinv cannot write to a drms series if input is not from a series\n");
      return 1;
    }
    filename_set = 1;
    infile = fopen(filename, "r");
    read_fit(infile, &npts, &n, &l, &f, &ef, &ux, &eux, &uy, &euy);
    fclose(infile);
    mask = (int*) malloc(npts*sizeof(int));
    interp(n, l, f, ef, ux, eux, uy, euy, npts);
    autoweed_vel(n, l, ux, uy, mask, npts);
    j = 0;
    for(i=0; i<npts; i++)	{
      if(mask[i])	{
	n[j] = n[i];
	l[j] = l[j];
	f[j] = f[i];
	ef[j] = ef[i];
	ux[j] = ux[i];
	eux[j] = eux[j];
	uy[j] = uy[i];
	euy[j] = euy[i];
	j++;
      }
    }
    npts = j;
    if(drms_output)	{
      drms_segment_filename(osegment_ux, outfilex);
      lenoutx = strlen(outfilex);
      drms_segment_filename(osegment_uy, outfiley);
      lenouty = strlen(outfiley);
    }
    ola_(l, n, f, ux, eux, uy, euy, &npts,
	kernel, outfilex, outfiley, ave, coef, &qave, &qcoef, &verbose, &ob, &oe, &num,
	&rb, &re, &amu, lenkern, lenoutx, lenouty, lenave, lencoef);
    free(n);
    free(mask);
    free(l);
    free(f);
    free(ef);
    free(ux);
    free(eux);
    free(uy);
    free(euy);
  }
  else	{
				  /* get a table of keyvalues for the search string ds */
    // Input is from 'ds' query
    /* Currently, the input loop is a bit strange:  the query is passed to 
       drms_record_getvector(), which passes back a table of keywords and 
       values, which will allow using arguments car, lat, lon to select 
       subsets from the query.  I may remove this.
    */
    keytable_coord = drms_record_getvector(drms_env, ds, keylist_coord, 
	DRMS_TYPE_DOUBLE, 0, &status);
    if(status != DRMS_SUCCESS)	{
      printf("Call to drms_record_getvector() failed, status code %i\n", status);
      return 1;
    }
    rec_ct = keytable_coord->axis[1];
    lends = strlen(ds);
    i=0;
    while(ds[i]!='[') i++;
    if(i<lends-1) {
      ds[i] = '\0';
      lends = i+1;
    }
    recnum_querry = (char*) malloc((lends+21)*sizeof(char));
    data_ptr = (double*) keytable_coord->data;
    for(rec_i=0; rec_i<rec_ct; rec_i++) {
      sprintf(recnum_querry, "%s[!recnum=%i!]", ds, (int) data_ptr[rec_i]);
      recordSet = drms_open_records(drms_env, recnum_querry, &status);
      record = recordSet->records[0];
      if (record->sunum != -1LL && record->su == NULL)	{
	record->su = drms_getunit (record->env, record->seriesinfo->seriesname,
	    record->sunum, 1, &status);
      }
      segment = drms_segment_lookup(record, seg);
      sprintf (filename, "%s/S%05i/%s", segment->record->su->sudir, segment->record->slotnum, seg);

      infile = fopen(filename, "r");
      status = read_fit(infile, &npts, &n, &l, &f, &ef, &ux, &eux, &uy, &euy);
      if( status )	{
	fprintf(stderr, "Error - no lines read from file %s\n", filename);
	drms_close_record(record, DRMS_FREE_RECORD);
	drms_close_records(recordSet, DRMS_FREE_RECORD);
	return 1;
      }
      fclose(infile);

      mask = (int*) malloc(npts*sizeof(int));
      interp(n, l, f, ef, ux, eux, uy, euy, npts);
      autoweed_vel(n, l, ux, uy, mask, npts);
      j = 0;
      for(i = 0; i < npts; i ++)	{
	if(mask[i])	{
	  n[j] = n[i];
	  l[j] = l[i];
	  f[j] = f[i];
	  ux[j] = ux[i];
	  uy[j] = uy[i];
	  eux[j] = eux[i];
	  euy[j] = euy[i];
	  j++;
	}
      }
      npts = j;
      if(drms_output)	{
	if(rec_i>0)	{
	  orecord = drms_create_record (drms_env, out, DRMS_PERMANENT, &status);
	  osegment_ux = drms_segment_lookup(orecord, "Ux_inv.out");
	  osegment_uy = drms_segment_lookup(orecord, "Uy_inv.out");
	}
	drms_segment_filename(osegment_ux, outfilex);
	suffix = strstr (outfilex, ".generic");
	if(suffix) *suffix = '\0';
	lenoutx = strlen(outfilex);
	drms_segment_filename(osegment_uy, outfiley);
	suffix = strstr (outfiley, ".generic");
	if(suffix) *suffix = '\0';
	lenouty = strlen(outfiley);
      }

      fprintf(stderr, "File Ux: %s  %i\n", outfilex, lenoutx);
      fprintf(stderr, "File Uy: %s  %i\n", outfiley, lenouty);

      ola_(l, n, f, ux, eux, uy, euy, &npts, 
	  kernel, outfilex, outfiley, ave, coef, &qave, &qcoef, &verbose, &ob, &oe, &num,
	  &rb, &re, &amu, lenkern, lenoutx, lenouty, lenave, lencoef);
      free(n);
      free(mask);
      free(l);
      free(f);
      free(ef);
      free(ux);
      free(eux);
      free(uy);
      free(euy);
      		/* If a drms output segment was opened, write the inversions to it */
      if(drms_output)	{
	int kstat = 0;
	char *key_str;
	double key_dbl;
	int key_int;
	TIME key_time;
	key_str = drms_getkey_string (record, "PrimeKeyString", &status);
	if(!status && check_and_set_key_str (orecord, "PrimeKeyString", key_str)) kstat = 1;
	key_str = drms_getkey_string (record, "CarrTime", &status);
	if(!status && check_and_set_key_str (orecord, "CarrTime", key_str)) kstat = 1;
	key_dbl = drms_getkey_double (record, "LonHG", &status);
	if(!status && check_and_set_key_float (orecord, "LonHG", key_dbl)) kstat = 1;
	key_dbl = drms_getkey_double (record, "LatHG", &status);
	if(!status && check_and_set_key_float (orecord, "LatHG", key_dbl)) kstat = 1;
	if(check_and_set_key_float(orecord, "amu", amu)) kstat = 1;
	if(check_and_set_key_float(orecord, "freqmin", ob)) kstat = 1;
	if(check_and_set_key_float(orecord, "freqmax", oe)) kstat = 1;
	if(check_and_set_key_float(orecord, "radmin", rb)) kstat = 1;
	if(check_and_set_key_float(orecord, "radmax", re)) kstat = 1;
	key_time = drms_getkey_time (record, "MidTime", &status);
	if(!status && check_and_set_key_time (orecord, "MidTime", key_time)) kstat = 1;
	key_int = drms_getkey_int (record, "Duration", &status);
	if(!status && check_and_set_key_int (orecord, "Duration", key_int)) kstat = 1;
	key_str = drms_getkey_string (record, "MapProj", &status);
	if(!status && check_and_set_key_str (orecord, "MapProj", key_str)) kstat = 1;
	key_dbl = drms_getkey_double (record, "MapScale", &status);
	if(!status && check_and_set_key_float (orecord, "MapScale", key_dbl)) kstat = 1;
	key_dbl = drms_getkey_double (record, "Map_PA", &status);
	if(!status && check_and_set_key_float (orecord, "Map_PA", key_dbl)) kstat = 1;
	key_dbl = drms_getkey_double (record, "Width", &status);
	if(!status && check_and_set_key_float (orecord, "Width", key_dbl)) kstat = 1;
	key_dbl = drms_getkey_double (record, "Height", &status);
	if(!status && check_and_set_key_float (orecord, "Height", key_dbl)) kstat = 1;
	key_dbl = drms_getkey_double (record, "ZonalTrk", &status);
	if(!status && check_and_set_key_float (orecord, "ZonalTrk", key_dbl)) kstat = 1;
	key_dbl = drms_getkey_double (record, "MeridTrk", &status);
	if(!status && check_and_set_key_float (orecord, "MeridTrk", key_dbl)) kstat = 1;
	if(check_and_set_key_str(orecord, "Module", module_name)) kstat = 1;
	key_str = drms_getkey_string(record, "Module", &status);
	if(check_and_set_key_str(orecord, "Source", key_str)) kstat = 1;
	drms_close_record(orecord, DRMS_INSERT_RECORD);
      }
      drms_close_record(record, DRMS_FREE_RECORD);
      drms_close_records(recordSet, DRMS_FREE_RECORD);
    }
  }

  return 0;
}

/*
   Revision History
   June 2009 - v 0.1 module created, no real drms functionality
   09.10.2009 - v 0.2 drms functionality added - input and output from drms records now 
   	possible, but input/output to rooted filenames also possible
*/
