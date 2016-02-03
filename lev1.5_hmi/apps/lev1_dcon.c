#include <jsoc_main.h>
#include <mkl.h>
#include <fresize.h>
#include <gapfill.h>
#include "limb_fit.h"

void warn(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
}

void die(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    exit(1);
}

void ck(MKL_LONG status)
{
    if (DftiErrorClass(status, DFTI_NO_ERROR))
	return;
    die(DftiErrorMessage(status));
}

//CORRECTION OF HEIGHT FORMATION (copied from ingest_dcon.c)
//returns 0 if corrections were successful, 1 otherwise
int heightformation(int FID, double OBSVR, float *CDELT1, float *RSUN, float *CRPIX1, float *CRPIX2, float CROTA2)
{
  int wl=0;
  int status=0;
  float correction=0.0,correction2=0.0;
  
  wl = (FID/10)%20;  //temp is now the filter index
  
  if( (wl >= 0) && (wl < 20) )
    {
      correction  = 0.445*exp(-(wl-10.-(float)OBSVR/(0.690/6173.*3.e8/20.)-0.25)*(wl-10.-(float)OBSVR/(0.690/6173.*3.e8/20.)-0.25)/7.1);
      correction2 = 0.39*(-2.0*(wl-10.- (float)OBSVR/(0.690/6173.*3.e8/20.)-0.35)/6.15)*exp(-(wl-10.-(float)OBSVR/(0.690/6173.*3.e8/20.)-0.35)*(wl-10.-(float)OBSVR/(0.690/6173.*3.e8/20.)-0.35)/6.15);
      //printf("A CDELT1= %f CROTA2= %f OBS_VR= %f RSUN= %f correction= %f \n",*CDELT1,CROTA2,OBSVR,*RSUN,correction);
      *CDELT1 = *CDELT1*(*RSUN)/((*RSUN)-correction);
      //printf("B CDELT1= %f CROTA2= %f OBS_VR= %f \n",*CDELT1,CROTA2,OBSVR);
      *RSUN   = *RSUN-correction;
      *CRPIX1 = *CRPIX1-cos(M_PI-CROTA2*M_PI/180.)*correction2;
      *CRPIX2 = *CRPIX2-sin(M_PI-CROTA2*M_PI/180.)*correction2;
	}
  else status=1;

  return status;
}

char *module_name = "lev1_dcon";
ModuleArgs_t module_args[] = {
    {ARG_STRING, "in", "", "input query"},
    {ARG_STRING, "out", "", "output series"},
    {ARG_STRING, "psf", "", "PSF FITS file path"},
    {ARG_INT, "iter", "25", "number of R-L iterations"},
    {ARG_END}
};

int DoIt() {
    const char *in = cmdparams_get_str(&cmdparams, "in", NULL);
    const char *out = cmdparams_get_str(&cmdparams, "out", NULL);
    const char *psf = cmdparams_get_str(&cmdparams, "psf", NULL);
    int iter = cmdparams_get_int(&cmdparams, "iter", NULL);
    int status;

    drms_series_exists(drms_env, out, &status);
    if (DRMS_SUCCESS != status)
	    die("Output series %s does not exit\n", out);

    // open and stage input records
    DRMS_RecordSet_t *rsin = drms_open_records(drms_env, in, &status);
    if (status || !rsin)
	    die("Can't do drms_open_records(%s)\n", in);
    else
	    warn("DRMS recordset %s opened\n", in);
    if (drms_stage_records(rsin, 1, 0) != DRMS_SUCCESS)
    {
        drms_close_records(rsin, DRMS_FREE_RECORD);
        die("Failure staging DRMS records (%s).\n", in);
    }

    int nrecs = rsin->n;
    if (!nrecs)
	    die("No records matching %s in found\n", in);
    else
	    warn("%d records found\n", nrecs);

    // read PSF
    DRMS_Array_t *arrpsf = drms_fitsrw_read(drms_env, psf, 0, NULL, &status);
    if (status || !arrpsf)
	    die("Can't read PSF file %s\n", psf);
    else if (arrpsf->naxis != 2 || arrpsf->axis[0] != 4096 || arrpsf->axis[1] != 4096 || arrpsf->type != DRMS_TYPE_DOUBLE)
	    die("PSF file %s does not contain 4096x4096 double precision array\n", psf);
    else
	    warn("PSF file %s read\n", psf);

    // pad PSF data array for FFT
    float *P = MKL_malloc(sizeof(float)*4096*4098, 64);
    if (!P)
	    die("Can't allocate memory for array P\n");
    for (int j=0; j<4096; ++j)
	    for (int i=0; i<4096; ++i)
	        P[4098*j+i] = ((double *)arrpsf->data)[4096*j+i];
    drms_free_array(arrpsf);

    // set up MKL DFTI stuff
    MKL_LONG lengths[2], strdfor[3], strdbak[3];
    DFTI_DESCRIPTOR *p, *q;
    lengths[0] = lengths[1] = 4096;
    strdfor[0] = 0; strdfor[1] = 4098; strdfor[2] = 1;
    strdbak[0] = 0; strdbak[1] = 2049; strdbak[2] = 1;
    ck(DftiCreateDescriptor(&p, DFTI_SINGLE, DFTI_REAL, 2, lengths));
    ck(DftiSetValue(p, DFTI_BACKWARD_SCALE, 1.0/(4096*4096)));
    ck(DftiSetValue(p, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX));
    ck(DftiSetValue(p, DFTI_INPUT_STRIDES, strdfor));
    ck(DftiSetValue(p, DFTI_OUTPUT_STRIDES, strdbak));
    ck(DftiCommitDescriptor(p));
    ck(DftiCreateDescriptor(&q, DFTI_SINGLE, DFTI_REAL, 2, lengths));
    ck(DftiSetValue(q, DFTI_BACKWARD_SCALE, 1.0/(4096*4096)));
    ck(DftiSetValue(q, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX));
    ck(DftiSetValue(q, DFTI_INPUT_STRIDES, strdbak));
    ck(DftiSetValue(q, DFTI_OUTPUT_STRIDES, strdfor));
    ck(DftiCommitDescriptor(q));

    // forward transform P
    ck(DftiComputeForward(p, P));

    // prepare other work arrays
    // O - original image (with all MISSING values replaced by zeors)
    // E - current estimate
    // Esav - copy of E
    // B - Blurred estimate
    // R - image divided by blurred estimate (shares same location with B)
    // C - correction for next iteraion (shares same location with B)
    float *O = MKL_malloc(sizeof(float)*4096*4098, 64);
    float *E = MKL_malloc(sizeof(float)*4096*4098, 64);
    float *Esav = MKL_malloc(sizeof(float)*4096*4098, 64);
    float *R = MKL_malloc(sizeof(float)*4096*4098, 64);
    float *B = R;
    float *C = R;
    if (!O || !E || !Esav || !R)
        die("Can't allocate memory for array O, E, Esav, or R\n");

    // create all output records (will fail for large nrecs)
    DRMS_RecordSet_t *rsout = drms_create_records(drms_env, nrecs, out, DRMS_PERMANENT, &status);
    if (status || !rsout)
        die("Can't create %d records in output series %s", nrecs, out);

    //
    // MAIN LOOP
    //
    for (int n=0; n<nrecs; ++n) {
        warn("Processing record #%d...\n", n);

        DRMS_Record_t *recin = rsin->records[n];
        DRMS_Record_t *recout = rsout->records[n];
        drms_copykeys(recout, recin, 1, kDRMS_KeyClass_Explicit);

        // copy segment bad_pixel_list
        DRMS_Segment_t *seg1in = drms_segment_lookup(recin, "bad_pixel_list");
        if (!seg1in)
        {
            warn("No bad_pixel_list input segment found.\n");
        } else {
	    DRMS_Array_t *arr1in = drms_segment_read(seg1in, DRMS_TYPE_INT, &status);
	    if (status || !arr1in)
		die("Can't read bad_pixel_list\n");
	    DRMS_Segment_t *seg1out = drms_segment_lookupnum(recout, 1);
	    if (!seg1out)
	    {
		die("Can't open bad_pixel_list output segment.\n");
	    }
	    if (DRMS_SUCCESS != drms_segment_write(seg1out, arr1in, 0))
	    {
		die("Can't write bad_pixel_list output segment.\n");
	    }
	    drms_free_array(arr1in);
	}

        // read segment image_lev1
        DRMS_Segment_t *seg0in = drms_segment_lookup(recin, "image_lev1");
        if (!seg0in)
        {
            warn("No image_lev1 input segment found.\n");
	    continue;
        }
        DRMS_Array_t *arr0in = drms_segment_read(seg0in, DRMS_TYPE_INT, &status);
        if (status || !arr0in)
            die("Can't read image_lev1\n");

        // pad image data array for FFT
        for (int j=0; j<4096; ++j)
            for (int i=0; i<4096; ++i) {
                int tmp = ((int *)arr0in->data)[4096*j+i];
                E[4098*j+i] = O[4098*j+i] = (tmp > 0) ? tmp : 0;
            }

        // do Richardson-Lucy iterations
        for (int k=0; k<iter; ++k) {
            memcpy(Esav, E, sizeof(float)*4096*4098);
            // B = conv(E,P)
            ck(DftiComputeForward(p, E));
            vcMul(2049*4096, E, P, B);
            ck(DftiComputeBackward(q, B));
            // R = O/B
            vsDiv(2*2049*4096, O, B, R);
            // C = corr(R,P) (or conv(R,P*))
            ck(DftiComputeForward(p, R));
            vcMulByConj(2049*4096, R, P, C);
            ck(DftiComputeBackward(q, C));
            // apply correction to get new estimate
            vsMul(2*2049*4096, Esav, C, E);
        }

        // copy final result, write segment
        int *imgout = malloc(sizeof(int)*4096*4096);
        if (!imgout)
            die("Can't allocate memory for output image\n");
        for (int j=0; j<4096; ++j)
            for (int i=0; i<4096; ++i) {
                int tmp = ((int *)arr0in->data)[4096*j+i];
                Esav[4096*j+i] = imgout[4096*j+i] = (tmp > 0) ? roundf(E[4098*j+i]) : tmp;
            }
        DRMS_Segment_t *seg0out = drms_segment_lookupnum(recout, 0);
        if (!seg0out)
        {
            die("Can't open image_lev1 output segment.\n");
        }
        if (arr0in->data) free(arr0in->data);
        arr0in->data = imgout;
        if (DRMS_SUCCESS != drms_segment_write(seg0out, arr0in, 0))
        {
            die("Can't write image_lev1 output segment.\n");
        }
        drms_free_array(arr0in);

        // do limb fit and height formation correction
        int FID = drms_getkey_int(recin, "FID", &status);
        float CROTA2 = drms_getkey_float(recin, "CROTA2", &status);
        float CDELT1 = drms_getkey_float(recin, "CDELT1", &status);
        double OBSVR = drms_getkey_double(recin, "OBS_VR", &status);
        double tmpX0, tmpY0, tmpRSUN;
        status = limb_fit(recin, Esav, &tmpRSUN, &tmpX0, &tmpY0, 4096, 4096, 0);
        float RSUN_LF = tmpRSUN;
        float X0_LF = tmpX0;
        float Y0_LF = tmpY0;
        status = heightformation(FID,OBSVR,&CDELT1,&RSUN_LF,&X0_LF,&Y0_LF,-CROTA2);
        drms_setkey_float(recout, "CRPIX1", X0_LF+1);
        drms_setkey_float(recout, "CRPIX2", Y0_LF+1);
        drms_setkey_float(recout, "R_SUN", RSUN_LF);
        drms_setkey_float(recout, "CDELT1", CDELT1);

	// set CALVER32 bit
	long long calver32 = drms_getkey_longlong(recin, "CALVER32", &status);
	calver32 += 0x200000;
	drms_setkey_longlong(recout, "CALVER32", calver32);
    }

    MKL_free(P);
    MKL_free(O);
    MKL_free(E);
    MKL_free(Esav);
    MKL_free(R);
    drms_close_records(rsin, DRMS_FREE_RECORD);
    drms_close_records(rsout, DRMS_INSERT_RECORD);

    return 0;
}
