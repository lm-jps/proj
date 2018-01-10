#include <jsoc_main.h>
#include <mkl.h>

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

#define DORL(X)\
for (int k=0; k<iter; ++k) {\
    cblas_scopy(4096*4098, E##X, 1, Esav##X, 1);\
    ck(DftiComputeForward(p, E##X));\
    vcMul(2049*4096, E##X, P, B##X);\
    ck(DftiComputeBackward(q, B##X));\
    vsDiv(2*2049*4096, O##X, B##X, R##X);\
    ck(DftiComputeForward(p, R##X));\
    vcMulByConj(2049*4096, R##X, P, C##X);\
    ck(DftiComputeBackward(q, C##X));\
    vsMul(2*2049*4096, Esav##X, C##X, E##X);\
}

char *module_name = "stokes_dcon";
ModuleArgs_t module_args[] = {
    {ARG_STRING, "in", "", "input query"},
    {ARG_STRING, "out", "", "output series"},
    {ARG_STRING, "psf", "", "PSF query"},
    {ARG_INT, "iter", "25", "number of R-L iterations"},
    {ARG_END}
};

int DoIt() {
    const char *in = cmdparams_get_str(&cmdparams, "in", NULL);
    const char *out = cmdparams_get_str(&cmdparams, "out", NULL);
    const char *psf = cmdparams_get_str(&cmdparams, "psf", NULL);
    char source[100];
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
    DRMS_RecordSet_t *rpsf = drms_open_records(drms_env, psf, &status);
    if (status || !rpsf)
	die("Can't do drms_open_records(%s)\n", psf);
    DRMS_Segment_t *segpsf = drms_segment_lookup(rpsf->records[0], "psf");
    if (!segpsf)
	die("No psf segment found.\n");
    DRMS_Array_t *arrpsf = drms_segment_read(segpsf, DRMS_TYPE_DOUBLE, &status);
    if (status || !arrpsf)
	die("Can't read psf record %s\n", psf);
    else if (arrpsf->naxis != 2 || arrpsf->axis[0] != 4096 || arrpsf->axis[1] != 4096 || arrpsf->type != DRMS_TYPE_DOUBLE)
	die("PSF record %s does not contain 4096x4096 double precision array\n", psf);
    warn("PSF record %s read\n", psf);

    // pad PSF data array for FFT
    float *P = MKL_malloc(sizeof(float)*4096*4098, 64);
    if (!P)
	    die("Can't allocate memory for array P\n");
    for (int j=0; j<4096; ++j)
	    for (int i=0; i<4096; ++i)
	        P[4098*j+i] = ((double *)arrpsf->data)[4096*j+i];
    drms_free_array(arrpsf);
    drms_close_records(rpsf, DRMS_FREE_RECORD);

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
    float *Oipq = MKL_malloc(sizeof(float)*4096*4098, 64);	// I+Q
    float *Oimq = MKL_malloc(sizeof(float)*4096*4098, 64);	// I-Q
    float *Oipu = MKL_malloc(sizeof(float)*4096*4098, 64);	// I+U
    float *Oimu = MKL_malloc(sizeof(float)*4096*4098, 64);	// I-U
    float *Oipv = MKL_malloc(sizeof(float)*4096*4098, 64);	// I+V
    float *Oimv = MKL_malloc(sizeof(float)*4096*4098, 64);	// I-V
    float *Eipq = MKL_malloc(sizeof(float)*4096*4098, 64);	// I+Q
    float *Eimq = MKL_malloc(sizeof(float)*4096*4098, 64);	// I-Q
    float *Eipu = MKL_malloc(sizeof(float)*4096*4098, 64);	// I+U
    float *Eimu = MKL_malloc(sizeof(float)*4096*4098, 64);	// I-U
    float *Eipv = MKL_malloc(sizeof(float)*4096*4098, 64);	// I+V
    float *Eimv = MKL_malloc(sizeof(float)*4096*4098, 64);	// I-V
    float *Esavipq = MKL_malloc(sizeof(float)*4096*4098, 64);	// I+Q
    float *Esavimq = MKL_malloc(sizeof(float)*4096*4098, 64);	// I-Q
    float *Esavipu = MKL_malloc(sizeof(float)*4096*4098, 64);	// I+U
    float *Esavimu = MKL_malloc(sizeof(float)*4096*4098, 64);	// I-U
    float *Esavipv = MKL_malloc(sizeof(float)*4096*4098, 64);	// I+V
    float *Esavimv = MKL_malloc(sizeof(float)*4096*4098, 64);	// I-V
    float *Ripq = MKL_malloc(sizeof(float)*4096*4098, 64);	// I+Q
    float *Rimq = MKL_malloc(sizeof(float)*4096*4098, 64);	// I-Q
    float *Ripu = MKL_malloc(sizeof(float)*4096*4098, 64);	// I+U
    float *Rimu = MKL_malloc(sizeof(float)*4096*4098, 64);	// I-U
    float *Ripv = MKL_malloc(sizeof(float)*4096*4098, 64);	// I+V
    float *Rimv = MKL_malloc(sizeof(float)*4096*4098, 64);	// I-V
    float *Bipq = Ripq, *Cipq = Ripq;
    float *Bimq = Rimq, *Cimq = Rimq;
    float *Bipu = Ripu, *Cipu = Ripu;
    float *Bimu = Rimu, *Cimu = Rimu;
    float *Bipv = Ripv, *Cipv = Ripv;
    float *Bimv = Rimv, *Cimv = Rimv;
    if (!Rimv)
        die("Can't allocate memory for array O, E, Esav, or R\n");

    // create all output records (will fail for large nrecs)
    DRMS_RecordSet_t *rsout = drms_create_records(drms_env, nrecs, out, DRMS_PERMANENT, &status);
    if (status || !rsout)
        die("Can't create %d records in output series %s", nrecs, out);

    //
    // LOOP OVER RECORDS
    //
    for (int n=0; n<nrecs; ++n) {
	warn("Processing record #%d...\n", n);
	DRMS_Record_t *recin = rsin->records[n];
	DRMS_Record_t *recout = rsout->records[n];
	drms_copykeys(recout, recin, 1, kDRMS_KeyClass_Explicit);

	//
	// LOOP OVER WAVELENGTH POSITIONS
	//
	for (int wl=0; wl<6; ++wl) {
	    DRMS_Segment_t *segIin = drms_segment_lookupnum(recin, 4*wl);
	    DRMS_Segment_t *segQin = drms_segment_lookupnum(recin, 4*wl+1);
	    DRMS_Segment_t *segUin = drms_segment_lookupnum(recin, 4*wl+2);
	    DRMS_Segment_t *segVin = drms_segment_lookupnum(recin, 4*wl+3);

	    DRMS_Array_t *arrI = drms_segment_read(segIin, DRMS_TYPE_FLOAT, &status);
	    if (status || !arrI)
		die("Can't read I at wavelength position %d\n", wl);
	    DRMS_Array_t *arrQ = drms_segment_read(segQin, DRMS_TYPE_FLOAT, &status);
	    if (status || !arrQ)
		die("Can't read Q at wavelength position %d\n", wl);
	    DRMS_Array_t *arrU = drms_segment_read(segUin, DRMS_TYPE_FLOAT, &status);
	    if (status || !arrU)
		die("Can't read U at wavelength position %d\n", wl);
	    DRMS_Array_t *arrV = drms_segment_read(segVin, DRMS_TYPE_FLOAT, &status);
	    if (status || !arrV)
		die("Can't read V at wavelength position %d\n", wl);

#pragma omp parallel for
	    // pad data arrays for FFT
	    for (int j=0; j<4096; ++j) {
		int idx = 4096*j;
		int id2 = 4098*j;
		for (int i=0; i<4096; ++i,++idx,++id2) {
		    float tmpi = ((float *)arrI->data)[idx];
		    float tmpq = ((float *)arrQ->data)[idx];
		    float tmpu = ((float *)arrU->data)[idx];
		    float tmpv = ((float *)arrV->data)[idx];
		    Eipq[id2] = Oipq[id2] = isnanf(tmpi+tmpq) ? 0 : tmpi+tmpq;
		    Eimq[id2] = Oimq[id2] = isnanf(tmpi-tmpq) ? 0 : tmpi-tmpq;
		    Eipu[id2] = Oipu[id2] = isnanf(tmpi+tmpu) ? 0 : tmpi+tmpu;
		    Eimu[id2] = Oimu[id2] = isnanf(tmpi-tmpu) ? 0 : tmpi-tmpu;
		    Eipv[id2] = Oipv[id2] = isnanf(tmpi+tmpv) ? 0 : tmpi+tmpv;
		    Eimv[id2] = Oimv[id2] = isnanf(tmpi-tmpv) ? 0 : tmpi-tmpv;
		}
	    }

	    // do Richardson-Lucy iterations
	    DORL(ipq);
	    DORL(imq);
	    DORL(ipu);
	    DORL(imu);
	    DORL(ipv);
	    DORL(imv);

#pragma omp parallel for
	    // reconstruct I,Q,U,V
	    for (int j=0; j<4096; ++j) {
		int idx = 4096*j;
		int id2 = 4098*j;
		for (int i=0; i<4096; ++i,++idx,++id2) {
		    ((float *)arrI->data)[idx] = (Eipq[id2]+Eimq[id2]+Eipu[id2]+Eimu[id2]+Eipv[id2]+Eimv[id2]) / 6.0;
		    ((float *)arrQ->data)[idx] = 0.5 * (Eipq[id2] - Eimq[id2]);
		    ((float *)arrU->data)[idx] = 0.5 * (Eipu[id2] - Eimu[id2]);
		    ((float *)arrV->data)[idx] = 0.5 * (Eipv[id2] - Eimv[id2]);
		}
	    }

	    // write segments
	    DRMS_Segment_t *segIout = drms_segment_lookupnum(recout, 4*wl);
	    if (DRMS_SUCCESS != drms_segment_write(segIout, arrI, 0))
		die("Can't write I%d segment.\n", wl);
	    drms_free_array(arrI);

	    DRMS_Segment_t *segQout = drms_segment_lookupnum(recout, 4*wl+1);
	    if (DRMS_SUCCESS != drms_segment_write(segQout, arrQ, 0))
		die("Can't write Q%d segment.\n", wl);
	    drms_free_array(arrQ);

	    DRMS_Segment_t *segUout = drms_segment_lookupnum(recout, 4*wl+2);
	    if (DRMS_SUCCESS != drms_segment_write(segUout, arrU, 0))
		die("Can't write U%d segment.\n", wl);
	    drms_free_array(arrU);

	    DRMS_Segment_t *segVout = drms_segment_lookupnum(recout, 4*wl+3);
	    if (DRMS_SUCCESS != drms_segment_write(segVout, arrV, 0))
		die("Can't write V%d segment.\n", wl);
	    drms_free_array(arrV);

	    drms_setkey_double(recout, "DATE", CURRENT_SYSTEM_TIME);
	    //drms_setkey_int(recout, "ITER", iter);
	    //drms_setkey_string(recout, "PSF", psf);
	    //snprintf(source, sizeof source, "%s[:#%lld]", recin->seriesinfo->seriesname, recin->recnum);
	    //drms_setkey_string(recout, "SOURCE", source);
	}
    }

    drms_close_records(rsin, DRMS_FREE_RECORD);
    drms_close_records(rsout, DRMS_INSERT_RECORD);

    return 0;
}
