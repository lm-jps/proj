#define ASDSERIES "sdo.lev0_asd_0004"

#define PT_OUT_OF_MEMORY (-1)
#define PT_DRMS_OPEN_FAILED (-2)

typedef struct {
    float sat_y0;
    float sat_z0;
    float sat_rot;
    char acs_mode[16];
    char acs_eclp[16];
    char acs_sunp[16];
    char acs_safe[16];
    char acs_cgt[16];
    char asd_rec[64];
} PTINFO;

static int find_closest(TIME *t, int n, TIME tt)
{
    TIME d, dmin = DBL_MAX;
    int i, idx = -1;

    if (!_DRMS_IS_T_MISSING(tt))
	for (i=0; i<n; ++i) {
	    if (_DRMS_IS_T_MISSING(t[i])) continue;
	    d = fabs(tt - t[i]);
	    if (d < dmin) {
		dmin = d; idx = i;
	    }
	}

    return idx;
}

// for a fixed vector x0 in some reference frame
// find its components x in another reference
// frame rotated from the first one by q
//
//      x = inv(q) x0 q
//
static void qrot(const double *q, const double *x0, double *x)
{
    double q11,q12,q13,q21,q22,q23,q31,q32,q33;
    q11 = 2 * (q[0]*q[0] + q[1]*q[1]) - 1.0;
    q12 = 2 * (q[1]*q[2] + q[0]*q[3]);
    q13 = 2 * (q[1]*q[3] - q[0]*q[2]);
    q21 = 2 * (q[1]*q[2] - q[0]*q[3]);
    q22 = 2 * (q[0]*q[0] + q[2]*q[2]) - 1.0;
    q23 = 2 * (q[2]*q[3] + q[0]*q[1]);
    q31 = 2 * (q[1]*q[3] + q[0]*q[2]);
    q32 = 2 * (q[2]*q[3] - q[0]*q[1]);
    q33 = 2 * (q[0]*q[0] + q[3]*q[3]) - 1.0;
    x[0] = q11*x0[0] + q12*x0[1] + q13*x0[2];
    x[1] = q21*x0[0] + q22*x0[1] + q23*x0[2];
    x[2] = q31*x0[0] + q32*x0[1] + q33*x0[2];
}

// quaternion product
//
//      r = p q
//
static void qmul(const double *p, const double *q, double *r)
{
    r[0] = p[0]*q[0] - p[1]*q[1] - p[2]*q[2] - p[3]*q[3];
    r[1] = p[0]*q[1] + p[1]*q[0] + p[2]*q[3] - p[3]*q[2];
    r[2] = p[0]*q[2] + p[2]*q[0] + p[3]*q[1] - p[1]*q[3];
    r[3] = p[0]*q[3] + p[3]*q[0] + p[1]*q[2] - p[2]*q[1];
}

// quaternion inverse
static void qinv(double *q)
{
    q[1] = -q[1]; q[2] = -q[2]; q[3] = -q[3];
}

// quaternion norm
static double qnorm(double *q)
{
    return sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
}

int get_pointing_info(DRMS_Env_t *drms_env, TIME *tobs, int nobs, PTINFO **ptinfo)
{
    PTINFO *p;
    DRMS_RecordSet_t *rset;
    TIME tstart, tend, t, *tasd;
    int *asd;
    char dsname[256];
    int status, nrec, i, j;

    asd = (int *) malloc(nobs*sizeof(int));
    p = *ptinfo = (PTINFO *) malloc(nobs*sizeof(PTINFO));
    if (!p) 
	return PT_OUT_OF_MEMORY;

    for (i=0; i<nobs; ++i) {
	asd[i] = -1;
	p[i].sat_y0 = p[i].sat_z0 = p[i].sat_rot = DRMS_MISSING_FLOAT;
	p[i].acs_eclp[0] = p[i].acs_sunp[0] = p[i].acs_safe[0] = 0;
	p[i].acs_mode[0] = p[i].acs_cgt[0] = p[i].asd_rec[0] = 0;
    }

    tstart = tobs[0] - 30.0;
    tend = tobs[nobs-1] + 30.0;
    sprintf(dsname, "%s[? PACKET_TIME > %.0f and PACKET_TIME < %.0f ?]", ASDSERIES, tstart, tend);
    rset = drms_open_records(drms_env, dsname, &status);
    if (!rset || !rset->n || status) return PT_DRMS_OPEN_FAILED;
    nrec = rset->n;

    tasd = (TIME *) malloc(nrec*sizeof(TIME));
    if (!tasd) return PT_OUT_OF_MEMORY;
    for (i=0; i<nrec; ++i)
	tasd[i] = drms_getkey_time(rset->records[i], "PACKET_TIME", &status);

    // find which ASD packet to use based on closest match between tobs and tasd
    for (i=0; i<nobs; ++i) {
	int done = 0;
	int idx;
	char *str;
	double qbcs[4], qsn[4], qq[4];
	double x[3], x0[3] = {1.0,0.0,0.0};
	double z[3], z0[3] = {0.0,0.0,1.0};
	double tmp;

	asd[i] = idx = find_closest(tasd, nrec, tobs[i]);

	for (j=0; j<i; ++j)
	    if (asd[i] == asd[j]) {	// have used this packet before
		memcpy(&p[i], &p[j], sizeof(PTINFO));
		done = 1;
		break;
	    }
	if (!done) {
	    if (asd[i] < 0)		// no match for whatever reason; give up
		continue;

	    str = drms_getkey_string(rset->records[idx], "ACS_AN_FLAG_CSS_ECLIPSE", &status);
	    if (str) {
		strncpy(p[i].acs_eclp, str, 16);
		free(str);
	    }
	    str = drms_getkey_string(rset->records[idx], "ACS_AN_FLAG_ACE_INSAFEHOLD", &status);
	    if (str) {
		strncpy(p[i].acs_safe, str, 16);
		free(str);
	    }
	    str = drms_getkey_string(rset->records[idx], "ACS_AN_FLAG_DSS_SUNPRES", &status);
	    if (str) {
		strncpy(p[i].acs_sunp, str, 16);
		free(str);
	    }
	    str = drms_getkey_string(rset->records[idx], "ACS_AN_NUM_CGT", &status);
	    if (str) {
		strncpy(p[i].acs_cgt, str, 16);
		free(str);
	    }
	    str = drms_getkey_string(rset->records[idx], "ACS_AN_ACS_MODE", &status);
	    if (str) {
		strncpy(p[i].acs_mode, str, 16);
		free(str);
	    }
	    snprintf(p[i].asd_rec, 64, "%s[:#%lld]", ASDSERIES, rset->records[idx]->recnum);

	    status = 0;
	    qbcs[0] = drms_getkey_float(rset->records[idx], "ACS_AN_QUAT_GCIFTOBCSF_EST_S", &status);
	    qbcs[1] = drms_getkey_float(rset->records[idx], "ACS_AN_QUAT_GCIFTOBCSF_EST_X", &status);
	    qbcs[2] = drms_getkey_float(rset->records[idx], "ACS_AN_QUAT_GCIFTOBCSF_EST_Y", &status);
	    qbcs[3] = drms_getkey_float(rset->records[idx], "ACS_AN_QUAT_GCIFTOBCSF_EST_Z", &status);
	    qsn[0] = drms_getkey_float(rset->records[idx], "ACS_AN_QUAT_GCIFTOSNRF_S", &status);
	    qsn[1] = drms_getkey_float(rset->records[idx], "ACS_AN_QUAT_GCIFTOSNRF_X", &status);
	    qsn[2] = drms_getkey_float(rset->records[idx], "ACS_AN_QUAT_GCIFTOSNRF_Y", &status);
	    qsn[3] = drms_getkey_float(rset->records[idx], "ACS_AN_QUAT_GCIFTOSNRF_Z", &status);
	    if (status) continue;

	    tmp = qnorm(qbcs);
	    if (tmp > 1.001 || tmp < 0.999) continue;
	    tmp = qnorm(qsn);
	    if (tmp > 1.001 || tmp < 0.999) continue;

	    qinv(qsn);
	    qmul(qsn, qbcs, qq);
	    
	    // x0: sun-pointing vector in SNR frame
	    // x:  sun-pointing vector in BCS frame
	    qrot(qq, x0, x);

	    // z0: solar north in SNR frame
	    // z:  solar north in BCS frame
	    qrot(qq,z0,z);

	    p[i].sat_rot = atan2(z[1],z[2])*180.0/M_PI;
	    p[i].sat_y0 = asin(x[1])*3600.0*180.0/M_PI;
	    p[i].sat_z0 = asin(x[2])*3600.0*180.0/M_PI;
	}
    }

    drms_close_records(rset, DRMS_FREE_RECORD);
    free(tasd);
    free(asd);
    return 0;
}
