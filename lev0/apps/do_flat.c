#define XXX_CANT_OPEN_HIMGCFGFILE -1
#define XXX_CORRUPT_HIMGCFGFILE -2
#define XXX_CANT_FIND_HIMGCFID -3

#define HIMGCFGFILE "/home/production/img_cnfg_ids"

#define MAXHIMGCFGS 4096
static int overscan_nrows[MAXHIMGCFGS];
static int overscan_ncols[MAXHIMGCFGS];
static int overscan_rowstr[MAXHIMGCFGS];
static int overscan_colstr[MAXHIMGCFGS];
static int overscan_valid[MAXHIMGCFGS];

///////////////////////////////////
int get_overscan_info(int himgcfid)
{
    static int firstcall = 1;
    static FILE *fp;
    static char line[256];
    int found = 0;
    int id;
    char j[64];

    if (firstcall) {
	fp = fopen(HIMGCFGFILE,"r");
	if (!fp)
	    return XXX_CANT_OPEN_HIMGCFGFILE;
	firstcall = 0;
    }

    rewind(fp);
    while (fgets(line, 256, fp)) {
	if (1 != sscanf(line, "%d", &id))
	    continue;
	if (id != himgcfid)
	    continue;
	if (13 != sscanf(line, "%s%s%s%s%s%s%s%s%s%d%d%d%d", j,j,j,j,j,j,j,j,j,
		&overscan_nrows[id], &overscan_ncols[id],
		&overscan_rowstr[id], &overscan_colstr[id]))
	    return XXX_CORRUPT_HIMGCFGFILE;
	found = 1;
	overscan_valid[himgcfid] = 1;
    }

    if (!found)
	return XXX_CANT_FIND_HIMGCFID;

    return 0;
}

///////////////////////////
int do_flat(LEV0LEV1 *info)
{
    int i,j,idx,IDX;
    int nr,nc,r1,r2,c1,c2;
    int status;
    static float dat[MAXPIXELS];
    float *fdat = info->adataff;
    float *ddat = info->adatadark;
    short *in = info->adata0;
    short *out = info->adata1;
    short tmp;

    if (!overscan_valid[info->himgcfid]) {
	status = get_overscan_info(info->himgcfid);
	if (status) return status;
    }

    nr = overscan_nrows[info->himgcfid];
    if (nr) {
	r1 = overscan_rowstr[info->himgcfid];
	r2 = r1 + nr;
    } else
	r1 = r2 = 4096;
    nr /= 2;

    nc = overscan_ncols[info->himgcfid];
    if (nc) {
	c1 = overscan_colstr[info->himgcfid];
	c2 = c1 + nc;
    } else
	c1 = c2 = 4096;
    nc /= 2;

    for (i=0; i<r1; ++i) {
	idx = 4096*i;
	IDX = idx + 4096*nr + nc;
	for (j=0; j<c1; ++j) {
	    tmp = in[idx++];
	    if (BLANK == tmp)
		out[IDX++] = BLANK;
	    else {
		out[IDX] = (tmp - ddat[IDX]) * fdat[IDX];
		++IDX;
	    }
	}
    }
    for (i=0; i<r1; ++i) {
	idx = 4096*i + c2;
	IDX = idx + 4096*nr - nc;
	for (j=c2; j<4096; ++j) {
	    tmp = in[idx++];
	    if (BLANK == tmp)
		out[IDX++] = BLANK;
	    else {
		out[IDX] = (tmp - ddat[IDX]) * fdat[IDX];
		++IDX;
	    }
	}
    }
    for (i=r2; i<4096; ++i) {
	idx = 4096*i;
	IDX = idx - 4096*nr + nc;
	for (j=0; j<c1; ++j) {
	    tmp = in[idx++];
	    if (BLANK == tmp)
		out[IDX++] = BLANK;
	    else {
		out[IDX] = (tmp - ddat[IDX]) * fdat[IDX];
		++IDX;
	    }
	}
    }
    for (i=r2; i<4096; ++i) {
	idx = 4096*i + c2;
	IDX = idx - 4096*nr - nc;
	for (j=c2; j<4096; ++j) {
	    tmp = in[idx++];
	    if (BLANK == tmp)
		out[IDX++] = BLANK;
	    else {
		out[IDX] = (tmp - ddat[IDX]) * fdat[IDX];
		++IDX;
	    }
	}
    }

    if (nr) {
	for (i=0; i<4096*nr; ++i) out[i] = BLANK;
	for (i=4096*(4096-nr); i<MAXPIXELS; ++i) out[i] = BLANK;
    }
    if (nc) {
	for (i=0; i<4096; ++i) {
	    for (j=0; j<nc; ++j)
		out[i*4096+j] = BLANK;
	    for (j=4096-nc; j<4096; ++j)
		out[i*4096+j] = BLANK;
	}
    }

    return 0;
}
