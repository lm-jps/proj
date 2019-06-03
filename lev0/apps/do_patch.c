// Camera 1 crop table corruption Dec 2011 - Jan 2012
// skip value for row 1065 in quadrant E (lower left)
// was 512 pixels too large.  Need to shift whole row
// right by 512 pixels and fill in MISSING on the left

void do_patch1(short *adata0)
{
    short *p = adata0+4096*1065;
    short tmp[2048];
    int i;
    for (i=0; i<512; ++i) tmp[i] = -32768;
    for (i=512; i<2048; ++i) tmp[i] = p[i-512];
    memcpy(p, tmp, 4096);
}


// Camera 1 lookup table corruption 2014.03.30 from 12:20 UT to 16:25 UT.
// Single bit hit to output value; 0x1912 (6418) changed to 0x3912 (14610)
// FSN range is 70397305 to 70405129
// Any value of 14610 should be changed to 6418

void do_patch2(short *adata0)
{
    int i;
    for (i=0; i<4096*4096; ++i)
	if (adata0[i] == 14610) adata0[i] = 6418;
}


// Camera 2 crop table corruption Dec 31 2016
// take value for row 3989 in quadrant G (upper right)
// was 512 pixels too few, causing all subsequent pixels
// (including those in quadrant H) to be shifted

void do_patch3(short *adata0)
{
    int i,j,idx;
    FILE *fp;
    static short skp[2048], tak[2048];
    static pix[6866600];
    static pix2[6866600];
    static int first_call = 1;

    if (first_call) {
	first_call = 0;
	fp = fopen("/home/production/cvs/TBL_JSOC/lev0/crop/hmi/crop7", "r");
	fscanf(fp, "%d", &i);
	fscanf(fp, "%d", &i);
	fscanf(fp, "%d", &i);
	for (i=0; i<2048; ++i)
	    fscanf(fp, "%hd %hd", &skp[i], &tak[i]);
	fclose(fp);
    }

    // crop upper half of 2D image to 1D
    idx = 0;
    for (i=4095; i>=2048; --i)
	for (j=4095-skp[4095-i]; j>=2048; --j)
	    pix[idx++] = adata0[4096*i+j];
    for (i=4095; i>=2048; --i)
	for (j=skp[4095-i]; j<2048; ++j)
	    pix[idx++] = adata0[4096*i+j];

    // insert 512 BLANKs after row 106 (counting from top) of quadrant G
    for (idx=0; idx<68736; ++idx)
	pix2[idx] = pix[idx];
    for (idx=68736; idx<68736+512; ++idx)
	pix2[idx] = -32768;
    for (idx=68736+512; idx<3433300; ++idx)
	pix2[idx] = pix[idx-512];

    // skip 940 pixels (two rows) at beginning of quadrant H
    for (idx=3433300+940; idx<6866600; ++idx)
	pix2[idx] = pix[idx-940];

    // uncrop upper half of image
    idx = 0;
    for (i=4095; i>=2048; --i)
	for (j=4095-skp[4095-i]; j>=2048; --j)
	    adata0[4096*i+j] = pix2[idx++];
    // skip 940 pixels
    idx += 940;
    // shift quadrant up by two rows
    for (i=4093; i>=2048; --i) {
	for (j=skp[4095-i]; j<2048; ++j)
	    adata0[4096*(i+2)+j] = pix2[idx++];
	// need to fix the left edge or else we will get negative MISSVALS
	for (j=skp[4095-i]; j<skp[4093-i]; ++j)
	    adata0[4096*(i+2)+j] = -32768;
    }
    // fill last two rows of quadrant H with BLANKs
    for (i=2048; i<2050; ++i)
	for (j=0; j<2048; ++j)
	    adata0[4096*i+j] = -32768;
}




// Camera 2 crop table corruption 2017.06.12
// skip value for row 2354 in quadrant G (upper right)
// was 16 pixels too large.  Need to shift whole row
// left by 16 pixels and fill in MISSING on the right

void do_patch4(short *adata0)
{
    int i;
    short *p = adata0+4096*2354;
    for (i=2048; i<4080; ++i)
	p[i] = p[i+16];
    for (i=4080; i<4096; ++i)
	p[i] = -32768;
}




// Camera 2 lookup table corruption 2019.05.30 09:21 - 16:13
// FSN range 157165643 - 157177859
// Any value of 0 should be changed to 5795

void do_patch5(short *adata0)
{
    int i;
    for (i=0; i<4096*4096; ++i)
	if (adata0[i] == 0) adata0[i] = 5795;
}
