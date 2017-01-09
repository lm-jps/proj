// crop table corruption Dec 2011 - Jan 2012
// skip value for row 1065 in quadrant E (lower left)
// was 512 pixels too large.  Need to shift whole row
// right by 512 pixels and fill in MISSING on the left

void do_patch1(short *adata0)
{
    short *p = adata0+4096*1065;
    short tmp[2048];
    int i;
    for (i=0; i<512; ++i) tmp[i] = DRMS_MISSING_SHORT;
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
