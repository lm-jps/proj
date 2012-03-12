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
