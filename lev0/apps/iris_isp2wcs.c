#include <drms.h>

#define GETKEY_ERROR_1 1
#define GETKEY_ERROR_2 2
#define GETKEY_ERROR_3 3
#define GETKEY_ERROR_4 4
#define GETKEY_ERROR_5 5
#define GETKEY_ERROR_6 6
#define GETKEY_ERROR_7 7
#define GETKEY_ERROR_8 8
#define GETKEY_ERROR_9 9
#define GETKEY_ERROR_10 10
#define GETKEY_ERROR_11 11
#define GETKEY_ERROR_12 12
#define GETKEY_ERROR_13 13
#define GETKEY_ERROR_14 14
#define GETKEY_ERROR_15 15
#define GETKEY_ERROR_16 16
#define GETKEY_ERROR_17 17
#define GETKEY_ERROR_18 18
#define GETKEY_ERROR_19 19
#define GETKEY_ERROR_20 20
#define SETKEY_ERROR_1 -1
#define SETKEY_ERROR_2 -2
#define SETKEY_ERROR_3 -3
#define SETKEY_ERROR_4 -4
#define SETKEY_ERROR_5 -5
#define SETKEY_ERROR_6 -6
#define SETKEY_ERROR_7 -7
#define SETKEY_ERROR_8 -8
#define SETKEY_ERROR_9 -9
#define SETKEY_ERROR_10 -10
#define SETKEY_ERROR_11 -11
#define SETKEY_ERROR_12 -12
#define SETKEY_ERROR_13 -13
#define SETKEY_ERROR_14 -14
#define SETKEY_ERROR_15 -15
#define SETKEY_ERROR_16 -16
#define SETKEY_ERROR_17 -17
#define SETKEY_ERROR_18 -18
#define SETKEY_ERROR_19 -19
#define SETKEY_ERROR_20 -20
#define BAD_INSTRUMENT_NAME -99
#define BAD_IMG_PATH -98

#define RADEG (180.0/M_PI)

int iris_isp2wcs(DRMS_Record_t *rs0, DRMS_Record_t *rs1, double scroll)
{
    int i, status;

// Set database values

int version   = 2;
double
al_fu1            =         -0.05000,
al_fu2            =         -0.43700,
al_nuv            =         -0.31800,
be_133            =         -0.20000,
be_140            =         -0.22400,
be_279            =          0.27400,
be_283            =          0.28600,
be_fsi            =          0.28000,
be_mir            =         -0.21000,
be_fu1            =         -1.5000,
be_fu2            =         -1.9700,
be_nuv            =          1.2600,
cdlt1_f1          =          0.01298,
cdlt1_f2          =          0.01272,
cdlt1_nu          =          0.02546,
cdlt2_f1          =          0.16632,
cdlt2_f2          =          0.16632,
cdlt2_nu          =          0.16637,
cdlt_133          =          0.16560,
cdlt_140          =          0.16560,
cdlt_279          =          0.16790,
cdlt_283          =          0.16790,
cdlt_fsi          =          0.16790,
cdlt_mir          =          0.16560,
cpx1_133          =        537.30,
cpx1_140          =        529.32,
cpx1_279          =        504.69,
cpx1_283          =        506.47,
cpx1_fsi          =        505.20,
cpx1_fu1          =        219.50,
cpx1_fu2          =        3807.3,
cpx1_mir          =        534.00,
cpx1_nuv          =        659.30,
cpx2_133          =        524.03,
cpx2_140          =        510.46,
cpx2_279          =        502.91,
cpx2_283          =        502.22,
cpx2_fsi          =        504.70,
cpx2_fu1          =        487.92,
cpx2_fu2          =        518.21,
cpx2_mir          =        516.70,
cpx2_nuv          =        530.35,
cvl1_fu1          =       1334.53,
cvl1_fu2          =       1402.77,
cvl1_nuv          =       2798.65,
focus_regcal      =       -100.00,
focus_xregr1      =         -0.0514701,
focus_xregr2      =          0.000248873,
focus_yregr1      =          0.0739689,
focus_yregr2      =         -0.00047510,
pzt_off_homea     =       -250.00,
pzt_off_homeb     =       -250.00,
pzt_off_homec     =       -250.00,
pzt_off_pixel     =          0.16790,
pzt_off_sjixa     =          0.001300,
pzt_off_sjixb     =         -0.14929,
pzt_off_sjixc     =          0.14783,
pzt_off_sjiya     =         -0.17203,
pzt_off_sjiyb     =          0.08709,
pzt_off_sjiyc     =          0.08687,
pzt_x_sign        =          1.0000,
pzt_y_sign        =          1.0000,
sji_ccd_roll_bias =         -0.28600,
slit_roll_bias    =          0.64600,
slit_rot          =          0.64600,
wm1_rmax          =        634.80,
wm2_rmax          =        628.50,
wm_focus_cal      =       -100.00,
wm_m1v            =        217.875,
wm_offset2        =        113.57,
wm_roll_bias      =         -0.646,
wm_x_sign         =         -1.000,
wm_xoff           =        103.00,
wm_y_sign         =          1.000,
wm_yoff           =        -92.000;

//Regression coefficients for NUV SJI X & Y motion of laser spot on
//image in pixels vs focus motor position in steps,
double
focus_xregr[] = {focus_xregr1, focus_xregr2},
focus_yregr[] = {focus_yregr1, focus_yregr2},
foc_rgcl = focus_regcal; 

//Derive the X & Y shifts (px) due to focus, for the position wm_focus_cal
double focus_temp[] = {wm_focus_cal-foc_rgcl, pow((wm_focus_cal-foc_rgcl),2)};
double wm_focus_cal_x = focus_xregr[0]*focus_temp[0]+focus_xregr[1]*focus_temp[1];
double wm_focus_cal_y = focus_yregr[0]*focus_temp[0]+focus_yregr[1]*focus_temp[1];

// *****************************************************************

//get relevant keywords from Image Status Packet (image header)

char *instrume, *img_path, *iissloop;
int isqpzta, isqpztb, isqpztc;
short iwbpzta, iwbpztb, iwbpztc, igtpoffx, igtpoffy;
float ifmpos, iwm1cpos, iwm2cpos;
short sum1, sum2; 
float off1, off2;

instrume = drms_getkey_string(rs0, "INSTRUME", &status);
if (status) return GETKEY_ERROR_1;
img_path = drms_getkey_string(rs0, "IMG_PATH", &status);
if (status) return GETKEY_ERROR_2;
iissloop = drms_getkey_string(rs0, "IISSLOOP", &status);
if (status) return GETKEY_ERROR_3;
isqpzta  = drms_getkey_int(rs0, "ISQPZTA", &status);
if (status) return GETKEY_ERROR_4;
isqpztb  = drms_getkey_int(rs0, "ISQPZTB", &status);
if (status) return GETKEY_ERROR_5;
isqpztc  = drms_getkey_int(rs0, "ISQPZTC", &status);
if (status) return GETKEY_ERROR_6;
iwbpzta  = drms_getkey_short(rs0, "IWBPZTA", &status);
if (status) return GETKEY_ERROR_7;
iwbpztb  = drms_getkey_short(rs0, "IWBPZTB", &status);
if (status) return GETKEY_ERROR_8;
iwbpztc  = drms_getkey_short(rs0, "IWBPZTC", &status);
if (status) return GETKEY_ERROR_9;
igtpoffx = drms_getkey_short(rs0, "IGTPOFFX", &status);
if (status) return GETKEY_ERROR_10;
igtpoffy = drms_getkey_short(rs0, "IGTPOFFY", &status);
if (status) return GETKEY_ERROR_11;
ifmpos   = drms_getkey_float(rs0, "IFMPOS", &status);
if (status) return GETKEY_ERROR_12;
iwm1cpos = drms_getkey_float(rs0, "IWM1CPOS", &status);
if (status) return GETKEY_ERROR_13;
iwm2cpos = drms_getkey_float(rs0, "IWM2CPOS", &status);
if (status) return GETKEY_ERROR_14;
sum1 = drms_getkey_short(rs0, "SUMSPTRL", &status);
if (status) return GETKEY_ERROR_15;
sum2 = drms_getkey_short(rs0, "SUMSPAT", &status);
if (status) return GETKEY_ERROR_16;
off1 = (sum1-1.0)/(sum1*2.0);
off2 = (sum2-1.0)/(sum2*2.0);


//  *****************************************************************

// Set spacecraft roll
//if n_elements(scroll) ne 1 then scroll = 0.0
double roll_angle = scroll;
//double garad = (slit_rot - scroll) * RADEG;
double garad = (-1.0*slit_rot - scroll)/RADEG; 

//  *****************************************************************

int pzts[3];
float pzt_off_sjix[3], pzt_off_sjiy[3];

pzts[0] = isqpzta - iwbpzta - pzt_off_homea;
pzts[1] = isqpztb - iwbpztb - pzt_off_homeb;
pzts[2] = isqpztc - iwbpztc - pzt_off_homec;

pzt_off_sjix[0] = pzt_off_sjixa;
pzt_off_sjix[1] = pzt_off_sjixb;
pzt_off_sjix[2] = pzt_off_sjixc;

pzt_off_sjiy[0] = pzt_off_sjiya;
pzt_off_sjiy[1] = pzt_off_sjiyb;
pzt_off_sjiy[2] = pzt_off_sjiyc;

// minus sign because solar coord change = -motion of image
double ccdx=0.0, ccdy=0.0;
for (i=0;i<3;++i) {
    ccdx -= (pzt_off_sjix[i]*pzts[i])*pzt_off_pixel*pzt_x_sign;
    ccdy -= (pzt_off_sjiy[i]*pzts[i])*pzt_off_pixel*pzt_y_sign;
}

// add shift of image by focus mechanism using current position focus
//    and subtract image shift by focus at the calibration pos
//focus_temp = [float(ifmpos)-foc_rgcl, (float(ifmpos)-foc_rgcl)^2]
focus_temp[0] = ifmpos - foc_rgcl;
focus_temp[1] = focus_temp[0] * focus_temp[0];
double focusx, focusy;
focusx = focus_temp[0]*focus_xregr[0]+focus_temp[1]*focus_xregr[1] - wm_focus_cal_x;
focusy = focus_temp[0]*focus_yregr[0]+focus_temp[1]*focus_yregr[1] - wm_focus_cal_y;

// minus sign because solar coord change = -motion of image
double ccdx2, ccdy2;
ccdx2 = ccdx - focusx*pzt_off_pixel;
ccdy2 = ccdy - focusy*pzt_off_pixel;

// rotate these to solar X, Y at roll = 0
double sth, cth;
sth = sin((slit_roll_bias+sji_ccd_roll_bias)/RADEG);
cth = cos((slit_roll_bias+sji_ccd_roll_bias)/RADEG);
ccdx =  ccdx2*cth + ccdy2*sth;
ccdy = -ccdx2*sth + ccdy2*cth;

// add wedge motor image shifts from formulae in iris_wedge2solar_tl
double fac = 2.*M_PI/240.;
double th1 =  fac*(iwm1cpos-wm_m1v);
double th2 =  fac*(iwm2cpos-wm_m1v+120.-wm_offset2);
double xwm2 = wm_x_sign*(wm1_rmax*sin(th1) + wm2_rmax*sin(th2));
double ywm2 = wm_y_sign*(wm1_rmax*cos(th1) + wm2_rmax*cos(th2));

// rotate by roll bias wrt solar North at roll=0 degrees
sth = sin((slit_roll_bias + wm_roll_bias)/RADEG);
cth = cos((slit_roll_bias + wm_roll_bias)/RADEG);
double xwm =  xwm2*cth + ywm2*sth;
double ywm = -xwm2*sth + ywm2*cth;

// add the center of the WM coordinate system:  note this is in solar
// X & Y, not in the wedge motor or GT diode or CCD X & Y
xwm = xwm + wm_xoff;
ywm = ywm + wm_yoff;

// co-add both contributions
double x2 = xwm + ccdx;
double y2 = ywm + ccdy;

// add GT bias offsets in ACS packets, only if ISS loop is not closed
// if loop is closed, then PZTs will zero this (assuming they have
// enough range)
// this is very rare and may never be used during science observing,
// so ignore the small roll bias between the GT diodes coordinate
// system and the others

if (strncmp(iissloop,"CLOSED",6)) {
// minus sign because solar coord change = -motion of image
    x2 = x2 - igtpoffx*0.01;
    y2 = y2 - igtpoffy*0.01;
}

// Now rotate everything to non-zero roll angle
double sroll = sin(roll_angle/RADEG);
double croll = cos(roll_angle/RADEG);
double x =  x2*croll + y2*sroll;
double y = -x2*sroll + y2*croll;

double crvalxy[] = {x,y};

//print, 'Calculated x and y', crvalxy

//  *****************************************************************

// set keywords: separate cases for FUV, NUV, SJI

float cdelt,alrad,berad,wcsaxes,crpix1,crpix2,crval1,crval2,crval3,cdelt1,cdelt2,cdelt3,berad1,alrad1,ab1cos,pc1_1,pc1_2,pc2_1,pc2_2,pc3_1,pc3_2,crpix1a,crpix2a,crval1a,crval2a,crval3a,cdelt1a,cdelt2a,cdelt3a,berad2,alrad2,ab2cos,pc1_1a,pc1_2a,pc2_1a,pc2_2a,pc3_1a,pc3_2a,lonpole,crpix3,crpix3a,xcen,ycen;
char *ctype1,*ctype2,*ctype3,*cunit1,*cunit2,*cunit3;
char *ctype1a,*ctype2a,*ctype3a,*cunit1a,*cunit2a,*cunit3a;
char *specsys1;
DRMS_Segment_t *segment0;
int naxis1, naxis2;
segment0 = drms_segment_lookupnum(rs0, 0);
naxis1 = segment0->axis[0];
naxis2 = segment0->axis[1];

if (strncmp(instrume,"FUV",3) == 0) {
    wcsaxes = 3;
    crpix1  = cpx1_fu1/sum1 + off1;
    crpix2  = cpx2_fu1/sum2 + off2;
    //crpix3 = 1.0;
    crval1  = cvl1_fu1;
    crval2  = crvalxy[1];
    crval3  = crvalxy[0];
    cdelt1  = cdlt1_f1*sum1;
    cdelt2  = cdlt2_f1*sum2;
    cdelt3  = cdelt2;
    ctype1  = "WAVE";
    ctype2  = "HPLT-TAN";
    ctype3  = "HPLN-TAN";
    cunit1  = "Angstrom";
    cunit2  = "arcsec";
    cunit3  = cunit2;
    berad1 = be_fu1 / RADEG;
    alrad1 = al_fu1 / RADEG;
    ab1cos  =  cos(alrad1-berad1);
    pc1_1   =  cos(berad1) / ab1cos;
    pc1_2   = -sin(berad1) / ab1cos * sum2/sum1;
    pc2_1   =  cos(garad) * sin(alrad1) / ab1cos * sum1/sum2;
    pc2_2   =  cos(garad) * cos(alrad1) / ab1cos;
    pc3_1   = -sin(garad) * sin(alrad1) / ab1cos * sum1/sum2;
    pc3_2   = -sin(garad) * cos(alrad1) / ab1cos;
    crpix1a = cpx1_fu2/sum1 + off1;
    crpix2a = cpx2_fu2/sum2 + off2;
    //crpix3a = 1.0;
    crval1a = cvl1_fu2;
    crval2a = crval2;
    crval3a = crval3;
    cdelt1a = cdlt1_f2*sum1;
    cdelt2a = cdlt2_f2*sum2;
    cdelt3a = cdelt2a;
    ctype1a = ctype1;
    ctype2a = ctype2;
    ctype3a = ctype3;
    cunit1a = cunit1;
    cunit2a = cunit2;
    cunit3a = cunit3;
    berad2 = be_fu2 / RADEG;
    alrad2 = al_fu2 / RADEG;
    ab2cos  =  cos(alrad2-berad2);
    pc1_1a  =  cos(berad2) / ab2cos;
    pc1_2a  = -sin(berad2) / ab2cos * sum2/sum1;
    pc2_1a  =  cos(garad) * sin(alrad2) / ab2cos * sum1/sum2;
    pc2_2a  =  cos(garad) * cos(alrad2) / ab2cos;
    pc3_1a  = -sin(garad) * sin(alrad2) / ab2cos * sum1/sum2;
    pc3_2a  = -sin(garad) * cos(alrad2) / ab2cos;
    xcen = crval3 + cdelt3*(pc3_1*((naxis1+1)/2. - crpix1) + pc3_2*((naxis2+1)/2. - crpix2));
    ycen = crval2 + cdelt2*(pc2_1*((naxis1+1)/2. - crpix1) + pc2_2*((naxis2+1)/2. - crpix2));
    //lonpole = 180.0;
    //specsys1= "HELIOCEN";
} else if (strncmp(instrume,"NUV", 3)==0) {
    wcsaxes = 3;
    berad = be_nuv / RADEG;
    alrad = al_nuv / RADEG;
    crpix1  = cpx1_nuv/sum1 + off1;
    crpix2  = cpx2_nuv/sum2 + off2;
    //crpix3  = 1.0;
    crval1  = cvl1_nuv;
    crval2  = crvalxy[1];
    crval3  = crvalxy[0];
    cdelt1  = cdlt1_nu*sum1;
    cdelt2  = cdlt2_nu*sum2;
    cdelt3  = cdelt2;
    ctype1  = "WAVE";
    ctype2  = "HPLT-TAN";
    ctype3  = "HPLN-TAN";
    cunit1  = "Angstrom";
    cunit2  = "arcsec";
    cunit3  = cunit2;
    pc1_1   =  cos(berad) / cos(alrad-berad);
    pc1_2   = -sin(berad) / cos(alrad-berad) * sum2/sum1;
    pc2_1   =  cos(garad) * sin(alrad) / cos(alrad-berad) * sum1/sum2;
    pc2_2   =  cos(garad) * cos(alrad) / cos(alrad-berad);
    pc3_1   = -sin(garad) * sin(alrad) / cos(alrad-berad) * sum1/sum2;
    pc3_2   = -sin(garad) * cos(alrad) / cos(alrad-berad);
    crpix1a = crpix1;
    crpix2a = crpix2;
    //crpix3a = 1.0;
    crval1a = crval1;
    crval2a = crval2;
    crval3a = crval3;
    cdelt1a = cdelt1;
    cdelt2a = cdelt2;
    cdelt3a = cdelt3;
    ctype1a = ctype1;
    ctype2a = ctype2;
    ctype3a = ctype3;
    cunit1a = cunit1;
    cunit2a = cunit2;
    cunit3a = cunit3;
    pc1_1a  = pc1_1;
    pc1_2a  = pc1_2;
    pc2_1a  = pc2_1;
    pc2_2a  = pc2_2;
    pc3_1a  = pc3_1;
    pc3_2a  = pc3_2;
    xcen = crval3 + cdelt3*(pc3_1*((naxis1+1)/2. - crpix1) + pc3_2*((naxis2+1)/2. - crpix2));
    ycen = crval2 + cdelt2*(pc2_1*((naxis1+1)/2. - crpix1) + pc2_2*((naxis2+1)/2. - crpix2));
    //lonpole = 180.0;
    //specsys1= "HELIOCEN";    // or 'SOURCE  ' if corrected for solar rot
} else if (strncmp(instrume,"SJI", 3)==0) {
    if(strncmp(img_path,"SJI_5000W",9)==0) {
        // berad = be_nsj / RADEG;
        // cdelt = cdlt_nsj;
        berad = be_fsi / RADEG;
        cdelt = cdlt_fsi;
        crpix1 = cpx1_fsi;
        crpix2 = cpx2_fsi;
    } else if(strncmp(img_path,"SJI_1330",8) == 0) {
        // berad = be_fsj / RADEG;
        // cdelt = cdlt_fsj;
        berad = be_133 / RADEG;
        cdelt = cdlt_133;
        crpix1 = cpx1_133;
        crpix2 = cpx2_133;
    } else if(strncmp(img_path,"SJI_2796",8) == 0) {
        // berad = be_nsj / RADEG;
        // cdelt = cdlt_nsj;
        berad = be_279 / RADEG;
        cdelt = cdlt_279;
        crpix1 = cpx1_279;
        crpix2 = cpx2_279;
    } else if(strncmp(img_path,"SJI_1400",8) == 0) {
        // berad = be_fsj / RADEG;
        // cdelt = cdlt_fsj;
        berad = be_140 / RADEG;
        cdelt = cdlt_140;
        crpix1 = cpx1_140;
        crpix2 = cpx2_140;
    } else if(strncmp(img_path,"SJI_2832",8) == 0) {
        // berad = be_nsj / RADEG;
        // cdelt = cdlt_nsj;
        berad = be_283 / RADEG;
        cdelt = cdlt_283;
        crpix1 = cpx1_283;
        crpix2 = cpx2_283;
    } else if(strncmp(img_path,"SJI_1600W",9) == 0) {
        // berad = be_fsj / RADEG;
        // cdelt = cdlt_fsj;
        berad = be_mir / RADEG;
        cdelt = cdlt_mir;
        crpix1 = cpx1_mir;
        crpix2 = cpx2_mir;
    } else
	return BAD_IMG_PATH;
    wcsaxes = 2;
    crpix1 = crpix1/sum1 + off1;
    crpix2 = crpix2/sum2 + off2;
    crval1  = crvalxy[0];
    crval2  = crvalxy[1];
    //crpix3 = 0.0;
   // crpix3a = 0.0;
    crval3  = 0.0;
    cdelt1  = cdelt*sum1;
    cdelt2  = cdelt*sum2;
    cdelt3  = 0.0;
    ctype1  = "HPLN-TAN";
    ctype2  = "HPLT-TAN";
    ctype3  = "none";
    cunit1  = "arcsec";
    cunit2  = "arcsec";
    cunit3  = "none";
    pc1_1   =  cos(berad+garad);
    pc1_2   = -sin(berad+garad) * sum2/sum1;
    pc2_1   =  sin(berad+garad) * sum1/sum2;
    pc2_2   =  cos(berad+garad);
    pc3_1   = 0.0;
    pc3_2   = 0.0;
    crpix1a = crpix1;
    crpix2a = crpix2;
    crval1a = crval1;
    crval2a = crval2;
    crval3a = crval3;
    cdelt1a = cdelt1;
    cdelt2a = cdelt2;
    cdelt3a = cdelt3;
    ctype1a = ctype1;
    ctype2a = ctype2;
    ctype3a = ctype3;
    cunit1a = cunit1;
    cunit2a = cunit2;
    cunit3a = cunit3;
    pc1_1a  = pc1_1;
    pc1_2a  = pc1_2;
    pc2_1a  = pc2_1;
    pc2_2a  = pc2_2;
    pc3_1a  = pc3_1;
    pc3_2a  = pc3_2;
    xcen = crval1 + cdelt1*(pc1_1*((naxis1+1)/2. - crpix1) + pc1_2*((naxis2+1)/2. - crpix2));
    ycen = crval2 + cdelt2*(pc2_1*((naxis1+1)/2. - crpix1) + pc2_2*((naxis2+1)/2. - crpix2));
    //lonpole = 180.0;
    //specsys1= "HELIOCEN";    // or 'SOURCE  ' if corrected for solar rot
} else
	return BAD_INSTRUMENT_NAME;

///////////////////////////////////////////////////
if (drms_setkey_int(rs1,"WCSDBVER",version)) return SETKEY_ERROR_1;
if (drms_setkey_int(rs1,"WCSAXES",wcsaxes)) return SETKEY_ERROR_1;
if (drms_setkey_float(rs1,"CRPIX1",crpix1)) return SETKEY_ERROR_2;
if (drms_setkey_float(rs1,"CRPIX2",crpix2)) return SETKEY_ERROR_2;
//if (drms_setkey_float(rs1,"CRPIX3",crpix3)) return SETKEY_ERROR_2;
if (drms_setkey_float(rs1,"CRPIX1A",crpix1a)) return SETKEY_ERROR_2;
if (drms_setkey_float(rs1,"CRPIX2A",crpix2a)) return SETKEY_ERROR_2;
//if (drms_setkey_float(rs1,"CRPIX3A",crpix3a)) return SETKEY_ERROR_2;
if (drms_setkey_float(rs1,"CDELT1",cdelt1)) return SETKEY_ERROR_3;
if (drms_setkey_float(rs1,"CDELT2",cdelt2)) return SETKEY_ERROR_3;
if (drms_setkey_float(rs1,"CDELT3",cdelt3)) return SETKEY_ERROR_3;
if (drms_setkey_float(rs1,"CDELT1A",cdelt1a)) return SETKEY_ERROR_3;
if (drms_setkey_float(rs1,"CDELT2A",cdelt2a)) return SETKEY_ERROR_3;
if (drms_setkey_float(rs1,"CDELT3A",cdelt3a)) return SETKEY_ERROR_3;
if (drms_setkey_string(rs1,"CUNIT1",cunit1)) return SETKEY_ERROR_4;
if (drms_setkey_string(rs1,"CUNIT2",cunit2)) return SETKEY_ERROR_4;
if (drms_setkey_string(rs1,"CUNIT3",cunit3)) return SETKEY_ERROR_4;
if (drms_setkey_string(rs1,"CUNIT1A",cunit1a)) return SETKEY_ERROR_4;
if (drms_setkey_string(rs1,"CUNIT2A",cunit2a)) return SETKEY_ERROR_4;
if (drms_setkey_string(rs1,"CUNIT3A",cunit3a)) return SETKEY_ERROR_4;
if (drms_setkey_string(rs1,"CTYPE1",ctype1)) return SETKEY_ERROR_5;
if (drms_setkey_string(rs1,"CTYPE2",ctype2)) return SETKEY_ERROR_5;
if (drms_setkey_string(rs1,"CTYPE3",ctype3)) return SETKEY_ERROR_5;
if (drms_setkey_string(rs1,"CTYPE1A",ctype1)) return SETKEY_ERROR_5;
if (drms_setkey_string(rs1,"CTYPE2A",ctype2a)) return SETKEY_ERROR_5;
if (drms_setkey_string(rs1,"CTYPE3A",ctype3a)) return SETKEY_ERROR_5;
if (drms_setkey_float(rs1,"PC1_1",pc1_1)) return SETKEY_ERROR_6;
if (drms_setkey_float(rs1,"PC2_1",pc2_1)) return SETKEY_ERROR_6;
if (drms_setkey_float(rs1,"PC3_1",pc3_1)) return SETKEY_ERROR_6;
if (drms_setkey_float(rs1,"PC1_2",pc1_2)) return SETKEY_ERROR_6;
if (drms_setkey_float(rs1,"PC2_2",pc2_2)) return SETKEY_ERROR_6;
if (drms_setkey_float(rs1,"PC3_2",pc3_2)) return SETKEY_ERROR_6;
if (drms_setkey_float(rs1,"PC1_1A",pc1_1a)) return SETKEY_ERROR_6;
if (drms_setkey_float(rs1,"PC2_1A",pc2_1a)) return SETKEY_ERROR_6;
if (drms_setkey_float(rs1,"PC3_1A",pc3_1a)) return SETKEY_ERROR_6;
if (drms_setkey_float(rs1,"PC1_2A",pc1_2a)) return SETKEY_ERROR_6;
if (drms_setkey_float(rs1,"PC2_2A",pc2_2a)) return SETKEY_ERROR_6;
if (drms_setkey_float(rs1,"PC3_2A",pc3_2a)) return SETKEY_ERROR_6;
if (drms_setkey_float(rs1,"CRVAL1",crval1)) return SETKEY_ERROR_7;
if (drms_setkey_float(rs1,"CRVAL2",crval2)) return SETKEY_ERROR_7;
if (drms_setkey_float(rs1,"CRVAL3",crval3)) return SETKEY_ERROR_7;
if (drms_setkey_float(rs1,"CRVAL1A",crval1a)) return SETKEY_ERROR_7;
if (drms_setkey_float(rs1,"CRVAL2A",crval2a)) return SETKEY_ERROR_7;
if (drms_setkey_float(rs1,"CRVAL3A",crval3a)) return SETKEY_ERROR_7; 
if (drms_setkey_float(rs1,"SAT_ROT",scroll)) return SETKEY_ERROR_8; 
if (drms_setkey_float(rs1,"XCEN",xcen)) return SETKEY_ERROR_9;
if (drms_setkey_float(rs1,"YCEN",ycen)) return SETKEY_ERROR_9;

	return 0;
}
