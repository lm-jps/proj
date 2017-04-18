/*
 *  pfss_q.c
 *
 *  This module performs smoothing on full-resolution synoptic maps with
 *  3 deg fwhm Gaussian and resample to 361x181, following Toth+ (2011).
 *  PFSS is computed to lmax=120, on a 361x181x50 grid.
 *  Additional Br is saved for 10 layers at 1441x721 where Q will be computed
 *
 *  Author:
 *      Xudong Sun, Based on Xuepu Zhao's Fortran code
 *
 *  Version
 *      v0.0    Apr 18 2017
 *
 *  Notes:
 *      v0.0
 *          Fixed grid
 *
 *
 *  Example call:
 *      pfss_q "in=su_xudong.synop_polfil[2099]" "out=su_xudong.pfss_q"
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>
#include "jsoc_main.h"
#include "radius.h"

#define PI              (M_PI)
#define TWOPI           (2.*M_PI)
#define RADSINDEG		(PI/180.)
#define RAD2ARCSEC		(648000./M_PI)

#define DIE(msg) {fflush(stdout); fprintf(stderr,"%s, status=%d\n", msg, status); return(status);}
#define SHOW(msg) {printf("%s", msg); fflush(stdout);}

#define kNotSpecified "Not Specified"

#define ARRLENGTH(ARR) (sizeof(ARR) / sizeof(ARR[0]))
#define NR_OUT      (ARRLENGTH(r_out))

/* ====================================================================================== */


char *module_name = "pfss_q";

ModuleArgs_t module_args[] =
{
    {ARG_STRING, "in", kNotSpecified, "Input synoptic map."},
    {ARG_STRING, "out", kNotSpecified, "Output PFSS series."},
    {ARG_END}
};

int DoIt(void)
{
    
    int status = DRMS_SUCCESS;

    //
    
    return DRMS_SUCCESS;
    
}

