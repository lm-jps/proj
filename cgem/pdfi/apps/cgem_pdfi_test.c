/*
 *  cgem_prep.c
 *
 *  Example:
 *  cgem_pdfi "in=hmi_test.cgem_pdfi_in[11158][2011.02.15_01:24/30m]" "out=hmi_test.cgem_pdfi_out"
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "jsoc_main.h"
#include "astro.h"

// Legacy macros

#define PI              (M_PI)
#define RADSINDEG		(PI/180.)
#define RAD2ARCSEC		(648000./M_PI)
#define SECINDAY		(86400.)
#define FOURK			(4096)
#define FOURK2          (16777216)

#define DIE(msg) {fflush(stdout); fprintf(stderr,"%s, status=%d\n", msg, status); return(status);}
#define SHOW(msg) {printf("%s", msg); fflush(stdout);}

#define kNotSpecified "Not Specified"

// FORTRAN function

extern void pdfi_wrapper4jsoc_ss_(int *m, int *n, double *a, double *b, double *c, double *d,
                                  double *bloncoe0, double *blatcoe0, double *brllcoe0,
                                  double *bloncoe1, double *blatcoe1, double *brllcoe1,
                                  double *vloncoe0, double *vlatcoe0, double *vrllcoe0,
                                  double *vloncoe1, double *vlatcoe1, double *vrllcoe1,
                                  double *lloncoe0, double *llatcoe0, double *lrllcoe0,
                                  double *lloncoe1, double *llatcoe1, double *lrllcoe1,
                                  double *t0, double *t1,
                                  double *blon0, double *blat0, double *brll0,
                                  double *blon1, double *blat1, double *brll1,
                                  double *elonpdfi, double *elatpdfi, double *erllpdfi);

// =====================================

char *module_name = "cgem_pdfi_test";

ModuleArgs_t module_args[] =
{
//    {ARG_STRING, "in", kNotSpecified, "Input data series."},
//    {ARG_STRING, "out", kNotSpecified, "Output data series."},
    {ARG_END}
};

int DoIt(void)
{
    
    int status = DRMS_SUCCESS;
    
    /* Input */
    
    int n, m;           // n column, m row
    double a, b, c, d;
    double t0, t1;
    
    double *bloncoe0, *blatcoe0, *brllcoe0;     // B vectors
    double *bloncoe1, *blatcoe1, *brllcoe1;
    
    double *vloncoe0, *vlatcoe0, *vrllcoe0;     // V vectors
    double *vloncoe1, *vlatcoe1, *vrllcoe1;
    
    double *lloncoe0, *llatcoe0, *lrllcoe0;     // los vectors
    double *lloncoe1, *llatcoe1, *lrllcoe1;
    
    // temp values
    
    n = 554;
    m = 532;
    a = 1.444096088409;
    b = 1.720951676369;
    c = 0.000000000000;
    d = 0.288267821074;
    t0 = 2455607.5582638895139098 * 3600. * 24.;
    t1 = 2455607.5749305561184883 * 3600. * 24.;

    
    bloncoe0 = (double *) (malloc((n + 1) * (m + 1) * sizeof(double)));
    blatcoe0 = (double *) (malloc((n + 1) * (m + 1) * sizeof(double)));
    brllcoe0 = (double *) (malloc((n + 1) * (m + 1) * sizeof(double)));
    bloncoe1 = (double *) (malloc((n + 1) * (m + 1) * sizeof(double)));
    blatcoe1 = (double *) (malloc((n + 1) * (m + 1) * sizeof(double)));
    brllcoe1 = (double *) (malloc((n + 1) * (m + 1) * sizeof(double)));
    
    vloncoe0 = (double *) (malloc((n + 1) * (m + 1) * sizeof(double)));
    vlatcoe0 = (double *) (malloc((n + 1) * (m + 1) * sizeof(double)));
    vrllcoe0 = (double *) (malloc((n + 1) * (m + 1) * sizeof(double)));
    vloncoe1 = (double *) (malloc((n + 1) * (m + 1) * sizeof(double)));
    vlatcoe1 = (double *) (malloc((n + 1) * (m + 1) * sizeof(double)));
    vrllcoe1 = (double *) (malloc((n + 1) * (m + 1) * sizeof(double)));
    
    lloncoe0 = (double *) (malloc((n + 1) * (m + 1) * sizeof(double)));
    llatcoe0 = (double *) (malloc((n + 1) * (m + 1) * sizeof(double)));
    lrllcoe0 = (double *) (malloc((n + 1) * (m + 1) * sizeof(double)));
    lloncoe1 = (double *) (malloc((n + 1) * (m + 1) * sizeof(double)));
    llatcoe1 = (double *) (malloc((n + 1) * (m + 1) * sizeof(double)));
    lrllcoe1 = (double *) (malloc((n + 1) * (m + 1) * sizeof(double)));
    
    /* Output */
    
    double *blon0, *blat0, *brll0;      // staggered B vector
    double *blon1, *blat1, *brll1;
    
    double *elonpdfi, *elatpdfi, *erllpdfi;     // output E vector
    
    // temp values
    
    blon0 = (double *) (malloc((n + 1) * m * sizeof(double)));
    blat0 = (double *) (malloc(n * (m + 1) * sizeof(double)));
    brll0 = (double *) (malloc(n * m * sizeof(double)));
    blon1 = (double *) (malloc((n + 1) * m * sizeof(double)));
    blat1 = (double *) (malloc(n * (m + 1) * sizeof(double)));
    brll1 = (double *) (malloc(n * m * sizeof(double)));
    
    elonpdfi = (double *) (malloc(n * (m + 1) * sizeof(double)));
    elatpdfi = (double *) (malloc((n + 1) * m * sizeof(double)));
    erllpdfi = (double *) (malloc((n + 1) * (m + 1) * sizeof(double)));
    
    /* Try calling FORTRAN wrapper */
    
    pdfi_wrapper4jsoc_ss_(&m, &n, &a, &b, &c, &d,
                          bloncoe0, blatcoe0, brllcoe0, bloncoe1, blatcoe1, brllcoe1,
                          vloncoe0, vlatcoe0, vrllcoe0, vloncoe1, vlatcoe1, vrllcoe1,
                          lloncoe0, llatcoe0, lrllcoe0, lloncoe1, llatcoe1, lrllcoe1,
                          &t0, &t1,
                          blon0, blat0, brll0, blon1, blat1, brll1,
                          elonpdfi, elatpdfi, erllpdfi);
    
    /* Some tests */
    
    
    
    /* Clean up */
    
    free(bloncoe0); free(blatcoe0); free(brllcoe0);
    free(bloncoe1); free(blatcoe1); free(brllcoe1);
    free(vloncoe0); free(vlatcoe0); free(vrllcoe0);
    free(vloncoe1); free(vlatcoe1); free(vrllcoe1);
    free(lloncoe0); free(llatcoe0); free(lrllcoe0);
    free(lloncoe1); free(llatcoe1); free(lrllcoe1);
    
    free(blon0); free(blat0); free(brll0);
    free(blon1); free(blat1); free(brll1);
    free(elonpdfi); free(elatpdfi); free(erllpdfi);
    
    //
    
    return DRMS_SUCCESS;
}