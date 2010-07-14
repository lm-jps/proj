/*
 * Name:		hccsss_pkg.h
 *
 * Description:		Header file that defines the prototypes for hccsss_pkg.c
 *
 * Original source:	IDL code developed by Xuepu Zhao (xuepu@sun.stanford.edu)
 *
 * Written by:		Xudong Sun (xudongs@stanford.edu)
 *
 * Version:		v1.0		Apr 01 2010
 *
 * Issues:
 *				v1.0
 *				Fixed all prototypes
 *				Moved from pfss.h
 *
 */
 
#ifndef _HCCSSS_PKG_H
#define _HCCSSS_PKG_H

#include "pfss.h"
#include "pfss_pkg.h"
#include "fieldline_pkg.h"



/* ############## Matrix inversion from LAPACK library ############## */

void dgetrf_(int *, int *, double *, int *, int *, int *);
void dgetri_(int *, double *, int *, int *, double *, int *, int *);



/* ################ Radial functions for pfss_pkg.h ################ */

/* Radial part of HCCSSS below cusp surface */
void hchf(float r, double *rrr, double *D_rrr, int lmax, float apar);

/* Radial part of HCCSSS above cusp surface */
void cshf(float r, double *rrr, double *D_rrr, int lc, float apar);



/* ################ Harmonic coefficients above RCP ################ */

/* Computing harmonic coefficients above cusp surface */
void gh_hccsss(float *g, float *h, struct Grid *grid, int lmax, 
        float *gc, float *hc, int lc, float apar);



/* ################ HCCSSS Field Vector ################ */

/* Assemble gh into B vectors on a single point, below cusp surface */
void hccsss0_hc(float *g, float *h, struct Point *pt, 
        float *Bvec, int lmax, float apar);

/* Assemble gh into B vectors on a single point, above cusp surface */
void hccsss0_cs(float *g, float *h, struct Point *pt, 
        float *Bvec, int lmax, float apar);

/* Assemble gh into B vectors for HCCSSS on a 3d grid */
void hccsss(float *g, float *h, float *gc, float *hc, struct Grid *grid, 
        float *Br, float *Bt, float *Bp, int lmax, int lc, float apar);


#endif
