/*
 * Name:		pfss_pkg.h
 *
 * Description:		Header file that defines the prototypes for pfss_pkg.c
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
 
#ifndef _PFSS_PKG_H
#define _PFSS_PKG_H

#include "pfss.h"



/* ########## Coefficients in three direction ######### */

/* Radial part for PFSS */
void rdr(float r, double *rrr, double *D_rrr, int lmax, float apar);

/* Theta part, Legendre function */
void pdpth(float th, double *leg, double *D_leg, int lmax);

/* Phi part, combine g and h with phi */
void csmph(float ph, double *mgh, double *D_mgh, float *g, float *h, int lmax);



/* ################# GH Coefficient ################# */

/* Use numerical integration to get gh */
void gh_pfss(float *map, float *g, float *h, struct Grid *grid, 
	int lmax, int sinlat);



/* ################# Field Strength ################# */

/* Compute Bvec at a point, geometries in CmbCoeff */
void Brtp(float r, float sint, float *Br0, float *Bt0, float *Bp0, 
	struct CombinedCoeff *CmbCoeff, int lmax);
	
/* Common field function for cube used by all PFSS-like models */
void Bcube(float *g, float *h, struct Grid *grid, float *Br,
	float *Bt, float *Bp, int lmax, float apar, void (*rfunc)(float, double *, double *, int, float));

/* Common field function for point used by all PFSS-like models */
void Bpoint(float *g, float *h, struct Point *pt, 
	float *Bvec, int lmax, float apar, void (*rfunc)(float, double *, double *, int, float));



/* ################ Calls for PFSS ################ */

/* Assemble gh into B vectors on a single point */	
void pfss0(float *g, float *h, struct Point *pt, 
	float *Bvec, int lmax, float apar);

/* Assemble gh into B vectors on the grid */
void pfss(float *g, float *h, struct Grid *grid, 
	float *Br, float *Bt, float *Bp, int lmax, float apar);

	
#endif
