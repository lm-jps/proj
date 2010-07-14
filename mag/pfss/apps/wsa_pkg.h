/*
 * Name:		wsa_pkg.h
 *
 * Description:		Header file that defines the prototypes for wsa_pkg.c
 *
 * Original source:	IDL code developed by Xuepu Zhao (xuepu@sun.stanford.edu)
 *
 * Written by:		Xudong Sun (xudongs@stanford.edu)
 *
 * Version:		v1.0		Apr 01 2010
 *				v1.1		Apr 19 2010
 *
 * Issues:
 *				v1.0
 *				Fixed all prototypes
 *				Moved from pfss.h
 *				v1.1
 *				Changed imf to int, date to double (sec)
 *
 */
 
#ifndef _WSA_PKG_H
#define _WSA_PKG_H

#include "pfss.h"
#include "pfss_pkg.h"
#include "fieldline_pkg.h"



/* ################ Empirical Function ################ */

/* WSA empirical function for speed */
float wsa_f(float f, float th);



/* ################ Angular separation ################ */

/* Great circle distance between pt1 and pt2 */
float gcd(float lon1, float lat1, float lon2, float lat2);

/* Smooth foot point map to avoid small closed structure */
void smoothfp(int *fpmap, int np, int nt);

/* Get open/close at each grid point based on downward */
void getbd(int *fpmap, float *fte, float *fpph, float *fpth, 
	int np, int nt);

/* Find minimum angular separation */
float ang_sep(float ph, float th, int *fpmap, int np, int nt);



/* ################ Getting the physical variables ################ */

/* Global field line tracing for flux tube expansion factor and foot point locations */
void fte_global(float *g, float *h, struct Grid *grid, float *fte, 
	float *fpph, float *fpth, float *brss, int lmax, 
	void (*integrator)(struct Point *, struct Point *, float *, float *, int, float, float, 
		void (*)(float *, float *, struct Point *, float *, int, float)));

/* Subearth field line tracing for flux tube expansion factor and foot point locations */
void fte_subearth(float *g, float *h, struct Gridsub *gridsub, float *fte_sub, 
	float *fpph_sub, float *fpth_sub, float *brss_sub, int lmax, 
	void (*integrator)(struct Point *, struct Point *, float *, float *, int, float, float, 
		void (*)(float *, float *, struct Point *, float *, int, float)));



/* ################ Propagating to 1AU ################ */

/* v, b: Series at source surface; v1AU, b1AU: Series at 1AU */
void propagate(double *t, float *v, float *b, double *t1AU, float *v1AU, 
	float *b1AU, int np);

/* Linear interpolation function */
void interpol(float *y0, double *x0, int n0, float *y, double *x, int n);

/* Takes the predicted series at 1AU average and interpolate onto a uniform grid */
void fixts1au(double *t, float *v, float *b, int n, 
	double *date, float *speed, int *imf, int num);



/* ################ Main WSA ################ */

/* Global and subearth WSA code at source surface */
void wsa_ss(float *g, float *h, int lmax, int rk4, 
	struct Grid *grid, struct Gridsub *gridsub, 
	float *v_g, float *brss_g, float *fte_g, float *fpph_g, float *fpth_g, 
	float *v_s, float *brss_s, float *fte_s, float *fpph_s, float *fpth_s);

/* Subearth WSA code at 1AU */
void wsa_1AU(struct Gridsub *gridsub, float *v_s, float *brss_s, double *date, 
	float *speed, int *imf, int tslen);
	
	
#endif
