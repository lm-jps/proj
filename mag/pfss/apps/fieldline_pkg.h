/*
 * Name:		fieldline_pkg.h
 *
 * Description:		Header file that defines the prototypes for fieldline_pkg.c
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
 
#ifndef _FIELDLINE_PKG_H
#define _FIELDLINE_PKG_H

#include "pfss.h"
#include "pfss_pkg.h"



/* ################## Utils ################## */

/* Copy point1 to point2 */
void copypoint(struct Point *pt1, struct Point *pt2);

/* Spherical to Cartesian */
void sph2cart(struct Point *pt);

/* Cartesian to Spherical */
void cart2sph(struct Point *pt);

/* Cartesian to Spherical, normalized, for 3D vector */
void vec_sph2car(float *vec_s, float *vec_c, struct Point *pt);



/* ################## Field Line Tracing Tools ################## */

/* Common function for derivative of field strength at a point */
void differential(float *dydx, float *Y, float *g, float *h, 
	int lmax, float apar, float dir, 
	void (*Bfunc)(float *, float *, struct Point *, float *, int, float));

/* Ad hoc RK4 for field line tracing from pt1 to pt2 */
void RK4(struct Point *pt1, struct Point *pt2, float *g, float *h, 
	int lmax, float apar, float dir, 
	void (*Bfunc)(float *, float *, struct Point *, float *, int, float));

/* Ad hoc Euler for field line tracing from pt1 to pt2 */
void Euler(struct Point *pt1, struct Point *pt2, float *g, float *h, 
	int lmax, float apar, float dir, 
	void (*Bfunc)(float *, float *, struct Point *, float *, int, float));



/* ################## Field Line Tracing Package ################## */

/* Compute flux tube expansion factor */
void compute_fte(float *g, float *h, struct FldLn *fl, int lmax, 
	float r_hi, float r_lo, float apar, 
	void (*Bfunc)(float *, float *, struct Point *, float *, int, float));

/* Ad hoc field line tracing code */
void trace(float *g, float *h, struct Point *pt0, struct FldLn *fl, int lmax, 
	void (*integrator)(struct Point *, struct Point *, float *, float *, int, float, float, 
		void (*)(float *, float *, struct Point *, float *, int, float)), 
	void (*Bfunc)(float *, float *, struct Point *, float *, int, float), 
	float r_hi, float r_lo, float apar);


#endif
