/*
 * Name:		fieldline_pkg.c
 *
 * Description:		This is the C package for field line tracing.
 *			Including method of open field line tracing (both globally and
 *			on subearth line), using Euler of RK4 integration.
 *			Also provided are 
 *                      Original code of the tracing is in IDL by Xuepu Zhao.
 *			Time and solar coordinate????
 *			Restructured and optimized by Xudong Sun.
 *
 * Function List:
 *			void copypoint(struct Point *pt1, struct Point *pt2)
 *			void sph2cart(struct Point *pt)
 *			void cart2sph(struct Point *pt)
 *			void vec_sph2car(float *vec_s, float *vec_c, struct Point *pt)
 *			void differential(float *dydx, float *Y, float *g, float *h, 
 *				int lmax, float apar, float dir, void (*Bfunc)())
 *			void RK4(struct Point *pt1, struct Point *pt2, float *g, float *h, 
 *				int lmax, float apar, float dir, void (*Bfunc)())
 *			void Euler(struct Point *pt1, struct Point *pt2, float *g, float *h, 
 *				int lmax, float apar, float dir, void (*Bfunc)())
 *			void compute_fte(float *g, float *h, struct FldLn *fl, int lmax, 
 *				float r_hi, float r_lo, float apar, void (*Bfunc)())
 *			void trace(float *g, float *h, struct Point *pt0, struct FldLn *fl, int lmax, 
 *				void (*integrator)(), void (*Bfunc)(), float r_hi, float r_lo, float apar)
 *
 *
 * Calling:		pfss_pkg.c
 * Called by:		jgh2wsa.c
 *
 * Original source:	IDL code developed by Xuepu Zhao (xuepu@sun.stanford.edu)
 *			IDL code developed for NLFFF by T. Wiegelmanm (wiegelmann@linmpi.mpg.de)
 * Adapted by:		Xudong Sun (xudongs@stanford.edu)
 *
 * Version:
 *			v1.0		Jul 29 2009
 *			v1.1		Aug 03 2009
 *			v1.1a		Aug 05 2009
 *			v1.1b		Oct 28 2009
 *			v1.2		Oct 29 2009
 *
 * Issues:
 *			v1.0
 *			Need to include "pfss_pkg.c" for data structure and constant definition?
 *			Could add a function that can start tracing from anypoint (works for loops then)
 *			For downward tracing the ending point are not on the same sphere (interpolation?)
 *			v1.1
 *			Add outside control of integrator in trace
 *			v1.1a
 *			Moved constants and data structures to pfss.h, called from there
 *			v1.1b
 *			Added r_hi, r_lo to trace() and compute_fte() for hccsss' use
 *			v1.2
 *			Added apar for most functions
 *			Substitude pfss0 for a more general Bfunc that could be point to different
 *			model calls. See fte_global() in wsa_pkg.c for calling example
 *			Corrected a call for compute_fte() in trace(), was missing apar (nov 05 2009)
 *
 */

#include "fieldline_pkg.h"


/* ################## Utils ################## */


/* Copy point1 to point2 */
void copypoint(struct Point *pt1, struct Point *pt2)
{
    pt2->r = pt1->r; pt2->t = pt1->t; pt2->p = pt1->p;
    pt2->x = pt1->x; pt2->y = pt1->y; pt2->z = pt1->z;
}


/* Spherical to Cartesian */
void sph2cart(struct Point *pt)
{
    pt->x = pt->r * sin(pt->t) * cos(pt->p);
    pt->y = pt->r * sin(pt->t) * sin(pt->p);
    pt->z = pt->r * cos(pt->t);
}


/* Cartesian to Spherical */
void cart2sph(struct Point *pt)
{
    pt->r = sqrt(pt->x * pt->x + pt->y * pt->y + pt->z * pt->z);
    pt->t = acos(pt->z / pt->r);
    pt->p = atan2(pt->y, pt->x);
    pt->p = (pt->p > 0) ? pt->p : (pt->p + 2.0 * M_PI);
}


/* Cartesian to Spherical, normalized, for 3D vector */
void vec_sph2car(float *vec_s, float *vec_c, struct Point *pt)
{
    double costh = cos(pt->t), sinth = sin(pt->t);
    double cosph = cos(pt->p), sinph = sin(pt->p);
    double vec_abs = sqrt(vec_s[0] * vec_s[0] + vec_s[1] * vec_s[1] + vec_s[2] * vec_s[2]);
    vec_c[0] = (vec_s[0] * sinth * cosph + vec_s[1] * costh * cosph - vec_s[2] * sinph) / vec_abs;
    vec_c[1] = (vec_s[0] * sinth * sinph + vec_s[1] * costh * sinph + vec_s[2] * cosph) / vec_abs;
    vec_c[2] = (vec_s[0] * costh - vec_s[1] * sinth) / vec_abs;
}



/* ################## Field Line Tracing Tools ################## */


/* This part could be rewritten for field line tracing in
 * other models: simply point Bfunc to another field
 * vector function. Vectors are converted into Cartesian coordinates.
 */
void differential(float *dydx, float *Y, float *g, float *h, 
	int lmax, float apar, float dir, 
	void (*Bfunc)(float *, float *, struct Point *, float *, int, float))
{
    struct Point *pt = (struct Point *)(malloc(sizeof(struct Point)));;
    pt->x = Y[0]; pt->y = Y[1]; pt->z = Y[2];
    cart2sph(pt);
    float Bvec[3];
    (*Bfunc)(g, h, pt, Bvec, lmax, apar);		// Field vector, spherical coordinate
    vec_sph2car(Bvec, dydx, pt);		// Convert B vector back to Cartesian
    for (int i = 0; i < 3; i++) 
        dydx[i] *= (dir * pt->r);		// Direction, proportional to radius
    free(pt);
}


// Ad hoc RK4 for field line tracing from pt1 to pt2
// Use "differential" for derivatives, step is constant STEP
// Dir controls the tracing direction
void RK4(struct Point *pt1, struct Point *pt2, float *g, float *h, 
	int lmax, float apar, float dir, 
	void (*Bfunc)(float *, float *, struct Point *, float *, int, float))
{
    float Y[3], Yt[3];
    float k1[3], k2[3], k3[3], k4[3];

    // Use Cartesian here
    Y[0] = pt1->x; Y[1] = pt1->y; Y[2] = pt1->z;
    differential(k1, Y, g, h, lmax, apar, dir, Bfunc);
    for (int i = 0; i < 3; i++) Yt[i] = Y[i] + STEP * k1[i] / 2.0;
    differential(k2, Yt, g, h, lmax, apar, dir, Bfunc);
    for (int i = 0; i < 3; i++) Yt[i] = Y[i] + STEP * k2[i] / 2.0;
    differential(k3, Yt, g, h, lmax, apar, dir, Bfunc);
    for (int i = 0; i < 3; i++) Yt[i] = Y[i] + STEP * k3[i];
    differential(k4, Yt, g, h, lmax, apar, dir, Bfunc);
    pt2->x = Y[0] + STEP * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6.0;
    pt2->y = Y[1] + STEP * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6.0;
    pt2->z = Y[2] + STEP * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]) / 6.0;
}

// Ad hoc Euler integrator, from pt1 to pt2
void Euler(struct Point *pt1, struct Point *pt2, float *g, float *h, 
	int lmax, float apar, float dir, 
	void (*Bfunc)(float *, float *, struct Point *, float *, int, float))
{
    float Y[3], k[3];
    Y[0] = pt1->x; Y[1] = pt1->y; Y[2] = pt1->z;
    differential(k, Y, g, h, lmax, apar, dir, Bfunc);
    pt2->x = Y[0] + STEP * k[0];
    pt2->y = Y[1] + STEP * k[1];
    pt2->z = Y[2] + STEP * k[2];
}



/* ################## Field Line Tracing Package ################## */


/* Compute flux tube expansion factor */
void compute_fte(float *g, float *h, struct FldLn *fl, int lmax, 
	float r_hi, float r_lo, float apar, 
	void (*Bfunc)(float *, float *, struct Point *, float *, int, float))
{
    struct Point *pt = (struct Point *)(malloc(sizeof(struct Point)));
    pt->r = r_lo;
    pt->t = fl->tt[fl->num];
    pt->p = fl->pp[fl->num];
    float Bvec[3];
    (*Bfunc)(g, h, pt, Bvec, lmax, apar);
    fl->fte = (Bvec[0] / fabs(fl->brss)) * (r_lo * r_lo) / (r_hi * r_hi);
    free(pt);
}



/* Ad hoc field line tracing code, tracing downward from Rss
 * Use either RK4 or Euler for integration
 * Note: this fuction only deals with open field line
 * a more general loop tracing can be done starting from any point
 */
void trace(float *g, float *h, struct Point *pt0, struct FldLn *fl, int lmax, 
	void (*integrator)(struct Point *, struct Point *, float *, float *, int, float, float, 
		void (*)(float *, float *, struct Point *, float *, int, float)), 
	void (*Bfunc)(float *, float *, struct Point *, float *, int, float), 
	float r_hi, float r_lo, float apar)
{
//    printf("%f, %f\n", r_hi, r_lo);
    // Auxiliary points
    struct Point *pt = (struct Point *)(malloc(sizeof(struct Point)));
    struct Point *ptt = (struct Point *)(malloc(sizeof(struct Point)));
//printf("b, ");
    // Preparation
    fl->op_cl = 1;	// Open field line
    float bvec[3];
//printf("a, ");
    (*Bfunc)(g, h, pt0, bvec, lmax, apar);	// Field at source surface
    fl->brss = bvec[0];
    float dir = (bvec[0] > 0) ? -1. : 1.;	// Reverse dir for downward tracing

    // Tracing
    sph2cart(pt0);
    copypoint(pt0, pt);	// Copy pt0 to pt
    int count = 0;
//    printf("%f, %f, %d\n", r_hi, r_lo, count);
    while (pt->r <= r_hi && pt->r >= r_lo && count < MAXSTP) {
//        printf("here ,");
        fl->rr[count] = pt->r;
        fl->tt[count] = pt->t;
        fl->pp[count] = pt->p;
        // From pt to ptt, ad hoc, Euler or RK4 integrator
        (*integrator)(pt, ptt, g, h, lmax, apar, dir, Bfunc);
        cart2sph(ptt);
        copypoint(ptt, pt);	// Copy ptt to pt
        count++;
    }
    if (count < MAXSTP) {
        fl->num = count - 1;	// If count equals 1 then num says 0 (error)
        compute_fte(g, h, fl, lmax, r_hi, r_lo, apar, Bfunc);
    }
    free(pt); free(ptt);
}

