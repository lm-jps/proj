/*
 * Name:		pfss.h
 *
 * Description:		Header file that defines the constants and special data 
 *			for structures basic functions in the PFSS package.
 *
 * Original source:	IDL code developed by Xuepu Zhao (xuepu@sun.stanford.edu)
 *
 * Written by:		Xudong Sun (xudongs@stanford.edu)
 *
 * Version:		v1.0		Aug 06 2009
 *				v1.1		Apr 01 2010
 *
 * Issues:
 *				v1.1
 *				Fixed all prototypes and moved them out to different header files
 *
 */

#ifndef _PFSS_H
#define _PFSS_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#ifndef DTOR
#define DTOR (M_PI / 180.)
#endif

/* #################### Constants #################### */


/* pfss_pkg */

#define R0	1.0		// Photosphere
#define RS	2.5		// Source surface
#define RCP	2.5		// Cusp surface for HCCSSS
#define RSS	15.0		// Source surface for HCCSSS
#define R0RS	(R0 / RS)
#define R0RS2	(R0RS * R0RS)

/* fieldline_pkg */

#define STEP	0.005		// Step size in tracing
#define MAXSTP	600		// Maximum step of field line tracing

/* wsa_pkg */

#define BOX	40.0		// Search boundary (deg) for angular separation
#define RSUN		6.955E5		// Solar radius in km
#define PROPDIST	(215.0 - RS)	// Propagation distance from source surface to 1AU
#define DAYSINCR	27.27526	// Days in a Carrington Rotation
#define SECSINDAY	8.64E4		// Seconds in a day
#define NSTEP		10	// Steps of propagation



/* ################## Data Structures ################## */


/* pfss_pkg */

/* Grid to compute on */
struct Grid {
    int np;
    int nt;
    int nr;
    float *ph;
    float *th;
    float *rr;
};
/* A temperary storage of combined coefficient */
struct CombinedCoeff {
    double *rrr;	// (lmax + 1)
    double *D_rrr;
    double *leg;	// (lmax + 1)(lmax + 2)/2
    double *D_leg;
    double *mgh;	// (lmax + 1)(lmax + 2)/2
    double *D_mgh;
};
/* Coordinate of a point, both Cartesian and spherical */
struct Point {
    float r;		// Spherical
    float t;
    float p;
    float x;		// Cartesian
    float y;
    float z;
};


/* fieldline_pkg */

/* Field line */
struct FldLn {
    int num;		// Number of points, 0 for failure
    int op_cl;		// Open field line (1), closed (0)
    float fte;		// Flux tube expansion factor
    float brss;		// Radial field at source surface
    float rr[MAXSTP];	// Point coordinates on field line
    float tt[MAXSTP];
    float pp[MAXSTP];
};


/* wsa_pkg */

/* Grid for subearth computation, for WSA */
struct Gridsub {
    int np;
    int nt;
    double *time;	// CMP time of each subearth point (since 1977.01.01) [NP]
    float *b0;		// b angle of corresponding ph [NP]
    float *ph;		// Longitude of subearth point [NP]
    float *th;		// Latitude of subearth point, N/S dth [NPx3]
};

#endif
