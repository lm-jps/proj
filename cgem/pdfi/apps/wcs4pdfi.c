/*
 *
 * Here are two functions that convert CGEM PDFI related geometry parameters to WCS coordinates.
 * For PDFI software, see http://cgem.ssl.berkeley.edu/cgi-bin/cgem/PDFI_SS/home
 *
 * Written by X. Sun (xudong@sun.stanford.edu)
 *
 * PDFI paramters:
 *  n: column # (input is n+1); m: row # (input is m+1)
 *  grid: enum staggeredGrids
 *  a: MINCOLAT; b: MAXCOLAT; c: MINLON; d: MAXLON (in rad)
 *
 * WCS parameters:
 *  crval1/2, crpix1/2, cdelt1/2 (in degree)
 *
 * Notes:
 *  This assumes uniform lat/lon grid. Reference lon is arbitrary
 *
 * History
 *  2016 Dec 02: v0
 *  2018 Mar 26: v1     a/b/c/d are limits of EDGE, fixed cdelt1/2
 *
 */

#include <math.h>

#define DTOR		(M_PI/180.)

// Grid names
enum staggeredGrids {
    COE, CE, CO, TE, PE
};

// Grid sizes, return n+1, m+1 if unspecified
void sizeofGrid(enum staggeredGrids grid, int n, int m, int *cols, int *rows)
{
    switch (grid) {
        case COE:       // B_obs, V_obs, E_r
            *cols = n + 1; *rows = m + 1;
            break;
        case CE:        // B_r, S_r
            *cols = n; *rows = m;
            break;
        case CO:        // div E_h, curl B_h
            *cols = n - 1; *rows = m - 1;
            break;
        case TE:        // B_t, V_t, E_p
            *cols = n; *rows = m + 1;
            break;
        case PE:        // B_p, V_p, E_t
            *cols = n + 1; *rows = m;
            break;
        default:        // Default: same as obs
            *cols = n + 1; *rows = m + 1;
            break;
    }
}

/*
 * From PDFI to WCS
 */

void pdfi2wcs(int n, int m,
              enum staggeredGrids grid,
              double a, double b, double c, double d,
              double *crval1, double *crval2,
              double *crpix1, double *crpix2,
              double *cdelt1, double *cdelt2)
{
    int cols, rows;
    sizeofGrid(grid, n, m, &cols, &rows);
    
    double minlat = 90. - b / DTOR;
    double maxlat = 90. - a / DTOR;
    double minlon = c / DTOR;
    double maxlon = d / DTOR;
    
    *crval1 = (minlon + maxlon) / 2.;
    *crval2 = (minlat + maxlat) / 2.;
    *crpix1 = (1. + cols) / 2.;
    *crpix2 = (1. + rows) / 2.;
    *cdelt1 = (maxlon - minlon) / cols;             // Fixed from cols-1
    *cdelt2 = (maxlat - minlat) / rows;
}

/*
 * From WCS to PDFI
 */

void wcs2pdfi(int n, int m,
              enum staggeredGrids grid,
              double crval1, double crval2,
              double crpix1, double crpix2,
              double cdelt1, double cdelt2,
              double *a, double *b, double *c, double *d)
{
    int cols, rows;
    sizeofGrid(grid, n, m, &cols, &rows);
    
    double minlat = crval2 + (0.5 - crpix2) * cdelt2;           // Edge is at 0.5 and rows+0.5
    double maxlat = crval2 + (rows + 0.5 - crpix2) * cdelt2;
    double minlon = crval1 + (0.5. - crpix1) * cdelt1;
    double maxlon = crval1 + (cols + 0.5 - crpix1) * cdelt1;
    
    *a = (90. - maxlat) * DTOR;
    *b = (90. - minlat) * DTOR;
    *c = minlon * DTOR;
    *d = maxlon * DTOR;
}
