/*
 * Module name:		polar_remap.c
 *
 * Description:
 *          
 *
 * Written by:		Xudong Sun (xudongs@stanford.edu)
 *
 */

#include <math.h>
#include "interp.c"

/* Stereographic projection converting a synoptic map to polar view
 * Reverse projection formula is used to find the map coordinates
 * of a certain image point. A bilinear interpolation is then used to get the
 * value. in the image. During projection, the pole is projected to the center
 * of the image, and the lower latitude limit is projected to an elipse that is
 * tangent to the edge of the image. Note there's no error checking here. Any
 * neighboring missing data during interpolation results in NaN.
 */

void synop2pole(int ns, int sinlat, double *image, int *imgDims,
                double *map, int *mapDims,
                double *lons, double *lats, double latLim)
{
    
    int nx = imgDims[0], ny = imgDims[1];       // image
    int xsz = mapDims[0], ysz = mapDims[1];     // map
    
    // Image center
    
    double xc = (nx - 1.) / 2.0;
    double yc = (ny - 1.) / 2.0;
    
    // Normalized radius
    
    double rLim = cos(latLim) / (1 + sin(fabs(latLim)));    // 1 for equator
    double rx = xc / rLim, ry = yc / rLim;
    
    // Do
    
    double x0, y0;      // normalized
    double r, th;       // Polar
    double lon, lat;
    double ind_lon, ind_lat;   // Indices for interpolation
    
    for (int row = 0; row < ny; row++) {
        for (int col = 0; col < nx; col++) {
            
            // Backward projection, find lon/lat for x0, y0
            
            x0 = (col - xc) / rx;
            y0 = (row - yc) / ry;
            r = sqrt(x0 * x0 + y0 * y0);
            th = atan2(y0, x0);
            lon = th;
            lat = M_PI / 2.0 - atan(r) * 2.0;
            if (ns) lat = -lat;
            
            // Within lat_lim? Removed since NAN's are being fitted too
            
            /*
            if (fabs(lat) < fabs(latLim)) {
                image[row * nx + col] = DRMS_MISSING_DOUBLE;
                continue;
            }*/
            
            // Get indices for interp
            
            ind_lon = get_index(0, xsz, lons, lon);
            ind_lat = get_index(sinlat, ysz, lats, lat);      // sinlat

            // Do interp, bilinear
            
            image[row * nx + col] = linintd(map, xsz, ysz, ind_lon, ind_lat);
            
        }
    }
    
    //
    
    return;
}

void pole2synop(int ns, int sinlat, double *image, int *imgDims,
                double *map, int *mapDims,
                double *lons, double *lats, double latLim)
{
    
    int nx = imgDims[0], ny = imgDims[1];       // image
    int xsz = mapDims[0], ysz = mapDims[1];     // map
    
    // Image center
    
    double xc = (imgDims[0] - 1.) / 2.0;
    double yc = (imgDims[1] - 1.) / 2.0;
    
    // Normalized radius
    
    double rLim = cos(latLim) / (1 + sin(fabs(latLim)));    // 1 for equator
    double rx = xc / rLim, ry = yc / rLim;

    // Find range
    
//    printf("  xsz=%d, ysz=%d\n", xsz, ysz);
    
    int row0, row1;     // lower, upper limit to fill
    if (ns) {
        row0 = 0;
        row1 = ceil(get_index(sinlat, ysz, lats, -fabs(latLim)));
    } else {
        row0 = floor(get_index(sinlat, ysz, lats, fabs(latLim)));
        row1 = ysz - 1;
    }
//    printf("  row0=%d, row1=%d\n", row0, row1);
    
    // Interp
    
    double x0, y0;	// Nomalized
    double x, y;
    double r, th;	// Polar coord
    double lon, lat;
    
    for (int row = row0; row <= row1; row++) {
        for (int col = 0; col < xsz; col++) {
            
            lon = lons[col];
            lat = lats[row];
            
            // Forward projection, find x, y for lon/lat
            
            r = cos(lat) / (1 + sin(fabs(lat)));
            th = lon;
            x0 = r * cos(th);
            y0 = r * sin(th);
            x = x0 * rx + xc;
            y = y0 * ry + yc;
            
            // Do interp
            
            map[row * xsz + col] = linintd(image, nx, ny, x, y);
            
        }
    }
    
    //
    
    return;
}