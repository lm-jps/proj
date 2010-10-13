
/*
 * THIS IS A GENERATED FILE; DO NOT EDIT.
 *
 * Declarations for calling `roi_stats_mag' a.k.a. `rsm' as a C library.
 *
 * Made by intermediate binary `roi_stats_mag.out' on:
 * 	Thu Oct  7 22:16:04 2010
 *
 * Code for include-generation driver `../Gen-include.c' last modified on:
 * 	Mon Jun  7 15:11:28 2010
 *
 */

// Original documentation block
/*
 * roi_stats_mag: accumulate statistics on regions
 * 
 *  [s,names]=roi_stats_mag(x,y,mag,geom,mode)
 *  * A set of per-region statistics is gathered based on regions
 *  encoded in the image inputs x (containing region tags, 1..Nr, or
 *  0 for no tag) and y (containing region indicators, 0/1).  The
 *  statistics are functions of region configuration, and of
 *  line-of-sight magnetic field mag.
 *  * x must be Nan, or a nonnegative integer.  Pixels in the
 *  range 1..Nr, where Nr is the number of regions, are treated as
 *  tags for being within the named region.
 *  Inputs x of 0 or NaN are treated as not belonging to any
 *  labeled region.
 *  * y must be NaN, 0, or 1 -- 1 indicates activity present.  Not all
 *  tagged pixels will be active; typically the tagged pixels are
 *  large blobs, and the active pixels are finer details within each
 *  blob.
 *  * For each region, a row of statistics is computed:
 *     1: rgnnum = # pixels tagged (all may not be active)
 *     2: rgnsize = projected (flat) area in microhemispheres (0..1e6)
 *     3: rgnarea = un-projected (solid-angle) area in microhemispheres (0..1e6)
 *     4: arnum = # active pixels (x == 1)
 *     5: arsize = projected (flat) active area in microhemispheres (0..1e6)
 *     6: ararea = unprojected (solid-angle) active area in microhemis (0..1e6)
 *     7: arlat,lon = average location, weighted by un-projected area
 *     9: arminlat,lon  lower corner of (lat,lon) bounding box
 *    11: armaxlat,lon  upper corner of (lat,lon) bounding box
 *    13: rgnbtot = sum of absolute LoS flux within the identified region
 *    14: rgnbnet = net LoS flux within the identified region
 *    15: rgnbpos = absolute value of total positive LoS flux
 *    16: rgnbneg = absolute value of total negative LoS flux
 *    17: arfwtlat,lon = flux-weighted center of active pixels
 *    19: arfwtpos_lat,lon = flux-weighted center of positive flux
 *    21: arfwtneg_lat,lon = flux-weighted center of negative flux
 *    23: daysgone = #days until bounding box vanishes from front of disk
 *    24: daysback = #days until bounding box first reappears on front
 *  These definitions are in a machine-readable include file,
 *  roi_stats_mag_def.h
 *  * The solar disk is at location given by center, observed at tip
 *  angle beta, in degrees.
 *  * Optionally returns the standard text name of each of the
 *  computed statistics.
 *  * This is implemented as a MEX file.
 * 
 *  Inputs:
 *    int  x(m,n);
 *    int  y(m,n);
 *    real mag(m,n);
 *    real geom(5);  -- [x0 y0 r_sun b p]
 *    string mode;
 * 
 *  Outputs:
 *    real s(nr,ns)
 *    opt string names(ns)
 * 
 *  See Also:
 * 
 *  turmon june 2010, sep 2010
 * 
 * 
 */

#ifndef _mexfn_roi_stats_mag_h_
#define _mexfn_roi_stats_mag_h_

// function entry point
mexfn_lib_t main_roi_stats_mag;

// argument counts
#define MXT_rsm_NARGIN_MIN 	5
#define MXT_rsm_NARGIN_MAX 	5
#define MXT_rsm_NARGOUT_MIN	0
#define MXT_rsm_NARGOUT_MAX	2

// input argument numbers
#define MXT_rsm_ARG_x	0
#define MXT_rsm_ARG_y	1
#define MXT_rsm_ARG_mag	2
#define MXT_rsm_ARG_geom	3
#define MXT_rsm_ARG_mode	4

// output argument numbers
#define MXT_rsm_ARG_s	0
#define MXT_rsm_ARG_names	1


#endif // _mexfn_roi_stats_mag_h_

// (file ends)
