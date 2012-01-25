/*
 * roi_stats_mag_defs: definitions for ROI statistics
 *
 * This list was a starting point determined by Michael Turmon, 
 * Todd Hoeksema, Xudong Sun, and others, in May 2010.
 * Updated October 2010.
 *
 * Intended to be included from various source files, so,
 * we use include guards to prevent double inclusion.
 */

#include <limits.h>  // for INT_MAX

#ifndef ROI_STATS_MAG_DEFS_H
#define ROI_STATS_MAG_DEFS_H 1

/*
 * RS_rgn_* are whole-region accumulators (all over the ROI)
 * RS_ar_*  are ar-only accumulators
 * the convention is to precede lat and lon with `_'
 */

typedef enum {
  // size
  RS_rgn_num = 0,   // # pixels tagged (all may not be active)
  RS_rgn_size,      // projected (flat) area in microhemispheres (0..1e6)
  RS_rgn_area,      // un-projected (solid-angle) area in microhemispheres (0..1e6)
  // extent
  RS_rgn_min_lat,   // lower corner of (lat,lon) bounding box for region
  RS_rgn_min_lon,
  RS_rgn_max_lat,   // upper corner of (lat,lon) bounding box for region
  RS_rgn_max_lon,
  // return time
  RS_rgn_daysgone,  // #days until ROI bounding box vanishes from disk front
  RS_rgn_daysback,  // #days until ROI bounding box first reappears on front
  // flux
  RS_rgn_btot,      // sum of absolute LoS flux within the identified region
  RS_rgn_bnet,      // net LoS flux within the identified region
  RS_rgn_bpos,      // absolute value of total positive ROI LoS flux
  RS_rgn_bneg,      // absolute value of total negative ROI LoS flux
  // flux moments
  RS_rgn_bsum1,     // sum of  LoS flux    within the identified region (== bnet)
  RS_rgn_bsum2,     // sum of (LoS flux)^2 within the identified region
  RS_rgn_bsum3,     // sum of (LoS flux)^3 within the identified region
  RS_rgn_bsum4,     // sum of (LoS flux)^4 within the identified region
  // flux moments, standardized
  RS_rgn_bmean,     // mean of LoS flux within the region
  RS_rgn_bsdev,     // standard deviation of LoS flux within the region
  RS_rgn_bskew,     // skewness of LoS flux within the region
  RS_rgn_bkurt,     // kurtosis (subtracting 3) of LoS flux within the region
  // size
  RS_ar_num,        // # active pixels
  RS_ar_size,       // projected (flat) active area in microhemispheres (0..1e6)
  RS_ar_area,       // unprojected (solid-angle) active area in microhemis (0..1e6)
  // flux
  RS_ar_btot,       // sum of absolute LoS flux within the active region
  RS_ar_bnet,       // net LoS flux within the active region
  RS_ar_bpos,       // absolute value of total positive AR LoS flux
  RS_ar_bneg,       // absolute value of total negative AR LoS flux
  // mean location
  RS_ar_area_lat,   // area-weighted center of active pixels (norm: RS_ar_area)
  RS_ar_area_lon,
  RS_ar_fwt_lat,    // flux-weighted center of active pixels (norm: RS_ar_btot)
  RS_ar_fwt_lon, 
  RS_ar_fwtpos_lat, // flux-weighted center of positive AR flux (norm: RS_ar_bpos)
  RS_ar_fwtpos_lon,
  RS_ar_fwtneg_lat, // flux-weighted center of negative AR flux (norm: RS_ar_bneg)
  RS_ar_fwtneg_lon,
  // bookkeeping
  RS_num_stats,     // number of statistics
  RS_MAX = INT_MAX  // ensure it's as big as an int
  } rs_stat_index_t;

// autogenerated from the above list with shell-command-on-region using:
// sed -e 's/RS_/\"/' -e 's/,.*/\",/'

static const char *RS_index2name[] = {
  // size
  "rgn_num",
  "rgn_size",
  "rgn_area",
  // extent
  "rgn_min_lat",
  "rgn_min_lon",
  "rgn_max_lat",
  "rgn_max_lon",
  // return time
  "rgn_daysgone",
  "rgn_daysback",
  // flux
  "rgn_btot",
  "rgn_bnet",
  "rgn_bpos",
  "rgn_bneg",
  // flux moments
  "rgn_bsum1",
  "rgn_bsum2",
  "rgn_bsum3",
  "rgn_bsum4",
  // flux moments, standardized
  "rgn_bmean",
  "rgn_bsdev",
  "rgn_bskew",
  "rgn_bkurt",
  // size
  "ar_num",
  "ar_size",
  "ar_area",
  // flux
  "ar_btot",
  "ar_bnet",
  "ar_bpos",
  "ar_bneg",
  // mean location
  "ar_area_lat",
  "ar_area_lon",
  "ar_fwt_lat",
  "ar_fwt_lon",
  "ar_fwtpos_lat",
  "ar_fwtpos_lon",
  "ar_fwtneg_lat",
  "ar_fwtneg_lon",
  // bookkeeping
  NULL,
  NULL,
};

// how to combine two or more summary statistics
// into a pooled statistic: sum, min, max, or
// weighted average

static const char *RS_index2combo[] = {
  // size
  "sum",          // rgn_num
  "sum",          // rgn_size
  "sum",          // rgn_area
  // extent
  "min",          // rgn_min_lat
  "min",          // rgn_min_lon
  "max",          // rgn_max_lat
  "max",          // rgn_max_lon
  // return time
  "max",          // rgn_daysgone
  "min",          // rgn_daysback
  // flux
  "sum",          // rgn_btot
  "sum",          // rgn_bnet
  "sum",          // rgn_bpos
  "sum",          // rgn_bneg
  // flux moments
  "sum",          // rgn_bsum1
  "sum",          // rgn_bsum2
  "sum",          // rgn_bsum3
  "sum",          // rgn_bsum4
  // flux moments, standardized
  "recomp",       // rgn_bmean
  "recomp",       // rgn_bsdev
  "recomp",       // rgn_bskew
  "recomp",       // rgn_bkurt
  // size
  "sum",          // ar_num
  "sum",          // ar_size
  "sum",          // ar_area
  // flux
  "sum",          // ar_btot
  "sum",          // ar_bnet
  "sum",          // ar_bpos
  "sum",          // ar_bneg
  // mean location
  "avg.ar_area",  // ar_area_lat
  "avg.ar_area",  // ar_area_lon
  "avg.ar_btot",  // ar_fwt_lat
  "avg.ar_btot",  // ar_fwt_lon
  "avg.ar_bpos",  // ar_fwtpos_lat
  "avg.ar_bpos",  // ar_fwtpos_lon
  "avg.ar_bneg",  // ar_fwtneg_lat
  "avg.ar_bneg",  // ar_fwtneg_lon
  // bookkeeping
  NULL,
  NULL,
};

// correspondence of statistics above to HMI "patch" keywords
// NULL means, no such keyword in patch series
// first char gives type of HMI keyword

static const char *RS_index2patch_keyname[] = {
  // size
  "iNPIX",        // RS_rgn_num
  "fSIZE",        // RS_rgn_size
  "fAREA",        // RS_rgn_area    
  // extent
  "fMIN_LAT",     // RS_rgn_min_lat   
  "fMIN_LON",     // RS_rgn_min_lon	 
  "fMAX_LAT",     // RS_rgn_max_lat   
  "fMAX_LON",     // RS_rgn_max_lon	 
  // return time
  NULL,           // RS_rgn_daysgone
  NULL,           // RS_rgn_daysback
  // flux
  "fMTOT",        // RS_rgn_btot    
  "fMNET",        // RS_rgn_bnet    
  "fMPOS_TOT",    // RS_rgn_bpos    
  "fMNEG_TOT",    // RS_rgn_bneg    
  // flux moments
  NULL,           // RS_rgn_bsum1
  NULL,           // RS_rgn_bsum2
  NULL,           // RS_rgn_bsum3
  NULL,           // RS_rgn_bsum4
  // flux moments, standardized
  "fMMEAN",       // RS_rgn_bmean
  "fMSTDEV",      // RS_rgn_bsdev
  "fMSKEW",       // RS_rgn_bskew
  "fMKURT",       // RS_rgn_bkurt
  // size
  "iNACR",        // RS_ar_num
  "fSIZE_ACR",    // RS_ar_size     
  "fAREA_ACR",    // RS_ar_area     
  // flux
  NULL,           // RS_ar_btot    
  NULL,           // RS_ar_bnet    
  NULL,           // RS_ar_bpos    
  NULL,           // RS_ar_bneg    
  // mean location
  NULL,           // RS_ar_area_lat
  NULL,           // RS_ar_area_lon	 
  "fFWT_LAT",     // RS_ar_fwt_lat   
  "fFWT_LON",     // RS_ar_fwt_lon 	 
  "fFWTPOS_LAT",  // RS_ar_fwtpos_lat
  "fFWTPOS_LON",  // RS_ar_fwtpos_lon
  "fFWTNEG_LAT",  // RS_ar_fwtneg_lat
  "fFWTNEG_LON",  // RS_ar_fwtneg_lon
  // end
};

// correspondence of statistics above to HMI "mask" keywords
// NULL means, no such keyword in mask series
// first char gives type of HMI keyword

static const char *RS_index2mask_keyname[] = {
  // size
  "iAR_NPIX",     // RS_rgn_num
  "fAR_SIZE",     // RS_rgn_size
  "fAR_AREA",     // RS_rgn_area    
  // extent
  NULL,           // RS_rgn_min_lat   
  NULL,           // RS_rgn_min_lon	 
  NULL,           // RS_rgn_max_lat   
  NULL,           // RS_rgn_max_lon	 
  // return time
  NULL,           // RS_rgn_daysgone
  NULL,           // RS_rgn_daysback
  // flux
  "fAR_MTOT",     // RS_rgn_btot    
  "fAR_MNET",     // RS_rgn_bnet    
  "fAR_MPOS",     // RS_rgn_bpos    
  "fAR_MNEG",     // RS_rgn_bneg    
  // flux moments
  NULL,           // RS_rgn_bsum1
  NULL,           // RS_rgn_bsum2
  NULL,           // RS_rgn_bsum3
  NULL,           // RS_rgn_bsum4
  // flux moments, standardized
  "fAR_MMEAN",    // RS_rgn_bmean
  "fAR_MSDEV",    // RS_rgn_bsdev
  "fAR_MSKEW",    // RS_rgn_bskew
  "fAR_MKURT",    // RS_rgn_bkurt
  // size
  NULL,           // RS_ar_num
  NULL,           // RS_ar_size     
  NULL,           // RS_ar_area     
  // flux
  NULL,           // RS_ar_btot    
  NULL,           // RS_ar_bnet    
  NULL,           // RS_ar_bpos    
  NULL,           // RS_ar_bneg    
  // mean location
  NULL,           // RS_ar_area_lat
  NULL,           // RS_ar_area_lon	 
  NULL,           // RS_ar_fwt_lat   
  NULL,           // RS_ar_fwt_lon 	 
  NULL,           // RS_ar_fwtpos_lat
  NULL,           // RS_ar_fwtpos_lon
  NULL,           // RS_ar_fwtneg_lat
  NULL,           // RS_ar_fwtneg_lon
  // end
};

#endif // include once
