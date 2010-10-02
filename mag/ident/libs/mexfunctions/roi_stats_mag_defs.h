/*
 * roi_stats_mag_defs: definitions for ROI statistics
 *
 * This list was a starting point determined by Michael Turmon, 
 * Todd Hoeksema, Xudong Sun, and others, in May 2010.
 *
 * intended to be included from various source files, so, 
 */

#include <limits.h>  // for INT_MAX

#ifndef ROI_STATS_MAG_DEFS_H
#define ROI_STATS_MAG_DEFS_H 1

typedef enum {
  RS_rgnnum = 0,  // # pixels tagged (all may not be active)
  RS_rgnsize,     // projected (flat) area in microhemispheres (0..1e6)
  RS_rgnarea,     // un-projected (solid-angle) area in microhemispheres (0..1e6)
  RS_arnum,       // # active pixels
  RS_arsize,      // projected (flat) active area in microhemispheres (0..1e6)
  RS_ararea,      // unprojected (solid-angle) active area in microhemis (0..1e6)
  RS_arlat,       // average location, weighted by un-projected area
  RS_arlon,
  RS_arminlat,    // lower corner of (lat,lon) bounding box
  RS_arminlon,
  RS_armaxlat,    // upper corner of (lat,lon) bounding box
  RS_armaxlon,
  RS_rgnbtot,     // sum of absolute LoS flux within the identified region
  RS_rgnbnet,     // net LoS flux within the identified region
  RS_rgnbpos,     // absolute value of total positive LoS flux
  RS_rgnbneg,     // absolute value of total negative LoS flux
  RS_arfwtlat,    // flux-weighted center of active pixels
  RS_arfwtlon, 
  RS_arfwtposlat, // flux-weighted center of positive flux
  RS_arfwtposlon, 
  RS_arfwtneglat, // flux-weighted center of negative flux
  RS_arfwtneglon,
  RS_daysgone,    // #days until bounding box vanishes from disk front
  RS_daysback,    // #days until bounding box first reappears on front
  RS_num_stats,   // number of statistics
  RS_MAX = INT_MAX  // ensure it's as big as an int
  } rs_stat_index_t;

// autogenerated from the above list with shell-command-on-region using:
// sed -e 's/RS_/\"/' -e 's/,.*/\",/'

static const char *RS_index2name[] = {
  "rgnnum",
  "rgnsize",
  "rgnarea",
  "arnum",
  "arsize",
  "ararea",
  "arlat",
  "arlon",
  "arminlat",
  "arminlon",
  "armaxlat",
  "armaxlon",
  "rgnbtot",
  "rgnbnet",
  "rgnbpos",
  "rgnbneg",
  "arfwtlat",
  "arfwtlon",
  "arfwtposlat",
  "arfwtposlon",
  "arfwtneglat",
  "arfwtneglon",
  "daysgone",
  "daysback",
  NULL,
  NULL,
};

// correspondence of statistics above to HMI keywords
// NULL means, no such HMI keyword
// first char gives type of HMI keyword

static const char *RS_index2keyname[] = {
  // region size
  "iNPIX",       // RS_rgnnum
  "fSIZE",       // RS_rgnsize
  "fAREA",       // RS_rgnarea    
  // ar size
  "iNACR",       // RS_arnum
  "fSIZE_ACR",   // RS_arsize     
  "fAREA_ACR",   // RS_ararea     
  // ar extent
  NULL,          // RS_arlat      
  NULL,          // RS_arlon	 
  NULL,          // RS_arminlat   
  NULL,          // RS_arminlon	 
  NULL,          // RS_armaxlat   
  NULL,          // RS_armaxlon	 
  // flux
  "fBTOT",       // RS_rgnbtot    
  "fBNET",       // RS_rgnbnet    
  "fBPOS_TOT",   // RS_rgnbpos    
  "fBNEG_TOT",   // RS_rgnbneg    
  // flux-weighted extent
  "fFWTLAT",     // RS_arfwtlat   
  "fFWTLON",     // RS_arfwtlon 	 
  "fFWTPOS_LAT", // RS_arfwtposlat
  "fFWTPOS_LON", // RS_arfwtposlon
  "fFWTNEG_LAT", // RS_arfwtneglat
  "fFWTNEG_LON", // RS_arfwtneglon
  NULL,          // RS_daysgone
  NULL,          // RS_daysback
  // end
};

#endif // include once
