/*
 *  soi_fun.h						~soi/(version)/include
 *
 *  Prototype decalarations required for SOI functions.
 *
 *  Source files are on directory ~soi/(version)/src/functions unless
 *    otherwise noted.
 *  Additional information is in the following man pages:
 *
 *  Responsible:  Rick Bogart			RBogart@solar.Stanford.EDU
 *
 *  Bugs:
 *
 *  Revision history is at the end of the file.
 */
#ifndef SOI_FUN_INCL
/****************************************************************************/
/**************************  INCLUDE STATEMENTS  ****************************/
/****************************************************************************/

#ifndef SOI_VERSION_INCL
#include "soi_version.h"
#endif

#ifndef SOI_SDS_INCL
#include "soi_sds.h"
#endif

/****************************************************************************/
/****************************  DEFINE STATEMENTS  ***************************/
/****************************************************************************/

#define BILINEAR		(0)
#define CUBIC_CONVOLUTION	(1)

#define VCOR_SIGN_ONLY		(0)

/* TBD by Schou or ??? - document meaning and supply more descriptive names */
#define VCOR_LEVEL1		(1)
#define VCOR_LEVEL2		(2)

#define MCOR_LEVEL0		(0)
							  /*  nocorrection  */
#define MCOR_LEVEL1		(1)
/*  line-of-sight correction for Br assuming observer at infinite distance  */

#define TILE_FMT        "V_pox_%d-%d_%d%c"
/*  example tile names:   V_pox_0-49_01h  V_pox_0-100_04h  */

#define LMAX            1500
#define MODES(L)        ((((L)+1)*(L))/2)
#define MINIMUM(X,Y)    (((X)<(Y))?(X):(Y))
#define MAXIMUM(X,Y)    (((X)>(Y))?(X):(Y))
#define TNEXT(X)        (((X)->tstart)+((X)->twidth))
#define LNEXT(X)        (((X)->lstart)+((X)->lwidth))

#define PLATFORM_UNKNOWN	(0)
#define PLATFORM_UNRECOGNIZED	(1)
#define PLATFORM_SOHO		(2)
#define PLATFORM_GONGPLUS	(3)
#define PLATFORM_MWO60		(4)
#define PLATFORM_BBSO		(5)
#define PLATFORM_TRACE		(6)
#define PLATFORM_SPOLE_JSO	(10)
#define PLATFORM_GONG		(30)
#define PLATFORM_OBSPM		(40)

#define INSTRUMENT_UNKNOWN	(0)
#define INSTRUMENT_UNRECOGNIZED	(1)
#define INSTRUMENT_SOHO_MDI	(10)
#define INSTRUMENT_SOHO_EIT	(11)
#define INSTRUMENT_GONG_TD	(20)
#define INSTRUMENT_GONG_CT	(21)
#define INSTRUMENT_GONG_TC	(22)
#define INSTRUMENT_GONG_BB	(23)
#define INSTRUMENT_GONG_ML	(24)
#define INSTRUMENT_GONG_LE	(25)
#define INSTRUMENT_GONG_UD	(26)
#define INSTRUMENT_GONG_MERGE	(29)
#define INSTRUMENT_MWO60_MOF	(30)
#define INSTRUMENT_BBSO_SINGER	(40)
#define INSTRUMENT_TRACE	(50)
#define INSTRUMENT_MOTH		(60)
#define INSTRUMENT_OBSPM_SPHG	(70)

#define NO_DATA_DICT	(0x0001)
#define NO_SEMIDIAMETER	(0x0002)
#define NO_XSCALE	(0x0004)
#define NO_YSCALE	(0x0008)
#define NO_XCENTERLOC	(0x0010)
#define NO_YCENTERLOC	(0x0020)
#define NO_HELIO_LATC	(0x0040)
#define NO_HELIO_LONC	(0x0080)
#define NO_HELIO_PA	(0x0100)
#define NO_OBSERVER_LAT	(0x0002)
#define NO_OBSERVER_LON	(0x0004)

#define SOI_FUN_VERSION_NUM	(4.8)
#define SOI_FUN_INCL	1

/****************************************************************************/
/*******************************  TYPEDEFS  *********************************/
/****************************************************************************/

typedef struct tile_info {
   struct tile_info *next;      /* just in case we want to make a list */
   char *name;  /* name of tile is series name */
   int  index;  /* index of tile is series number */
   int  tstart; /* start time in minutes */
   int  twidth; /* width of tile in minutes */
   int  lstart; /* start value of l */
   int  lwidth; /* width/height of tile in l */
   int  modes;  /* modes in lwidth starting at lstart */
} TILE;

typedef struct gradinfo {
   float a[3];
   int np;
   double p;
   double e;
   double V;
} GRADINFO;


/****************************************************************************/
/***************************  FUNCTION PROTOTYPES  **************************/
/****************************************************************************/

						/* source file: apodize.c */
/*  Note:  If this is a general function it belongs in lib_M.d  */
/*  Probably the input array should be in an SDS  */
/* extern int apodize(float *data, int cols, int rows, int apod); */
extern int apodize(
float *data,		/* input/output data array */
double B0,		/* heliographic latitude of disk center */
int cols, int rows,	/* width and height of input array */
double Lmin,		/* start of longitude range */
double Ldelta,		/* increment in longitude pre col */
double sinBdelta,	/* increment in sin latitude per row */
int apodlevel,		/* type of apodization */
			/* 0 to do no apodization */
			/* 1 to apodize in true solar coordinates */
			/* 2 to apodize in ideal solar coordinates */
double apinner,		/* fractional radius to start apodization at */
double apwidth,		/* width of apodization */
int apel,		/* do elliptical apodization as described by */
                        /* apx and apy */
double apx,		/* divide the x position by this before applying */
			/* apodization */
double apy		/* divide the y position by this before applying */
			/* apodization */
);

					     /*  source file: data_stats.c  */
extern int data_stats	(SDS *in, SDS_STATS *in_out);
extern void fget_stats (float *data, int ntot, int *nv, double *min,
    double *max, double *mean, double *stdv, double *skew, double *kurt,
    double *med, double *min_abs);

					       /*  source file: keywords.c  */
extern double GONG_correct_pa (SDS *img);
	/*  note: the above belongs in soi_GONG.h  */
extern int observer_loc_params (SDS *rec, double *obs_lat, double *obs_lon);
extern int observing_source (SDS *rec, int *platform, int *instrument);
extern TIME observing_time (SDS *rec);
extern int solar_image_params (SDS *rec, double *xscl, double *yscl,
    double *ctrx, double *ctry, double *latc, double *lonc,
    double *apsd, double *pang, double *eecc, double *eang,
    int *xinv, int *yinv);
extern char *platform_name (int platform);
extern char *instrument_name (int instrument);

						/*  source file: masking.c  */
extern int ordered_mask_compare (unsigned int a, unsigned int b,
    int nmasks, ...);

					     /*  source file: matrix_sum.c  */
extern int matrix_sum (double a, SDS *am, double b, SDS *bm, SDS *cm);


typedef struct RotRate_struct
{
  float lat;
  float r;
} RotRate_t;
					      /*  source file: obs2helio.c  */
extern int obs2helio (
	float   *V,      /* Input velocity array                         */
                 	/* assumed declared V[ypixels][xpixels]          */
	float   *U,      /* Output projected array                       */
	int     xpixels, /* x width of input array                       */
	int     ypixels, /* y width of input array                       */
	double  x0,	/* x pixel address of disk center               */
	double  y0,	/* y pixel address of disk center               */
	double  B0,      /* Heliographic Latitude of disk center         */
	double  P,       /* Angle between CCD y-axis and solar vertical  */
                 	/* positive to the east (CCW)                    */
	double  S,       /* Subtended radius of the sun (in radians)     */
	double  rsun,	/* pixels */
	double  Rmax,    /* Maximum disk radius to use (e.g. 0.95)       */

	int	interpolation,
	int     cols,    /* width of output array                        */
	int     rows,    /* height of output array                       */
	double  Lmin,    /* start of longitude range                     */
	double  Ldelta,  /* increment of longitude per col               */
	double  Ladjust,  
	double  sinBdelta, /* increment of sin latitude per row 	*/ 

	double smajor, double sminor, double sangle, double xscale,
        double yscale, char *orientation, int mag_correction,
	int velocity_correction,
	double obs_vr, double obs_vw, double obs_vn, double vsign,
	int NaN_beyond_rmax,
	int carrStretch,
	float   diffrotA,
	float   diffrotB,
	float   diffrotC,
	RotRate_t *rRates,
	int        size);

						   /*  source file: plm1.f  */
extern int setplm_ (int	*lmin, int *lmax, int *m, int *nx, int *indx, 
	double *x, int *nplm, double *plm);

					/*  source file: sds_bin_extract.c  */
extern int sds_bin_average (SDS *sds, int *binsize);
extern int sds_bin_sum (SDS *sds, int *binsize);
extern int sds_extract_subset (SDS *sds, int *offset, int *size);
extern SDS *SDS_bin_average (SDS *sds, int *binsize);
extern SDS *SDS_bin_sum (SDS *sds, int *binsize);

					     /*  source file: sds_interp.c  */
extern int sds_regrid (SDS *sds, int *length, int scheme);
float SDSimaginterp (SDS *image, double x, double y);

					  /*  source file: sds_rowfft842.c  */
extern int sds_rowfft842 ( int isign, SDS * real, SDS * imag);

					  /*  source file: shc_inner_pro.c  */
extern int shc_inner_product (
	SDS *even_real,		      /*  Folded row FFTs of remapped data  */
	SDS *even_imag,		      /*  Even parity                       */
	SDS *odd_real, SDS *odd_imag, /*  Odd parity              	    */
	SDS *masks, int l_min, int l_max, int delta_m,
	double d_m,			  /*  ratio of column of data to m  */
	SDS *i_p_real, SDS *i_p_imag);

						   /*  source file: tile.c  */
extern int lchunks (int lwidth);
extern int mchunks (int mwidth);
extern int linfo (int lindex, int mwidth, int *start, int *width);
extern TILE *get_tile (KEY *params, char *root, int tstep);
extern TILE *intersect_tiles (TILE *tile1, TILE *tile2);
extern int sds_calc_stats (SDS *sds);

#endif

/*
 *  Revision History
 *  v 1.5  96.08.07	R Bogart	added prototype for fget_stats()
 *  v 1.6  96.09.10	R Bogart	added min_abs argument to fget_stats()
 *  v 2.0  96.11.21	R Bogart	removed sds_gongprep declaration
 *	   96.11.22	R Bogart	added declarations from sds_bin_extract
 *	   96.11.27	R Bogart	added temporary declaration for
 *	mdi_bin_avg ()
 *  v 2.1  96.12.13	R Bogart	added declarations from sds_bin_extract
 *	for sds_bin_sum, SDS_bin_average, SDS_bin_sum; removed temporary
 *	   97.04.16	K Chu		added include of <soi_version.h>
 *	   97.05.24	P Scherrer	removed show_tile() spec
 *	   98.02.04	K Chu		added NaN_beyond_rmax to obs2helio()
 *	parameter list
 *  v 2.10 98.03.25	R Bogart	added declaration from masking; removed
 *	history prior to v.1; removed unused sds_remap and obsolete
 *	sds_makemasks declarations and declarations of functions prototyped
 *	in MDI.h
 *	   99.04.27	J Schou		added GONG style elliptical apodization
 *	   99.12.09	K Chu		added argument rm_limbdark to
 *	obs2helio()
 *	   99.12.16	J Aloise	removed limbdark from obs2helio()
 *	   00.05.15	K Chu		added mag_correction to obs2helio()
 *  v 4.8  01.12.06	R Bogart	added declarations for keywords.c
 *	02.09.10	R Bogart	added more "
 *	05.01.05	R Bogart	added more "
 *	05.09.19	R Bogart	added more "
 */

/*
$Id: soi_fun.h,v 1.1 2009/04/24 21:52:49 production Exp $
$Source: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/jpe/apps/soi_fun.h,v $
$Author: production $
*/
/*
 * $Log: soi_fun.h,v $
 * Revision 1.1  2009/04/24 21:52:49  production
 * *** empty log message ***
 *
 * Revision 1.36  2007/09/24  16:22:27  arta
 * Add structure definition to hold rotation rate v. latitude data, change signature of obs2helio to pass rotation rate data back to caller.
 *
 * Revision 1.35  2007/08/22 18:15:04  arta
 * Fix the build - wrong declaration of obs2helio(), also missing semicolon in obs2helio.c
 *
 * Revision 1.34  2007/02/02 22:11:15  arta
 * Adjust longitude for differential rotation.
 *
 * Revision 1.33  2005/09/19 23:21:20  rick
 * added Meudon identifiers
 *
 * Revision 1.32  2005/01/06 00:44:38  rick
 * see above
 *
 * Revision 1.31  2004/06/11 18:43:22  rick
 * added MOTH identifiers
 *
 * Revision 1.30  2003/11/06 01:45:02  rick
 * added argument to sds_regrid
 *
 * Revision 1.29  2003/03/03 18:51:24  rick
 * see above
 *
 * Revision 1.28  2002/10/08 17:41:53  rick
 * see above
 *
 * Revision 1.27  2002/05/24  00:09:29  rick
 * see above
 *
 * Revision 1.26  2002/02/21  20:04:05  rick
 * see above
 *
 * Revision 1.25  2001/12/10  18:49:16  rick
 * see above
 *
 * Revision 1.24  2000/05/15  17:19:12  kehcheng
 * added mag_correction to obs2helio()
 *
 * Revision 1.23  1999/12/16 19:44:41  CM
 * remove limbdark in obs2helio
 *
 * Revision 1.22  1999/12/09  23:26:21  kehcheng
 * added argument rm_limbdark to obs2helio()
 *
 * Revision 1.21  1999/04/27 18:43:40  schou
 * Added GONG style elliptical apodization.
 *
 * Revision 1.20  1998/04/02  22:02:26  rick
 * see above
 *
 * Revision 1.19  1998/02/24  18:45:40  kehcheng
 * added NaN_beyond_rmax to obs2helio() parameter list
 *
 * Revision 1.18  1997/05/24 01:28:14  phil
 * removed show_tile spec
 *
 * Revision 1.17  1997/04/16  21:52:35  kehcheng
 * added #include <soi_version.h>
 *
 * Revision 1.16  1997/01/13  18:44:56  rick
 * see above
 *
 * Revision 1.15  1996/11/27  19:05:17  rick
 * as above
 *
 * Revision 1.14  1996/11/22  22:25:07  rick
 * see above
 *
 * Revision 1.13  1996/11/02  00:48:48  phil
 * added sds_calc_stats
 *
 * Revision 1.12  1996/09/12  09:49:18  rick
 * see above
 *
 */
