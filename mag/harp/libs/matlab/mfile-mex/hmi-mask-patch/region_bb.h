
/*
 * THIS IS A GENERATED FILE; DO NOT EDIT.
 *
 * Declarations for calling `region_bb' a.k.a. `rbb' as a C library.
 *
 * Made by intermediate binary region_bb.out on Mon May  3 17:24:41 2010
 *
 * Source file ../Gen-include.c last modified on Mon May  3 17:04:13 2010
 */

// Original documentation block
/*
 * region_bb: find bounding boxes for regions
 * 
 *  [bb]=region_bb(x,coord)
 *  * The list of bounding boxes of regions in the input map x is
 *  returned in bb.
 *  * x must be Nan, or a nonnegative integer.  Pixels in the
 *  range 1..Nr, where Nr is the number of regions, are treated as
 *  belonging to the corresponding region.
 *  Inputs x of 0 or NaN are treated as off-region.
 *  * Each bounding box row is of the form [Lm Ln Um Un] where
 *  L is the lower corner (in both directions) and U is the
 *  upper corner diagonal from L.  The "m" (first coordinate) and
 *  "n" coordinates are given in that order for L and U.  By
 *  convention, the smallest L = (0,0) and the largest U = (m-1,n-1),
 *  but see the coord input.
 *  * If no pixel in a given region numbered between 1..Nr is found,
 *  NaN is returned for that row's bounding box.
 *  * The coord input allows for some coordinate transforms.
 *  coord(1) is an offset, and coord(2) is a block size.  To get
 *  standard Matlab coordinates, use [1 1] to add coord(1)=1 to the
 *  edges.  This is the default: the "conventional" L and U are
 *  adjusted before output.
 *  * The block size is for situations where the pixels were
 *  averaged from a larger image.  In this case the corners are
 *  essentially multiplied by coord(2).  The lowest pixel in the
 *  L block and the highest pixel in the U block are given as the
 *  bounding box coordinates.
 *  * This is implemented as a MEX file.
 * 
 *  Inputs:
 *    real x(m,n); -- 0, NaN, or 1..Nr
 *    opt int coord(1) or coord(2) = [1 1];
 * 
 *  Outputs:
 *    int bb(Nr,4);
 * 
 *  See Also:
 * 
 *  turmon may 2001
 * 
 * 
 */

// function entry point
mexfn_lib_t main_region_bb;

// argument counts
#define MXT_rbb_NARGIN_MIN 	1
#define MXT_rbb_NARGIN_MAX 	2
#define MXT_rbb_NARGOUT_MIN	0
#define MXT_rbb_NARGOUT_MAX	1

// input argument numbers
#define MXT_rbb_ARG_x	0
#define MXT_rbb_ARG_coord	1

// output argument numbers
#define MXT_rbb_ARG_bb	0


// (file ends)
