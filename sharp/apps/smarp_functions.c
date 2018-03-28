
/*===========================================
 
 The following 3 functions calculate the following spaceweather indices:
 
 USFLUX Total unsigned flux in Maxwells
 MEANGBZ Mean value of the vertical field gradient, in Gauss/Mm
 R_VALUE Karel Schrijver's R parameter
 
 The indices are calculated on the pixels in which the conf_disambig segment is greater than 70 and
 pixels in which the bitmap segment is greater than 42. These ranges are selected because the CCD
 coordinate bitmaps are interpolated for certain data (at the time of this CVS submit, all data
 prior to 2013.08.21_17:24:00_TAI contain interpolated bitmaps; data post-2013.08.21_17:24:00_TAI
 contain nearest-neighbor bitmaps).
 
 In the CCD coordinates, this means that we are selecting the pixels that that equal 33 or 34 in bitmap. Here are the definitions of the pixel values:
 
 
 For bitmap:
 1  : weak field outside smooth bounding curve
 2  : strong field outside smooth bounding curve
 33 : weak field inside smooth bounding curve
 34 : strong field inside smooth bounding curve
 
 Written by Monica Bobra 15 August 2012
 Potential Field code (appended) written by Keiji Hayashi
 Error analysis modification 21 October 2013
 
 ===========================================*/
#include <math.h>
#include <mkl.h>

#define PI       (M_PI)
#define MUNAUGHT (0.0000012566370614) /* magnetic constant */

/*===========================================*/

/* Example function 1: Compute total unsigned flux in units of G/cm^2 */

//  To compute the unsigned flux, we simply calculate
//  flux = surface integral [(vector Bz) dot (normal vector)],
//       = surface integral [(magnitude Bz)*(magnitude normal)*(cos theta)].
//  However, since the field is radial, we will assume cos theta = 1.
//  Therefore the pixels only need to be corrected for the projection.

//  To convert G to G*cm^2, simply multiply by the number of square centimeters per pixel.
//  As an order of magnitude estimate, we can assign 0.5 to CDELT1 and 722500m/arcsec to (RSUN_REF/RSUN_OBS).
//  (Gauss/pix^2)(CDELT1)^2(RSUN_REF/RSUN_OBS)^2(100.cm/m)^2
//  =Gauss*cm^2

int computeAbsFlux(float *bz, int *dims, float *absFlux,
                   float *mean_vf_ptr, float *count_mask_ptr,
                   int *bitmask, float cdelt1, double rsun_ref, double rsun_obs)

{
    
    int nx = dims[0];
    int ny = dims[1];
    int i = 0;
    int j = 0;
    int count_mask = 0;
    double sum = 0.0;
    *absFlux = 0.0;
    *mean_vf_ptr = 0.0;
    
     
    if (nx <= 0 || ny <= 0) return 1;
    
	for (i = 0; i < nx; i++)
	{
	   for (j = 0; j < ny; j++)
	   {
	    if ( bitmask[j * nx + i] < 42 ) continue;
            if isnan(bz[j * nx + i]) continue;
            sum += (fabs(bz[j * nx + i]));
            count_mask++;
	   }
	}
    
    *mean_vf_ptr     = sum*cdelt1*cdelt1*(rsun_ref/rsun_obs)*(rsun_ref/rsun_obs)*100.0*100.0;
    *count_mask_ptr  = count_mask;
    //printf("mean_vf_ptr=%f\n",*mean_vf_ptr);
    //printf("count_mask_ptr=%f\n",*count_mask_ptr);
    //printf("cdelt1=%f\n",cdelt1);
    //printf("rsun_ref=%f\n",rsun_ref);
    //printf("rsun_obs=%f\n",rsun_obs);

    return 0;
}

/*===========================================*/
/* Example function 2:  Derivative of B_vertical SQRT( (dBz/dx)^2 + (dBz/dy)^2 ) */

int computeBzderivative(float *bz, int *dims, float *mean_derivative_bz_ptr, int *bitmask, float *derx_bz, float *dery_bz)
{
    
    int nx = dims[0];
    int ny = dims[1];
    int i = 0;
    int j = 0;
    int count_mask = 0;
    double sum = 0.0;
    *mean_derivative_bz_ptr = 0.0;
    
    if (nx <= 0 || ny <= 0) return 1;
    
    /* brute force method of calculating the derivative (no consideration for edges) */
    for (i = 1; i <= nx-2; i++)
    {
	for (j = 0; j <= ny-1; j++)
        {
           derx_bz[j * nx + i] = (bz[j * nx + i+1] - bz[j * nx + i-1])*0.5;
        }
    }
    
    /* brute force method of calculating the derivative (no consideration for edges) */
    for (i = 0; i <= nx-1; i++)
    {
        for (j = 1; j <= ny-2; j++)
        {
           dery_bz[j * nx + i] = (bz[(j+1) * nx + i] - bz[(j-1) * nx + i])*0.5;
        }
    }
    
    /* consider the edges for the arrays that contribute to the variable "sum" in the computation below.
    ignore the edges for the error terms as those arrays have been initialized to zero. 
    this is okay because the error term will ultimately not include the edge pixels as they are selected out by the mask and bitmask arrays.*/
    i=0;
    for (j = 0; j <= ny-1; j++)
    {
        derx_bz[j * nx + i] = ( (-3*bz[j * nx + i]) + (4*bz[j * nx + (i+1)]) - (bz[j * nx + (i+2)]) )*0.5;
    }
    
    i=nx-1;
    for (j = 0; j <= ny-1; j++)
    {
        derx_bz[j * nx + i] = ( (3*bz[j * nx + i]) + (-4*bz[j * nx + (i-1)]) - (-bz[j * nx + (i-2)]) )*0.5;
    }
    
    j=0;
    for (i = 0; i <= nx-1; i++)
    {
        dery_bz[j * nx + i] = ( (-3*bz[j*nx + i]) + (4*bz[(j+1) * nx + i]) - (bz[(j+2) * nx + i]) )*0.5;
    }
    
    j=ny-1;
    for (i = 0; i <= nx-1; i++)
    {
        dery_bz[j * nx + i] = ( (3*bz[j * nx + i]) + (-4*bz[(j-1) * nx + i]) - (-bz[(j-2) * nx + i]) )*0.5;
    }
    
    
    for (i = 0; i <= nx-1; i++)
    {
        for (j = 0; j <= ny-1; j++)
        {
            if ( bitmask[j * nx + i] < 42 ) continue;
            if ( (derx_bz[j * nx + i] + dery_bz[j * nx + i]) == 0) continue;
            if isnan(bz[j * nx + i])      continue;
            if isnan(bz[(j+1) * nx + i])  continue;
            if isnan(bz[(j-1) * nx + i])  continue;
            if isnan(bz[j * nx + i-1])    continue;
            if isnan(bz[j * nx + i+1])    continue;
            if isnan(derx_bz[j * nx + i]) continue;
            if isnan(dery_bz[j * nx + i]) continue;
            sum += sqrt( derx_bz[j * nx + i]*derx_bz[j * nx + i]  + dery_bz[j * nx + i]*dery_bz[j * nx + i] ); /* Units of Gauss */
            count_mask++;
        }
    }
    
    *mean_derivative_bz_ptr = (sum)/(count_mask); // would be divided by ((nx-2)*(ny-2)) if shape of count_mask = shape of magnetogram
    //printf("mean_derivative_bz_ptr=%f\n",*mean_derivative_bz_ptr);
    
	return 0;
}

/*===========================================*/
/* Example function 3: R parameter as defined in Schrijver, 2007 */
//
// Note that there is a restriction on the function fsample()
// If the following occurs:
//      nx_out > floor((ny_in-1)/scale + 1) 
//      ny_out > floor((ny_in-1)/scale + 1),
// where n*_out are the dimensions of the output array and n*_in 
// are the dimensions of the input array, fsample() will usually result 
// in a segfault (though not always, depending on how the segfault was accessed.) 

int computeR(float *los, int *dims, float *Rparam, float cdelt1,
             float *rim, float *p1p0, float *p1n0, float *p1p, float *p1n, float *p1,
             float *pmap, int nx1, int ny1, 
             int scale, float *p1pad, int nxp, int nyp, float *pmapn)

{ 
    int nx = dims[0];
    int ny = dims[1];
    int i = 0;
    int j = 0;
    int index, index1;
    double sum = 0.0;
    *Rparam = 0.0;
    struct fresize_struct fresboxcar, fresgauss;
    struct fint_struct fints;
    float sigma = 10.0/2.3548;
    
    // set up convolution kernels
    init_fresize_boxcar(&fresboxcar,1,1);
    init_fresize_gaussian(&fresgauss,sigma,20,1);

    // =============== [STEP 1] =============== 
    // bin the line-of-sight magnetogram down by a factor of scale 
    fsample(los, rim, nx, ny, nx, nx1, ny1, nx1, scale, 0, 0, 0.0);

    // =============== [STEP 2] =============== 
    // identify positive and negative pixels greater than +/- 150 gauss
    // and label those pixels with a 1.0 in arrays p1p0 and p1n0
    for (i = 0; i < nx1; i++)
    {
        for (j = 0; j < ny1; j++)
        {
            index = j * nx1 + i;
            if (rim[index] > 150)
                p1p0[index]=1.0;
            else
                p1p0[index]=0.0;
            if (rim[index] < -150)
                p1n0[index]=1.0;
            else
                p1n0[index]=0.0;
        }
    }

    // =============== [STEP 3] =============== 
    // smooth each of the negative and positive pixel bitmaps      
    fresize(&fresboxcar, p1p0, p1p, nx1, ny1, nx1, nx1, ny1, nx1, 0, 0, 0.0);
    fresize(&fresboxcar, p1n0, p1n, nx1, ny1, nx1, nx1, ny1, nx1, 0, 0, 0.0);

    // =============== [STEP 4] =============== 
    // find the pixels for which p1p and p1n are both equal to 1. 
    // this defines the polarity inversion line
    for (i = 0; i < nx1; i++)
    {
        for (j = 0; j < ny1; j++)
        {
            index = j * nx1 + i;
            if ((p1p[index] > 0.0) && (p1n[index] > 0.0))
                p1[index]=1.0;
            else
                p1[index]=0.0;
        }
    }

    // pad p1 with zeroes so that the gaussian colvolution in step 5
    // does not cut off data within hwidth of the edge
   
    // step i: zero p1pad
    for (i = 0; i < nxp; i++)
    {
        for (j = 0; j < nyp; j++)
        {
            index = j * nxp + i;
            p1pad[index]=0.0;
        }
    }

    // step ii: place p1 at the center of p1pad
    for (i = 0; i < nx1; i++)
    {
       for (j = 0; j < ny1; j++)
       {
            index  = j * nx1 + i; 
            index1 = (j+20) * nxp + (i+20);
            p1pad[index1]=p1[index];
        }
    }

    // =============== [STEP 5] =============== 
    // convolve the polarity inversion line map with a gaussian
    // to identify the region near the plarity inversion line
    // the resultant array is called pmap
    fresize(&fresgauss, p1pad, pmap, nxp, nyp, nxp, nxp, nyp, nxp, 0, 0, 0.0);


   // select out the nx1 x ny1 non-padded array  within pmap
    for (i = 0; i < nx1; i++)
    {
       for (j = 0; j < ny1; j++)
       {
            index  = j * nx1 + i; 
            index1 = (j+20) * nxp + (i+20);
            pmapn[index]=pmap[index1];
        }
    }

    // =============== [STEP 6] =============== 
    // the R parameter is calculated
    for (i = 0; i < nx1; i++)
    {
        for (j = 0; j < ny1; j++)
        {
            index = j * nx1 + i;
            if isnan(pmapn[index]) continue;
            if isnan(rim[index]) continue;
            sum += pmapn[index]*abs(rim[index]);
        }
    }
    
    if (sum < 1.0)
        *Rparam = 0.0;
    else
        *Rparam = log10(sum);

    //printf("R_VALUE=%f\n",*Rparam);

    free_fresize(&fresboxcar);
    free_fresize(&fresgauss);

    return 0;

}

/*===========================================*/

char *smarp_functions_version() // Returns CVS version of smarp_functions.c
{
    return strdup("$Id");
}
