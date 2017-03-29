/*
 * FLCT
 * Written byFLCT: http://solarmuri.ssl.berkeley.edu/overview/publicdownloads/software.html
 * Copyright (C) 2007-2009 Regents of the University of California
 * Modified from v1.01-1 for HMI/JSOC pipeline by X. Sun
 * Input:
 *  dims[2]:            array dimensions
 *  *bz0, *bz1:         input arrays (dims[0]*dims[1])
 *  deltas, deltast, sigma:
                        pixel size, time diff, window size
 *  useThresh/doFilter: flag for threshold/filtering
 *  threshb/kr:          threshold, filter wavelength
 * Output:
 *  vx/vy/vm:           velocity x/y, mask
 *
 */

# include <stdio.h>
# include <string.h>
# include <stdlib.h>
# include <math.h>
# include <complex.h>
# include <fftw3.h>

typedef int i4;
typedef float f4;

// Prototypes

i4 cross_cor (i4 init, i4 hires, i4 expand, double *arr, double *barr,
              double **absccor, i4 nx, i4 ny, double *shiftx, double *shifty,
              i4 filterflag, double kr);
i4 interpcc2d (double *fdata, double xmiss, i4 nx, i4 ny,
               double *xwant, i4 nxinterp, double *ywant, i4 nyinterp, double **finterp);
i4 gaussfilt(double *filter, double *kx, double *ky, i4 nx, i4 ny, double kr);
i4 shift2d (double *arr, i4 nx, i4 ny, i4 ishift, i4 jshift);
i4 maxloc (double *arr, i4 xsize);
i4 minloc (double *arr, i4 xsize);
i4 make_freq(double *k, i4 ndim);
i4 is_big_endian ();


// Main function
// Skip is not implemented

int flct4jsoc(int *dims, double *bz0, double *bz1,
              double deltas, double deltat, double sigma,
              int useThresh, double threshb, int doFilter, double kr,
              int quiet, int absthresh,
              double *vx, double *vy, double *vm)
{
    int status = 0;
    char *version ="1.01   ";
 
    // Input parameters
    // Syntax: %s ifile ofile deltat deltas sigma -t thr -k kr -s N[pP][qQ][i] -h -q\n\n"
    // deltat, deltas, sigma, kr, quiet are passed in
    
    double deltinv = 1. / deltat;
    i4 sigmaeq0 = (sigma > 0.) ? 0 : 1;
    double sigminv = (sigmaeq0) ? -999. : (1. / sigma);
    i4 hires = -1, verbose = ((quiet) ? 0 : 1), expand = 1;
    double thresh = (useThresh) ? fabs(threshb) : 0.;
    i4 filter = doFilter;
    i4 poffset = 0, qoffset = 0, skip = 0, skipon = 0;      // No skip info for now
    i4 skipxy, noskipx, noskipy, noskipxy;
    i4 xoffset = 0, yoffset = 0, interpolate = 0;
    skipon = skip + abs(qoffset) + abs(poffset);
    i4 ibe = is_big_endian();
    char *aloc = NULL;
    if (absthresh) aloc = (char *) malloc(sizeof(char));
    char *ploc = NULL, *qloc = NULL, *intloc = NULL;
    double xmiss = -999999.;
    
    if (verbose) {          // Print some info
        printf("flct: Version %s Copyright: 2007,2008 University of California\n",
               version);
        if (ibe)
            printf ("flct: large endian machine; i/o not byteswapped\n");
        else
            printf ("flct: small endian machine; i/o will be byteswapped\n");
        printf ("flct: deltat = %g\n", deltat);
        printf ("flct: deltas = %g\n", deltas);
        printf ("flct: sigma = %g\n", sigma);
        printf ("flct: threshold image value for LCT is %g\n", thresh);
        if (aloc)
            printf ("flct: threshold forced to be in abs. units\n");
        if (skipon)
            printf ("flct: skip = %d pixels with p=%d, q=%d\n",skip, poffset, qoffset);
        if (skipon && ((poffset < 0) || (qoffset < 0)))
            printf ("flct: p=%d, q=%d: negative p,q will be reset to skip-|p,q|\n",
                    poffset, qoffset);
        if (interpolate)
            printf ("flct: skipped pixels interpolated with cubic convolution\n");
        if (filter)
            printf ("flct: filter rolloff value for input images is %g\n", kr);
        if (hires == 0) printf ("flct: hires option turned on\n");
    }
    if (poffset < 0) poffset = skip - abs(poffset);
    if (qoffset < 0) qoffset = skip - abs(qoffset);
    
    // Input arrays
    
    i4 ier = 0, transp = 1;     // nonzero to transpose input/output arrays
    int nx = dims[1], ny = dims[0], nxny = nx * ny;     // nx: #rows, ny: #columns, tested
    i4 nxorig = nx, nyorig = ny;
    double *f1 = bz0, *f2 = bz1;
    
    if (verbose) {      // Print some info
        printf ("flct: from input file, nx = %d, ny = %d\n", nx, ny);
        if ((skip >= nx) || (skip >= ny))
        {
            printf("flct: skip = %d is too big compared to nx or ny, fatal\n", skip);
            return 1;
        }
        printf("bz0[0]=%f\n", f1[0]);
        printf("bz0[33]=%f\n", f1[33]);
    }
    
    // Figure out size of sliding box in which FFTs will be done
    
    i4 nt, ndmin, nsize = 0;
    double tol = 1e-2;
    if (!sigmaeq0) {
        nt = (i4)sigma * sqrt(log (1. / tol));
        ndmin = (nx < ny) ? (((nx / 3) / 2)) * 2 : (((ny / 3) / 2)) * 2;
        nsize = ((nt / 2) * 2 < ndmin) ? (nt / 2) * 2 : ndmin;
        if (verbose) printf ("flct: nominal sliding box size = %d\n",
                             2 * nsize);
        if(nsize <= 0) {
            printf("flct: error - illegal box size, exiting\n");
            return 1;
        }
    } else {
        nx = 1; ny = 1;     // sigma = 0 means we'll only compute one point
    }

    /* figure out if threshold is in absolute or fractional units
     * and if fractional, convert to absolute units.  If thresh is between
     * zero and 1 (not inclusive) then it's assumed to be fractional,
     * (unless the threshold string ends with 'a') and must be scaled.
     * if aloc == NULL, there's no 'a' in the threshold string.
     */
    
    if (!sigmaeq0) {
        if ((thresh > 0.) && (thresh < 1.) && (aloc == NULL)) {
            double *f1temp = (double *) malloc (sizeof (double) * nx * ny);
            double *f2temp = (double *) malloc (sizeof (double) * nx * ny);
            for (int i = 0; i < nx * ny; i++) {
                /* compute abs value of f1,f2 arrays as f1temp,
                 * f2temp arrays */
                *(f1temp + i) = (double) fabs (*(f1 + i));
                *(f2temp + i) = (double) fabs (*(f2 + i));
            }
            /* now find maximum absolute value of both images */
            i4 iloc1 = maxloc (f1temp, nx * ny);
            i4 iloc2 = maxloc (f2temp, nx * ny);
            double f1max = *(f1temp + iloc1);
            double f2max = *(f2temp + iloc2);
            double fmax = (f1max > f2max) ? f1max : f2max;
            /* now convert relative thresh to absolute threshhold */
            thresh *= fmax;
            if (verbose) 
                printf ("flct: relative threshold in abs. units = %g\n", thresh);
            free (f1temp); free (f2temp);
        }
    }
    if (aloc) free(aloc);
    
    /* Create velocity arrays vx,vy and the velocity mask array vm */
    /* the vm array (velocity mask array, not to be confused with the
     * gaussian mask array that is used to modulate the images)
     * will later be set to 1.0 where the velocity is computed,
     * and to 0.0 where the velocity is not computed -- because the image
     * pixels are below the noise threshold value "thresh" input from
     * the command line. */
    // NOTE: vx, vy, vm already created and passed in

    /* Now create master gaussian image mask: */
    
    double *gaussdata = (double *) malloc (sizeof (double) * (2 * nxorig) * (2 * nyorig));
    if (!sigmaeq0) {    /* this case for sigma > 0 */
        for (int i = 0; i < 2 * nxorig; i++) {
            double argx = sigminv * (double) (i - nxorig);
            for (int j = 0; j < 2 * nyorig; j++) {
                double argy = sigminv * (double) (j - nyorig);
                *(gaussdata + i * (2 * ny) + j) = exp (-argx * argx - argy * argy);
            }
        }
    } else {    /* this case for sigma = 0. ie set gaussian to 1.0 */
        for (int i = 0; i < 2 * nxorig; i++) {
            for (int j = 0; j < 2 * nyorig; j++) {
                *(gaussdata + i * (2 * nyorig) + j) = (double) 1.;
            }
        }
    }

    /* Now do the master loop over i,j for computing velocity field: */
    
    i4 init, icc;
    i4 hardworkneeded, belowthresh;
    i4 imin0, jmax0, jmin0, imax0, imin, jmin, imax, jmax, isize, jsize;
    double fabsbar, f1bar=0., f2bar=0., vxx=0., vyy=0., vmask;
    double shiftx, shifty;
    double *g1, *g2, *absccor;
    
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            
            if ((nx == 1) && (ny == 1)) {
                init = 2; /* special case: must initialize AND destroy plans */
            } else if ((i == 0) && (j == 0) && ((nx + ny) > 2)) {
                /* 1st time through, set init to 1 so that
                 * fftw FFT plans are initialized */
                init = 1;
            }
            else if ((i == (nx - 1)) && (j == (ny - 1)) && ((nx+ny) > 2)) {
                /* last time through, set init to -1 so that
                 * fftw static variables are freed and plans
                 * destroyed */
                init = -1;
            } else {
                /* the rest of the time just chunk along */
                init = 0;
            }
            
            /* Now, figure out if image value is below
             * threshold: */
            /* the data is considered below theshhold if the
             * absolute value of average of the pixel value from the 2 images
             * is below threshold */
            
            fabsbar = 0.5 * (fabs (*(f1 + i * ny + j) + *(f2 + i * ny + j)));
            belowthresh = (fabsbar < thresh);
            
            /* Or alternatively:
             belowthresh = ((fabs1 < thresh) || (fabs2 < thresh));
             */
            /* all the hard work of doing the cross-correlation
             * needs to be done if the avg data is above the
             * threshold OR if init != 0 */
            
            /* added skip logic here */
            // No use
            
            if (skipon) {
                if (transp) {
                    xoffset=qoffset;
                    yoffset=poffset;
                } else {
                    xoffset=poffset;
                    yoffset=qoffset;
                }
                noskipx = !((i - xoffset) % skip);
                noskipy = !((j - yoffset) % skip);
                noskipxy = noskipx * noskipy;
            } else {
                noskipxy=1;
            }
            
            hardworkneeded = (((!belowthresh) && (noskipxy)) || (init != 0));
            
            if (hardworkneeded) {
                
                /* the hard work for this particular pixel starts
                 * now */
                
                /* Now find where the gaussian modulated image
                 * is
                 * chopped off by the sliding box.  The sliding
                 * box is centered at i,j, unless the edges of
                 * the box would go outside the array --
                 * then the
                 * sliding box just sits at edges and/or corners
                 * of the array   */
                
                if (!sigmaeq0) { /* for sigma > 0 */
                    
                    imin0 = (0 > (i - (nsize - 1))) ? 0 : i - (nsize - 1);
                    imax0 = ((nx - 1) < (i + nsize)) ? nx - 1 : i + nsize;
                    imin = (imax0 == nx - 1) ? nx - 1 - (2 * nsize - 1) : imin0;
                    imax = (imin0 == 0) ? 0 + (2 * nsize - 1) : imax0;
                    
                    jmin0 = (0 > (j - (nsize - 1))) ? 0 : j - (nsize - 1);
                    jmax0 = ((ny - 1) < (j + nsize)) ? ny - 1 : j + nsize;
                    jmin = (jmax0 == ny - 1) ? ny - 1 - (2 * nsize - 1) : jmin0;
                    jmax = (jmin0 == 0) ? 0 + (2 * nsize - 1) : jmax0;
                    
                    isize = imax - imin + 1;
                    jsize = jmax - jmin + 1;
                    
                    /* If the following tests for isize, jsize fail,
                     this is very bad:  exit */
                    
                    if (isize != 2 * nsize) {
                        printf ("flct: exiting, bad isize = %d\n", isize);
                        return 1;
                    }
                    if (jsize != 2 * nsize)  {
                        printf ("flct: exiting, bad jsize = %d\n", jsize);
                        return 1;
                    }
                    
                } else {/* if sigma = 0. just set isize=nxorig, jsize=nyorig */
                    
                    isize = nxorig;
                    jsize = nyorig;
                    imin = 0;
                    jmin = 0;
                    
                }
                
                /* Compute sub-image means of f1 and f2: */
                
                f1bar = 0.;
                f2bar = 0.;
                for (int ii = 0; ii < isize; ii++) {
                    for (int jj = 0; jj < jsize; jj++) {
                        f1bar=f1bar + *(f1 + (ii+imin)*nyorig + (jj+jmin));
                        f2bar=f2bar + *(f2 + (ii+imin)*nyorig + (jj+jmin));
                    }
                }
                
                f1bar = f1bar / ((double)isize*jsize);
                f2bar = f2bar / ((double)isize*jsize);
                
                g1 = (double *) malloc (sizeof (double) * isize * jsize);
                g2 = (double *) malloc (sizeof (double) * isize * jsize);
                
                /* Now fill the reduced size arrays (sub-images) with the
                 * appropriate values from the
                 * full-sized arrays: */
                
                for (int ii = 0; ii < isize; ii++) {
                    for (int jj = 0; jj < jsize; jj++) {
                        *(g1 + ii * jsize + jj) =
                            *(gaussdata + (nxorig-i+(ii+imin))*2*nyorig
                              + nyorig-j+(jj+jmin)) *
                            (*(f1 + (ii + imin) * nyorig + (jj + jmin))-f1bar) ;
                        *(g2 + ii * jsize + jj) =
                            *(gaussdata + (nxorig-i+(ii+imin))*2*nyorig
                              + nyorig-j+(jj+jmin)) *
                            (*(f2 + (ii + imin) * nyorig + (jj + jmin))-f2bar) ;
                    }
                }
                
                /* Call to cross_cor is used to find the
                 * relative
                 * shift of image g2 relative to g1: */
                
                icc = cross_cor (init, hires, expand, g1, g2, &absccor,
                                 isize, jsize, &shiftx, &shifty, filter, kr);

                /* Now free up all the temporary arrays created
                 * during the loop */
                
                free (g1);
                free (g2);
                free (absccor);
                
                /* Now convert shifts to units of velocity using
                 * deltas and deltat */
                
                /* Note: if (transp), then the meaning of 
                 * velocities
                 * has to be switched between x and y */
                
                if (transp) {
                    vxx = shifty * deltinv * deltas;
                    vyy = shiftx * deltinv * deltas;
                } else {
                    vxx = shiftx * deltinv * deltas;
                    vyy = shifty * deltinv * deltas;
                }
                
                /* all the hard work for this pixel is now done */
                
            } // hardworkneeded
            
            /* default value for vmask is 1. */
            
            vmask = 1.;
            
            if ((belowthresh || !noskipxy) && !sigmaeq0) {
                
            /* If data below threshold, set vxx, vyy to xmiss 
             * and vmask to 0, meaning vel not computed. */
            /* If sigma=0 just ignore the threshold and compute anyway */

                vxx = xmiss;
                vyy = xmiss;
                vmask = 0.;
            }
            
            if ((j == 0) && (verbose))  {
                printf ("flct: progress  i = %d out of %d\r", i, nx - 1);
                fflush (stdout);
            }
            
            *(vx + i * ny + j) = vxx;
            *(vy + i * ny + j) = vyy;
            *(vm + i * ny + j) = vmask;
            if (verbose && sigmaeq0) {
                printf("\nflct: vx = %g vy = %g \n",vxx,vyy);
                fflush(stdout);
            }
            
        } // j
    } // i
    
    // Skipped all skipon related loops
    // L960-L1124
    
    /* Finally, reset any values of vx,vy that are equal to xmiss to zero, and
     make sure corresponding value of vm is also zero. */
    
    for (int i=0; i < nx; i++)  {
        for(int j=0; j < ny; j++)  {
            if (*(vx+i*ny+j) == xmiss) *(vm+i*ny+j) = 0.;
            if (*(vx+i*ny+j) == xmiss) *(vx+i*ny+j) = 0.;
            if (*(vy+i*ny+j) == xmiss) *(vy+i*ny+j) = 0.;
        }
    }

    /* free the gaussian mask array */
    free (gaussdata);

    if (verbose)
        printf ("\nflct: finished\n");
    
    return 0;

}   // flct4jsoc



// ========================
// Cross correlation?
i4 cross_cor (i4 init, i4 hires, i4 expand, double *arr, double *barr,
              double **absccor, i4 nx, i4 ny, double *shiftx, double *shifty,
              i4 filterflag, double kr)
{
    /*  To use C99 complex arithmetic, uncomment next line */
    /* #include <complex.h> */
    
    /* #include <fftw3.h> */
    i4 i, j, ixx, iyy, maxind, ixmax, iymax, ishft, maxfine, absccmax;
    i4 nxinterp, nyinterp, nfgppergp;
    double normfac, rangex, rangey, shiftx0, shifty0, shiftxx, shiftyy;
    double *xwant, *ywant, *peakarea;
    double shiftsubx, shiftsuby, fx, fy, fxx, fyy, fxy;
    double xmiss=0.;
    
    /* following variables must be saved between calls; declared static */
    
    static double *ina, *inb, *ccor;
    static double *filter, *kx, *ky;
    static fftw_complex *outa, *outb, *ccorconj;
    static fftw_plan pa, pb, pback;
    
    /* absccor is a double pointer containing abs. value of cc function */
    
    /* debug:
     printf("cross_cor: nx = %d, ny = %d\n",nx,ny);
     */
    
    *absccor = malloc (sizeof (double) * nx * ny);
    
    /* Set up interpolation grid depending on whether hires set or not */
    
    if (hires == 1)
    {
        nxinterp = 101;
        nyinterp = 101;
        nfgppergp = 50;
    }
    else
    {
        nxinterp = 21;
        nyinterp = 21;
        nfgppergp = 10;
    }
    /*	printf("initialization stuff done in cross_cor\n"); */
    if ((init == 1) || (init == 2))
    {
        /* First time through: */
        /* Initialization of FFT variables and FFTW plans. */
        /* NOTE -- empirically had to add 1 to "y" dimensions of outa,
         * outb, and ccorconj to
         * avoid a memory leak and a seg fault at fftw_free */
        
        /* should check to see if still a problem */
        
        outa = (fftw_complex *) fftw_malloc (sizeof (fftw_complex) *
                                             nx * ((ny / 2) + 2));
        outb = (fftw_complex *) fftw_malloc (sizeof (fftw_complex) *
                                             nx * ((ny / 2) + 2));	/* valgrind sometimes complains */
        ccorconj = (fftw_complex *) fftw_malloc (sizeof (fftw_complex) *
                                                 nx * ((ny / 2) + 2));
        
        ina = (double *) fftw_malloc (sizeof (double) * nx * ny);
        inb = (double *) fftw_malloc (sizeof (double) * nx * ny);
        ccor = (double *) fftw_malloc (sizeof (double) * nx * ny);
        filter = (double *) fftw_malloc (sizeof (double) * nx * ny);
        kx=(double *) fftw_malloc (sizeof(double)*nx);
        ky=(double *) fftw_malloc (sizeof(double)*ny);
        if(filterflag)
        {
            make_freq(kx,nx);
            make_freq(ky,ny);
            gaussfilt(filter,kx,ky,nx,ny,kr);
        }
        
        for (i = 0; i < nx * ny; i++)
        {
            *(ina + i) = (double) 0.;
            *(inb + i) = (double) 0.;
        }
        for (i = 0; i < nx * ((ny / 2 + 1)); i++)
        {
            *(ccorconj+i)=0.;
            
            /* above line works if complex.h is invoked, but we instead
             * follow the advice on complex arithmetic described in
             * the fftw3 tutorial */
            /*
             ccorconj[i][0] = 0.;
             ccorconj[i][1] = 0.;
             */
        }
        pa = fftw_plan_dft_r2c_2d (nx, ny, ina, outa, FFTW_MEASURE);
        pb = fftw_plan_dft_r2c_2d (nx, ny, inb, outb, FFTW_MEASURE);
        pback = fftw_plan_dft_c2r_2d (nx, ny, ccorconj, ccor, FFTW_MEASURE);
    }
    
    for (i = 0; i < nx * ny; i++)
        
    {
        /*		printf("1st loop: i = %d, *(arr+i)= %g, *(barr+i) = %g\n",
         i,*(arr+i),*(barr+i)); */
        
        /* copy from input doubles to fftw variables */
        
        *(ina + i) = (double) (*(arr + i));
        *(inb + i) = (double) (*(barr + i));
    }
    
    /* actually do the forward FFTs: */
    
    fftw_execute (pa);
    fftw_execute (pb);
    
    /* calculate normalization factor */
    
    normfac = (1. / ((double) nx * ny));
    normfac *= normfac; /* square of above line */
    
    /* Now apply the gaussian filter to the FFT'd data in frequency space */
    
    if(filterflag)
    {
        for (i=0;i<nx;i++)
        {
            for (j=0;j<(ny/2)+1;j++)
            {
                /*
                 outa[i*((ny/2)+1)+j][0]=outa[i*((ny/2)+1)+j][0]*filter[i*ny+j];
                 outa[i*((ny/2)+1)+j][1]=outa[i*((ny/2)+1)+j][1]*filter[i*ny+j];
                 outb[i*((ny/2)+1)+j][0]=outb[i*((ny/2)+1)+j][0]*filter[i*ny+j];
                 outb[i*((ny/2)+1)+j][1]=outb[i*((ny/2)+1)+j][1]*filter[i*ny+j];
                 */
                outa[i*((ny/2)+1)+j]=outa[i*((ny/2)+1)+j]*filter[i*ny+j];
                outb[i*((ny/2)+1)+j]=outb[i*((ny/2)+1)+j]*filter[i*ny+j];
                
            }
        }
    }
    
    /* Now calculate product of conj(outa) * outb */
    
    for (i = 0; i < nx * ((ny/2) + 1); i++)
    {
        
        *(ccorconj+i)=(conj(*(outa+i))*(*(outb+i)))*normfac;
        
        /* the above commented out line can be used if complex.h
         * is invoked, but this is only accepted by C99 compliant
         * compilers.  For now, we used the suggested approach of
         * of the fftw tutorial.  This appears to be required if
         * one uses old e.g. 2.95 or 2.96 versions of gcc or
         * mingw/msys in windows systems */
        /*
         ccorconj[i][0] = (outa[i][0] * outb[i][0] + outa[i][1] * outb[i][1])
         * normfac;
         ccorconj[i][1] = (outa[i][0] * outb[i][1] - outa[i][1] * outb[i][0])
         * normfac;
         */
    }
    
    /* now do the inverse transform to get cc function */
    
    fftw_execute (pback);
    
    /* now calculate the absolute value of cc function */
    
    for (i = 0; i < nx * ny; i++)
    {
        *(*absccor + i) = (double) fabs(*(ccor+i));
    }
    
    if ((init == -1) || (init == 2))
    {
        /* Last time through: free all the plans and static variables */
        
        fftw_free (outa);
        fftw_free (outb);
        fftw_free (ccorconj);
        fftw_free (ccor);
        fftw_free (filter);
        fftw_free (kx);
        fftw_free (ky);
        fftw_free (ina);
        fftw_free (inb);
        fftw_destroy_plan (pa);
        fftw_destroy_plan (pback);
        fftw_destroy_plan (pb);
    }
    
    /* Now shift the absccor array by nx/2, ny/2 to avoid aliasing problems */
    
    ishft = shift2d (*absccor, nx, ny, nx / 2, ny / 2);
    
    /* Now find maximum of the shifted cross-correlation function to 1 pixel
     accuracy:  */
    
    absccmax=1;
    maxind = maxloc (*absccor, nx * ny);
    if( *(*absccor+maxind) == (double)0.)
    {
        absccmax=0;
    }
    if(absccmax == 1)
    {
        ixmax = maxind / ny;
        iymax = maxind % ny;
    }
    else
    {
        ixmax = nx/2;
        iymax = ny/2;
    }
    shiftx0 = ixmax;
    shifty0 = iymax;
    shiftsubx=0.;
    shiftsuby=0.;
    
    if(( expand == 1) && (hires == -1) && (ixmax > 0) && (ixmax < (nx-1))
       && (iymax > 0) && (iymax < (ny-1)) && (absccmax == 1))
    {
        fx=0.5* ( *(*absccor+(ixmax+1)*ny+iymax) -
                 *(*absccor+(ixmax-1)*ny+iymax) );
        fy=0.5* ( *(*absccor+ixmax*ny+iymax+1) - *(*absccor+ixmax*ny+iymax-1) );
        fxx = ( *(*absccor+(ixmax+1)*ny+iymax)+ *(*absccor+(ixmax-1)*ny+iymax)
               -2.*( *(*absccor+ixmax*ny+iymax))  );
        fyy = ( *(*absccor+ixmax*ny+iymax+1) + *(*absccor+ixmax*ny+iymax-1)
               -2.*( *(*absccor+ixmax*ny+iymax)) );
        fxy = 0.25*( *(*absccor+(ixmax+1)*ny+iymax+1) +
                    *(*absccor+(ixmax-1)*ny+iymax-1) -
                    *(*absccor+(ixmax+1)*ny+iymax-1) -
                    *(*absccor+(ixmax-1)*ny+iymax+1) );
        /* In following expressions for subshifts, shift is in units of pixel length */
        shiftsubx=(fyy*fx-fy*fxy)/(fxy*fxy-fxx*fyy);
        shiftsuby=(fxx*fy-fx*fxy)/(fxy*fxy-fxx*fyy);
    }
    
    shiftxx=shiftx0 + shiftsubx;
    shiftyy=shifty0 + shiftsuby;
    /*
     printf("shiftx0-nx/2 = %g\n",(shiftx0-(double)(nx/2)));
     printf("shifty0-ny/2 = %g\n",(shifty0-(double)(ny/2)));
     
     */
    
    /* Now create x, y arrays to define desired interpolation grid: */
    if(hires != -1)
    {
        
        rangex = (double) (nxinterp - 1) / nfgppergp;
        rangey = (double) (nyinterp - 1) / nfgppergp;
        
        xwant = (double *) malloc (sizeof (double) * nxinterp);
        ywant = (double *) malloc (sizeof (double) * nyinterp);
        
        for (i = 0; i < nxinterp; i++)
        {
            *(xwant + i) = ((((double) i) * rangex) / ((double) (nxinterp - 1)))
            - 0.5 * rangex + shiftx0;
            /*                 printf("xwant[%d] = %g\n",i,*(xwant+i)); */
        }
        for (j = 0; j < nyinterp; j++)
        {
            *(ywant + j) = ((((double) j) * rangey) / ((double) (nyinterp - 1)))
            - 0.5 * rangey + shifty0;
            /*                 printf("ywant[%d] = %g\n",j,*(ywant+j)); */
        }
        
        /* Now, do the interpolation of the region of the peak of the cc fn */
        
        interpcc2d (*absccor, xmiss, nx, ny, xwant, nxinterp, ywant,
                    nyinterp, &peakarea);
        
        /* Following writeimage stmt is available if you need to examine the
         * peakarea array for debugging - note transpose of x,y for IDL  read */
        
        /*
         transp=1;
         writeimage("peakarea.dat",peakarea,nxinterp,nyinterp,transp);
         */
        
        /* Now find the peak of the interpolated function */
        
        maxfine = maxloc (peakarea, nxinterp * nyinterp);
        ixx = maxfine / nyinterp;
        iyy = maxfine % nyinterp;
        /* Here is where to compute sub-pixel shifts in peakarea if they're wanted */
        shiftsubx=0.;
        shiftsuby=0.;
        if((expand == 1) && (ixx > 0) && (ixx < (nxinterp-1)) && (iyy > 0)
           && (iyy < (nyinterp-1)))
        {
            fx=0.5* ( *(peakarea+(ixx+1)*nyinterp+iyy) -
                     *(peakarea+(ixx-1)*nyinterp+iyy) );
            fy=0.5* ( *(peakarea+ixx*nyinterp+iyy+1) -
                     *(peakarea+ixx*nyinterp+iyy-1) );
            fxx = ( *(peakarea+(ixx+1)*nyinterp+iyy)+
                   *(peakarea+(ixx-1)*nyinterp+iyy)
                   -2.*( *(peakarea+ixx*nyinterp+iyy))  );
            fyy = ( *(peakarea+ixx*nyinterp+iyy+1) + *(peakarea+ixx*nyinterp+iyy-1)
                   -2.*( *(peakarea+ixx*nyinterp+iyy)) );
            fxy = 0.25*( *(peakarea+(ixx+1)*nyinterp+iyy+1) +
                        *(peakarea+(ixx-1)*nyinterp+iyy-1) -
                        *(peakarea+(ixx+1)*nyinterp+iyy-1) -
                        *(peakarea+(ixx-1)*nyinterp+iyy+1) );
            /* In following expressions for subshifts, must mpy by unit of pixel length */
            shiftsubx=((fyy*fx-fy*fxy)/(fxy*fxy-fxx*fyy))*
            rangex/((double) (nxinterp -1));
            shiftsuby=((fxx*fy-fx*fxy)/(fxy*fxy-fxx*fyy))*
            rangey/((double) (nyinterp -1));
        }
        shiftxx = *(xwant + ixx) + shiftsubx;
        shiftyy = *(ywant + iyy) + shiftsuby;
        /* Free the variables created during interpolation */
        free (xwant);
        free (ywant);
        free (peakarea);
    }
    
    /* Now, assign values to shiftx, shifty to return to calling program */
    
    *shiftx = shiftxx - (double) (nx / 2);
    *shifty = shiftyy - (double) (ny / 2);
    
    /* Following expressions used if only 1 pixel accuracy needed from absccor 
     *
     *shiftx=((double)ixmax)-(double)(nx/2);
     *shifty=((double)iymax)-(double)(ny/2);
     */
    
    return 0;
}


// ===========================
/* Assumes kx of size nx, ky of size ny, and filter of size (nx,ny) */
i4 gaussfilt(double *filter, double *kx, double *ky, i4 nx, i4 ny, double kr)
{
    double kxmax,kymax,kxroll,kyroll,smxinv,smyinv,argx,argy;
    i4 i,j;
    kxmax=(double)kx[nx/2];
    kymax=(double)ky[ny/2];
    kxroll=kr*kxmax;
    kyroll=kr*kymax;
    smxinv=(double)1./kxroll;
    smyinv=(double)1./kyroll;
    for (i=0;i<nx;i++)
    {
        argx=kx[i]*smxinv;
        for(j=0;j<ny;j++)
        {
            argy=ky[j]*smyinv;
            filter[i*ny+j]=exp( -(argx*argx + argy*argy) );
        }
    }
    return 0;
}

// ===========================
/* Circular shift of the x,y indices of array *arr by ishift,jshift */
/* This function is similar to the shift function in IDL.  nx, ny
 * are the assumed dimensions of the array */

i4 shift2d (double *arr, i4 nx, i4 ny, i4 ishift, i4 jshift)
{
    double *temp;
    i4 i, j, ii, jj;
    temp = (double *) malloc (sizeof (double) * nx * ny);
    for (i = 0; i < nx; i++)
    {
        ii = (i + ishift) % nx;	/* ii = (i + ishift) modulo nx */
        
        for (j = 0; j < ny; j++)
        {
            jj = (j + jshift) % ny;	/* jj = (j+jshift) modulo ny */
            
            /* Now members of temp array get shifted: */
            
            *(temp + ii * ny + jj) = *(arr + i * ny + j);
        }
    }
    
    /* Now copy temp array back into arr, then destroy temp and return */
    
    memcpy ((void *) arr, (void *) temp, nx * ny * sizeof (double));
    free (temp);
    return 0;
}

// ===========================
/* finds the location of the maximum of the double array *arr and returns it. */

i4 maxloc (double *arr, i4 xsize)
{
    i4 i, location;
    double amax;
    /* initialize amax and location to 0th element */
    amax = *(arr + 0);
    location = 0;
    for (i = 1; i < xsize; i++)
    {
        if (*(arr + i) > amax)
        {
            amax = *(arr + i);
            location = i;
        }
    }
    return location;
}

// ===========================
/* finds the location of the minimum of the double array *arr and returns it. */
i4 minloc (double *arr, i4 xsize)
{
    i4 i, location;
    double amin;
    /* initialize amin and location to 0th element */
    amin = *(arr + 0);
    location = 0;
    for (i = 1; i < xsize; i++)
    {
        if (*(arr + i) < amin)
        {
            amin = *(arr + i);
            location = i;
        }
    }
    return location;
}

// ===========================
/* k is assumed already allocated in main program, with dimension ndim */

i4 make_freq(double *k, i4 ndim)
{
    /* k is assumed already allocated in main program, with dimension ndim */
    i4 n21,i,inext;
    n21=(ndim/2)-1;
    for (i=0;i<n21+1;i++)
    {
        k[i]=(double)i;
    }
    
    inext=n21+1;
    if((ndim/2)*2 != ndim)
    {
        k[inext]=(double)(ndim/2);
        inext++;
        k[inext]=-(double)(ndim/2);
        inext++;
    }
    
    else
    {
        k[inext]=(double)(ndim/2);
        inext++;
    }
    
    for (i=inext;i<ndim;i++)
    {
        k[i]=-(double)(n21-(i-inext));
    }
    /* debug */
    
    /*
     for (i=0;i<ndim;i++)
     {
     printf("i = %d, k = %g\n",i,k[i]);
     }
     */
    
    /* end debug */
    
    return 0;
}


//================================
i4 interpcc2d (double *fdata, double xmiss, i4 nx, i4 ny,
               double *xwant, i4 nxinterp, double *ywant, i4 nyinterp, double **finterp)
{
    /*
     * This function does cubic convolution interpolation onto an array
     * finterp from data defined in array fdata.  nx, ny are the
     * assumed dimensions of fdata, and nxinterp, nyinterp are the
     * assumed dimensions of finterp.  The values of x,y at which
     * the interpolation is desired are passed in through the arrays
     * xwant and ywant, which are dimensioned nxinterp and nyinterp,
     * respectively.  It is assumed that xwant, ywant are in units of
     * the indices of the original data array (fdata),
     * treated as floating point (double precision, actually)
     * numbers. Arrays fdata, xwant, and ywant are passed in
     * as pointers; The array finterp is defined in this function
     * as a "double" pointer and the array is created and passed back to
     * the calling function.  In the calling function, finterp is declared
     * as a pointer, but when it is passed into this function as
     * an argument, the address of the pointer is used.
     *
     * if any of the datapoints within a kernel weighting distance of
     * xwant and ywant are equal to xmiss,
     * the returned value of finterp is also set to xmiss.  xmiss is a user-
     * defineable calling argument.
     */
    
    double *cdata;
    /*  double txt, tyt, xint, yint, ftmp, xmiss = 0.; */
    double txt, tyt, xint, yint, ftmp;
    
    /* Logic for a user-defined value of xmiss has been added.  Previously
     * was just set to 0 as a local variable */
    
    double tx, ty, rx, ry;
    i4 i, ii, j, jj, itemp, jtemp, izero, jzero, databad;
    /*  i4 transp; */
    
    /* Now, create the cdata array, bigger by 1 gp than fdata
     * all around the borders: */
    
    cdata = (double *) malloc (sizeof (double) * (nx + 2) * (ny + 2));
    
    /* Now fill the interior of cdata with fdata */
    
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            *(cdata + (i + 1)*(ny + 2) + (j + 1)) = *(fdata + i*ny + j);
        }
    }
    
    /*
     * The basic concept for filling in edges and corners of cdata is this:
     * The edge point is equal to 3*(value of adjacent point)
     * -3*value(next to adjacent point) + 1*value(3rd point).  This
     * prescription yields an extrapolation which is consistent with
     * a 3rd (or is it 4th?) order Taylor expansion of the function
     * evaluated at the last real gridpoint, and extrapolated to the
     * edge point.  This procedure is followed
     * thoughout here, though I think it isn't really correct for the
     * corner points because there I think an expansion from both
     * both directions should be done.  But no harm seems to be done
     * to the results.
     */
    
    /* Fill in the edges of cdata: */
    
    for (j = 0; j < ny; j++)
    {
        
        /* left and right edges: */
        
        *(cdata + 0*(ny + 2) + (j+1)) = *(fdata + 2*ny + j)
        - 3. * (*(fdata + 1*ny + j)) + 3. * (*(fdata + 0*ny + j));
        
        *(cdata + (nx + 1)*(ny + 2) + (j + 1)) = *(fdata + (nx - 3)*ny + j)
        - 3. * (*(fdata + (nx - 2)*ny + j)) + 3. * (*(fdata + (nx - 1)*ny + j));
    }
    for (i = 0; i < nx; i++)
    {
        
        /* bottom and top edges: */
        
        *(cdata + (i + 1)*(ny + 2) + 0) = *(fdata + i*ny + 2)
        - 3. * (*(fdata + i*ny + 1)) + 3. * (*(fdata + i*ny + 0));
        
        *(cdata + (i + 1)*(ny + 2) + ny + 1) = *(fdata + i*ny + ny - 3)
        - 3. * (*(fdata + i*ny + ny - 2)) + 3. * (*(fdata + i*ny + ny - 1));
    }
    
    /* Now fill in the 4 corners: */
    
    *(cdata + 0*(nx + 2) + 0) =
    3. * (*(cdata + 1*(ny + 2) + 0)) -
    3. * (*(cdata + 2*(ny + 2) + 0)) + *(cdata + 3*(ny + 2) + 0);
    
    *(cdata + (nx + 1)*(ny + 2) + 0) =
    3. * (*(cdata + nx*(ny + 2) + 0)) -
    3. * (*(cdata + (nx - 1)*(ny + 2) + 0)) + *(cdata +
                                                (nx - 2)*(ny + 2) + 0);
    
    *(cdata + 0*(ny + 2) + ny + 1) =
    3. * (*(cdata + 0*(ny + 2) + ny)) -
    3. * (*(cdata + 0*(ny + 2) + ny - 1)) + *(cdata + 0*(ny + 2) + ny - 2);
    
    *(cdata + (nx + 1)*(ny + 2) + ny + 1) =
    3. * (*(cdata + nx*(ny + 2) + ny + 1)) -
    3. * (*(cdata + (nx - 1)*(ny + 2) + ny + 1)) + *(cdata +
                                                     (nx - 2)*(ny + 2) +
                                                     ny + 1);
    
    /* Now create the space for finterp */
    
    *finterp = (double *) malloc (sizeof (double) * nxinterp * nyinterp);
    
    /* Now interpolate onto the desired grid */
    
    for (i = 0; i < nxinterp; i++)
    {
        /* starting the outer loop over x */
        
        xint = *(xwant + i);
        
        /* make sure izero is in bounds */
        
        itemp = ((i4) xint > 0) ? (i4) xint : 0;
        izero = (itemp < (nx - 2)) ? itemp : nx - 2;
        for (j = 0; j < nyinterp; j++)
        {
            /* starting the outer loop over y */
            
            yint = *(ywant + j);
            if ((yint < 0.) || (yint > (double) (ny - 1))
                || ((xint < 0) || (xint > (double) (nx - 1))))
            {
                /* if data missing, set interp to xmiss */
                
                /* Debug
                 printf("interpccd2: i=%d,j=%d gets finterp[i,j] set to xmiss\n",
                 i,j);
                 */
                *(*finterp + i * nyinterp + j) = xmiss;
            }
            else
            {
                /* make sure jzero is in bounds */
                
                jtemp = ((i4) yint > 0) ? (i4) yint : 0;
                jzero = (jtemp < (ny - 2)) ? jtemp : ny - 2;
                
                /* initialize the temporary finterp value */
                
                ftmp = (double) 0.;
                
                /* start the innermost loops over neighboring
                 * data points*/
                
                databad=0;
                for (ii = -1; ii < 3; ii++)
                {
                    txt = xint - (double) (izero + ii);
                    tx = (double) fabs (txt);
                    
                    /* evaluate kernel wt function r(tx): */
                    
                    /* Note no testing for out of range
                     * values of |tx| or |ty| > 2 --
                     * we assume the tx, ty are properly
                     * computed such that their absolute
                     * value never exceeds 2. */
                    
                    rx = (tx >= (double) 1.0) ?
                    (((((double) (-0.5)) * tx +
                       ((double) 2.5)) * tx) -
                     (double) 4.) * tx + (double) 2. :
                    (((double) (1.5)) * tx -
                     ((double) (2.5))) * tx * tx + (double) 1.;
                    
                    for (jj = -1; jj < 3; jj++)
                    {
                        
                        tyt = yint - (double) (jzero + jj);
                        ty = (double) fabs (tyt);
                        
                        /* evaluate kernel weighting
                         * function r(ty): */
                        
                        ry = (ty >= (double) 1.0) ?
                        (((((double) (-0.5)) * ty +
                           ((double) 2.5)) * ty) -
                         (double) 4.) * ty + (double) 2. :
                        (((double) (1.5)) * ty -
                         ((double) (2.5))) * ty * ty + (double) 1.;
                        
                        /* do the cubic convolution
                         * over the neighboring data
                         * points, using the x and
                         * y evaluated kernel weighting
                         * functions rx and ry: */
                        
                        ftmp = ftmp +
                        *(cdata + (izero + 1 + ii)*(ny + 2)
                          + jzero + 1 + jj) * rx*ry;
                        if( *(cdata+(izero+1+ii)*(ny+2)+jzero+1+jj) == xmiss)
                            databad=1;
                    }
                }
                /* now assign this value to interpolated
                 array, unless one of the values was xmiss: */
                if(databad)
                {
                    /* Debug
                     printf("interpcc2d: i=%d,j=%d gives databad\n",i,j);
                     */
                    *(*finterp + i*nyinterp + j) = xmiss;
                }
                else
                {
                    *(*finterp + i*nyinterp + j) = ftmp;
                }
            }
        }
    }
    
    
    /* DEBUG
     transp=1;
     writeimage("cdata.dat",cdata,nx+2,ny+2,transp);
     */
    
    /* free the cdata space */
    free (cdata);
    
    /* we're done */
    
    return 0;
}


// ===========================
/* This function returns 1 if it is a big endian machine, 0 otherwise */

i4 is_big_endian ()
{
    const unsigned char fakeword[4] = { 0xFF, 0x00, 0xFF, 0x00 };
    i4 realword = *((i4 *) fakeword);
    if (realword == 0xFF00FF00)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}
