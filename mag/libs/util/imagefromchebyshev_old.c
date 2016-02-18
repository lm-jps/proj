/*
 * Name:		imagefromchebyshev.c
 *
 * Description:		Recover 2D image from the Chebyshev coefs
 *
 * Function List:
 *			cheby_basis(int n, int order, double x[], cheb[], int ldc)
 *			fitimage(double *image, int m, int n, int order, int fill)
 *
 * Called by:		jpolfil.c
 *
 * Original source:	Rasmus Larsen (rmlarsen@gmail.com)
 * Adapted by:		Xudong Sun (xudongs@stanford.edu)
 *
 * Version:		v1.0		Sep 25 2009
 *
 * Issues:
 *			
 */

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif


/* Prototypes for Fortran BLAS and LAPACK routines */
void dgels_(const char *, int *, int *, int *, double *, int *, double *, int *, double *, int *, int *);



/* ################ Chebyshev Polynomials ################ */


/* Generate matrix with Chebychev polynomials T_0,T_1,...,T_{order-1}
 * evaluated on the gridpoints in x.
 */
static void cheby_basis(n, order, x, cheb, ldc)
    int n, order;
    double x[], cheb[];
    int ldc;
{
    int i, s;

    for (i = 0; i < n; i++) {
        cheb[i] = 1;
        cheb[i + ldc] = x[i];
    }  
    for (s = 2; s < order; s++) {
        for (i = 0; i < n; i++) {
            cheb[s * ldc + i] = 
                2 * x[i] * cheb[(s - 1) * ldc + i] - cheb[(s - 2) * ldc + i];
        }
    }
}  



/* ################ Image Fitting ################ */


void imagefromchebyshev(image, m, n, order, coef, xdelta, ydelta)
    double *image, xdelta, ydelta;
    int m, n, order;
    double coef[];
{
    int i, j, k, l, s, t, npix, lda, mn;
    int nterms, info, lwork, *map[2], nimgs;
    double *A, *B, *x, *y, *work, *chebx, *cheby;
//    double *A, *x, *y, *work;
    double pix;

    if (m <= 0 || n <= 0 || order <= 1)
        return;
    nimgs = 1;

    /* Allocate memory for matrix, right-hand sides etc. */
    lda = m * n;
//    A = (double *)malloc(m * n * order * order * sizeof(double));
//    B = (double *)malloc(m * n * sizeof(double));
//    image = (float *)malloc(m * n * sizeof(float));
//    work = (double *)malloc(m * n * sizeof(double));
    map[0] = (int *)malloc(n * m * sizeof(int));
    map[1] = (int *)malloc(n * m * sizeof(int));
    x = (double *)malloc(n * sizeof(double));
    y = (double *)malloc(m * sizeof(double));
    chebx = (double *)malloc(order * n * sizeof(double));
    cheby = (double *)malloc(order * m * sizeof(double));
    nterms = order * order;

    /* Set up abscissa grid on [-1;1]x[-1;1] and Chebychev polynomials. */
    for (i = 0; i < n; i++) 
        x[i] = (double)(2.0 * (i + xdelta)) / (n - 1) - 1;
    cheby_basis(n, order, x, chebx, n);
    for (j = 0; j < m; j++)
        y[j] = (double)(2.0 * (j + ydelta)) / (m - 1) - 1;
    cheby_basis(m, order, y, cheby, m);

    /* Extract "good" pixels and set up matrix with bivariate polynomials
       in the corresponding points. */
/*
    npix = 0;
    for (i = 0; i < n; i++)
    for (j = 0; j < m; j++) {
        pix = image[i * m + j];
        if (!isnan(pix)) {
            B[npix] = pix;
            map[0][npix] = i;
            map[1][npix] = j;
            npix++;
        }
    }
    for (s = 0; s < order; s++)
    for (t = 0; t < order; t++)
    for (l = 0; l < npix; l++) {
        i = map[0][l];
        j = map[1][l];
        A[(m * n) * (s * order + t) + l] = chebx[n * s + i] * cheby[m * t + j];
    }
*/

    /* Extract remaining right-hand sides */
/*
    for (i = 0; i < npix; i++) {
        j = map[0][i] * m + map[1][i];
        B[i] = image[j];
    }
*/
    /* Solve the least-squares problems. */
/*
    mn = MIN(npix, nterms);
    lwork = MAX(1, mn + MAX(npix, MAX(nterms, 1)) * 32);
    work = (double *)malloc(lwork * sizeof(double));

    dgels_("n", &npix, &nterms, &nimgs, A, &lda, B, &lda, work, &lwork, &info);

    if (info < 0) {
        fprintf(stderr, "DGELS: The %dth argument had an illegal value.\n", -info);
        exit(1);
    }
*/
    /* Overwrite image with the smoothed images corresponding to the
       least-squares solutions. */

//    if (fill) {
        for (i = 0; i < n; i++)
	for (j = 0; j < m; j++) {
	    image[i * m + j] = 0;
        }
        for (s = 0; s < order; s++)
	for (t = 0; t < order; t++)
	for (i = 0; i < n; i++)
	for (j = 0; j < m; j++) {
	    image[i * m + j] +=
                coef[s * order + t] * chebx[n * s + i] * cheby[m * t + j];
//printf("%lf,%lf,%lf\n", B[s * order + t], chebx[n * s + i], cheby[m * t + j]);
        }

//printf("%lf,%lf,%lf\n", coef[12 * 15 + 9], chebx[2048 * 10 + 2048], cheby[2048 * 10 + 2048]);

/*
    } else {
        for (l = 0; l < npix; l++) {
	    i = map[0][l];
	    j = map[1][l];
	    image[i * m + j] = 0;
        }
        for (s = 0; s < order; s++) 
	for (t = 0; t < order; t++)
	for (l = 0; l < npix; l++) {
	    i = map[0][l];
	    j = map[1][l];
	    image[i * m + j] +=
                B[s * order + t] * chebx[n * s + i] * cheby[m * t + j];
	}
    }
*/

//    for (i = 0; i < order * order; i++) {
//        coef[i] = B[i];
//        printf("%lf\n", B[i]);
//        }

    /* Free memory */
    free(map[0]); free(map[1]);
    free(x); free(y);
    free(chebx); free(cheby);
//    free(A); free(B); free(work);
//    return A;
}


#undef MAX
#undef MIN
