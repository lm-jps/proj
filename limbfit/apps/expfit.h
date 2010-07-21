/*--------------------------------------------------------------------------*/
/* expfit.c -- model functions for exponential + quadratic */
/* Marcelo Emilio (c) Jan 2009 Modified from GNU Library manual */

#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_rng.h>                                      
#include <gsl/gsl_randist.h>                                  
#include <gsl/gsl_vector.h>                                   
#include <gsl/gsl_blas.h>                                     
#include <gsl/gsl_multifit_nlin.h>                            
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_min.h>    

int expb_f (const gsl_vector * x, void *data, gsl_vector * f);
int expb_df (const gsl_vector * x, void *data, gsl_matrix * J);
int expb_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J);

     
struct dataF {
       size_t n;
       double *t;
       double * y;
       double * sigma;
     };
     

