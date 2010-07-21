#include "expfit.h"

/*--------------------------------------------------------------------------*/
/* expfit.c -- model functions for exponential + quadratic */
/* Marcelo Emilio (c) Jan 2009 Modified from GNU Library manual */
     
int expb_f (const gsl_vector * x, void *data, 
             gsl_vector * f)
     {
       size_t n = ((struct dataF *)data)->n;
       double *y = ((struct dataF *)data)->y;
       double *t = ((struct dataF *)data)->t;
       double *sigma = ((struct dataF *) data)->sigma;
     
       double A0 = gsl_vector_get (x, 0);
       double A1 = gsl_vector_get (x, 1);
       double A2 = gsl_vector_get (x, 2);
       double A3 = gsl_vector_get (x, 3);
       double A4 = gsl_vector_get (x, 4);
       double A5 = gsl_vector_get (x, 5);
     
       size_t i;
     
       for (i = 0; i < n; i++)
         {
           /* Model Yi = A0 * exp(-z^2/2) + A3 + A4*i +A5*i*i  */
           /* z=(i-A1)/A2                                     */
           /*double t = i;*/
           double z = (t[i]-A1)/A2;
           double Yi = A0 * exp (-z*z/2.) + A3 +A4*t[i] + A5*t[i]*t[i];
           gsl_vector_set (f, i, (Yi - y[i])/sigma[i]);
         }
     
       return GSL_SUCCESS;
     }

/*--------------------------------------------------------------------------*/

int expb_df (const gsl_vector * x, void *data, 
              gsl_matrix * J)
     {
       size_t n = ((struct dataF *)data)->n;
       double *sigma = ((struct dataF *) data)->sigma;
       double *t = ((struct dataF *)data)->t;

       double A0 = gsl_vector_get (x, 0);
       double A1 = gsl_vector_get (x, 1);
       double A2 = gsl_vector_get (x, 2);
       //double A3 = gsl_vector_get (x, 3);
       //double A4 = gsl_vector_get (x, 4);
       //double A5 = gsl_vector_get (x, 5);

       size_t i;
     
       for (i = 0; i < n; i++)
         {
           /* Jacobian matrix J(i,j) = dfi / dxj, */
           /* where fi = (Yi - yi)/sigma[i],      */
           /* Model Yi = A0 * exp(-z^2/2) + A3 + A4*i +A5*i*i  */
           /* z=(i-A1)/A2                         */
           /* and the xj are the parameters (A0,A1,A2,A3,A4,A5) */
           /*double t = i;*/
           double s = sigma[i];
           double z = (t[i]-A1)/A2;
           double e = exp(-z*z/2.);
           gsl_matrix_set (J, i, 0, e/s); 
           gsl_matrix_set (J, i, 1, A0/A2 * z * e/s);
           gsl_matrix_set (J, i, 2, A0/A2 * z * z * e/s);
           gsl_matrix_set (J, i, 3, 1/s);
           gsl_matrix_set (J, i, 4, t[i]/s);
           gsl_matrix_set (J, i, 5, t[i] * t[i]/s);

         }
       return GSL_SUCCESS;
     }
     
int expb_fdf (const gsl_vector * x, void *data,
               gsl_vector * f, gsl_matrix * J)
     {
       expb_f (x, data, f);
       expb_df (x, data, J);
     
       return GSL_SUCCESS;
     }


