#include "expmax.h"

double exp_plus_quadratic_function (double x, void *params)
     {
       struct exp_plus_quadratic_params *p 
         = (struct exp_plus_quadratic_params *) params;
     
       double A0 = p->A0;
       double A1 = p->A1;
       double A2 = p->A2;
       double A3 = p->A3;
       double A4 = p->A4;
       double A5 = p->A5;
     
       double z = (x - A1) / A2;
       return -(A0 * exp (-z * z / 2.) +A3 + A4 *x + A5 * x * x);
     }

