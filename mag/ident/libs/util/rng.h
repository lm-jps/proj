/*
 * Include file for
 *
 * Random number generators implemented by Mike Turmon between 1995-1997.
 * 
 * See also rng.c
 *  
 */

#ifndef _rng_h_
#define _rng_h_

#ifdef __cplusplus
extern "C" {
#ifdef NOT_DEFINED
} /* fool emacs */
#endif
#endif

/* status utilities */
char *rng_name(void);             /* plain text name of generator */
int   rng_numeric_name(void);     /* code (small int) for generator */

/* initialization */
int   rng_init_state_seeded(int); /* set state with supplied seed */
                                  /* returns supplied seed */
int   rng_init_state(void);       /* as above, will invent own seed */
                                  /* returns invented seed */
int   rng_new_seed(void);         /* invent a seed from clock and pid */

/* get and set state information */
int   rng_state_length(void);     /* sizeof state info in bytes */
int   rng_get_state(char *, int); /* put <= arg2 bytes of state in arg1 */
                                  /* returns #bytes placed, 
				     or 0 if too much state for arg1 */
int   rng_set_state(char *, int); /* set <= arg2 bytes of state from arg1 */
                                  /* returns #bytes placed, 
				     or 0 if too much state for arg1 */
// real random numbers
double rng_uniform(void);          /* return U([0,1]) */

// normal
double rng_normal(void);           /* return N(0,1) */
double rng_normal_old(void);       /* return N(0,1) */
double rng_normal_sdev(double mu, double sigma); /* N(mu,sigma^2) */
double rng_normal_var(double mu, double sigma2); /* N(mu,sigma2) */

// continuous, non-normal
double rng_gamma(double, double); /* return gamma random number */
double rng_beta(double, double);  /* return beta random number */
double rng_chi2(double r);        /* return chi^2 with r DOF */
double rng_exponential(void);     /* return exponential(1) */

// discrete
double rng_poisson(double);       /* return poisson */
double rng_geometric(double);     /* return geometric */
double rng_binomial(double, double); /* return binomial */

// vector/matrix
int rng_normal_vector(double* x, double* mu, double* sigma, int n, int d, 
		      int stridex1, int stridex2, int stridemu);
int rng_wishart(double* sigma, double dof, double scale, int N, double* w);

#ifdef __cplusplus
#ifdef NOT_DEFINED
{ /* fool emacs */
#endif
}	/* extern "C" */
#endif

#endif /* _rng_h_ */
