/*
* Random number generators implemented by Mike Turmon between 1995-1997.
* 
* Normal RNG improved for better speed, 2003.
* 
* Some are LCG's, some are GFSR's.  Also includes a recent 'twisted'
* GFSR, which is not mine.  The TGFSR may not be 64-bit ready.
* 
* Notes: 
*  -- random is the standard unix shift register algorithm, with 256
*  bytes of state (64 int's).
*  -- drand48 is a linear congruential generator (LCG) with multiplier
*  a = 25214903917, offset b = 11, and modulus M = 2^48.
*  -- both minimal generators are LCG's with a = 16807, b = 0, M = 2^31-1.
*  See S.K.Park and K.W.Miller, "Random number generators: good ones are hard
*  to find", Comm. ACM, October 1988, vol 31, number 10, pp 1192-1201, and 
*  many other refs.
*  -- GFSR521 is the 521-bit shift-register generator described in:
*  B.D. Ripley, "Thoughts on pseudorandom number generators," 
*  J. Comput. Appl. Math., 31:153-163, 1990. 
*  -- TGFSR is the very recent `twisted' shift-register of Matsumoto 
*  and Nishimura
* Not here:
*  -- EICG, the particular modulus p=2^31-1 inversive generator described
*  by P. Hellekalek in "Inversive Pseudorandom Number Generators:
*  Concepts, Results, and Links,"  In: C. Alexopoulos et al.,
*  editors, Proceedings of the 1995 Winter Simulation Conference, 
*  pages 255-262, 1995.  The EICG iteration is defined by:
*              y_n = inv(a*(n_0 + n) + b) (mod p)      n >= 0
*  where y=inv(z) means y*z = 1 (mod p), which is unique for prime p.
*/

/* Timings on a sun ultrasparc170E, with the sparcworks v4 compiler and
* option -fast, for 10^6 random numbers, are given below:
*  <Generator>                                 <time>    <time/0.26>
*  random from UNIX:				0.70 sec  (2.69)
*  drand48 from UNIX:				2.26 sec  (8.69)
*  LCG of Park/Miller (int version):		0.86 sec  (3.31)
*  LCG of Park/Miller (real version):		2.17 sec  (8.35)
*  GFSR of length-521 of Ripley 1990:		0.26 sec  (1.00)
*  Matsumoto `twisted' GFSR:			0.36 sec  (1.38)
*/

/*LINTLIBRARY*/

#include <sys/types.h>  /* for getpid() */
#include <sys/time.h>   /* for gettimeofday() */
#include <stdio.h>      /* for fprintf() */
#include <stdlib.h>     /* for calloc() */
#include <unistd.h>     /* for getpid() */
#include <string.h>     /* for memcpy() */
#include <math.h>       /* for M_PI, sin, log */
#include <time.h>       /* clock(), time() */
#include "rng.h"        /* for consistency */

#define RANDOM       1 /* random()  random number generator, a SR */
#define DRAND48      2 /* drand48() random number generator, an LCG */
#define MINIMAL_INT  3 /* "minimal" LCG of Park and Miller, integer version */
#define MINIMAL_REAL 4 /* "minimal" LCG of Park and Miller, real version */
#define GFSR521      5 /* 521-bit SR of Ripley 1990 */
#define TGFSR        6 /* `twisted' SR of Matsumoto and Nishimura 1997 */

#ifndef GENERATOR
#define GENERATOR    GFSR521
#endif

/* == GLOBAL DEFINITIONS ============================================== */

#define GenIs(g) (GENERATOR == g)  /* are we using generator g? */

#include <strings.h>
#if GenIs(MINIMAL_REAL)
#include <math.h> 
#elif GenIs(RANDOM)
long  random(void);
int   srandom(unsigned);
char *initstate(unsigned, char *, int);
char *setstate(char *);
#elif GenIs(DRAND48)
#include <stdlib.h>
#endif

/* == UTILITIES ======================================================= */

static int 
rng_is_bigendian(void)
{
  union {
    long l;
    char c[sizeof(long)];
  } u;

  u.l = 1L;
  return(u.c[sizeof(long) - 1] == 1);
}


/* 
 * generate a double nan.  routine is local to this file
 */
static double 
rng_getnand(void)
{
  /* set the inf/nan buffer up only once */
  static int ieee_nan_inited = 0;
  /* the union ensures the underlying double is aligned properly */
  static union {
    double d;
    unsigned char c[8];
  } nan;

  if (!ieee_nan_inited) {
    int i;
    for (i = 0; i < 8; i++)
      nan.c[i] = 1;  /* have seen 0 here also */
    if (rng_is_bigendian()) {
      nan.c[0] = 0x7F;
      nan.c[1] = 0xF0; /* with 0 above goes 0xF8 here */
    } else {
      nan.c[7] = 0x7F;
      nan.c[6] = 0xF0; /* ditto */
    }
    ieee_nan_inited = 1; /* flag: buffer is set up now */
  }
  return(nan.d);
}



/*
 * short text description of the current generator
 */
char *
rng_name(void)
{
#if GenIs(RANDOM)
return("random from UNIX");

#elif GenIs(MINIMAL_REAL)
return("LCG of Park/Miller (real version)");

#elif GenIs(MINIMAL_INT)
  return("LCG of Park/Miller (int version)");

#elif GenIs(DRAND48)
  return("drand48 from UNIX");

#elif GenIs(GFSR521)
  return("GFSR of length-521 of Ripley 1990");

#elif GenIs(TGFSR)
  return("Matsumoto `twisted' GFSR");

#else
#error "empty function"
#endif
}

/*
 * short text description of the current generator
 */
int
rng_numeric_name(void)
{
  return(GENERATOR);
}

/* == STATIC STORAGE FOR STATE ======================================== */

#if GenIs(GFSR521)
#define GFSR_P 521 /* how far back to look for first bit to XOR */
#define GFSR_Q 32  /* distance from first bit to second bit  */
#define GFSR_L 32  /* word length */
static unsigned state[GFSR_P]; /* holds the shift register */
static unsigned *state_dst; /* oldest value */
static unsigned *state_src; /* GFSR_Q in front of oldest value */

#elif GenIs(RANDOM)
#define RANDOM_REGLENGTH 64
static long state[RANDOM_REGLENGTH];

#elif GenIs(MINIMAL_REAL)
static double state; /* but this double only holds integers */

#elif GenIs(MINIMAL_INT)
static int state;

#elif GenIs(TGFSR)
#define TGFSR_P 624
static unsigned long state[TGFSR_P]; /*  624 words */

#endif
/* == INITIALIZING THE STATE ============================================ */

/* special for the GFSR521 */
#if GenIs(GFSR521)
#define MODULUS_I 2147483647  /* 2^31 - 1*/
#define MULTIPLIER_I 16807    /* 7**5 */
#define MODoverMULT_I 127773  /* MODULUS / MULTIPLIER */
#define MODmodMULT_I 2836     /* MODULUS % MULTIPLIER */

static
void
init_gfsr(int seed)
{
  unsigned char x0[2*GFSR_P]; /* seed bits : 0/1 */
  int  x0_seed;
  int i, j;

  /* set up first GFSR_P entries via linear generator */
  for (x0_seed = seed, i = 0; i < GFSR_P; i++) {
    x0_seed = 
      MULTIPLIER_I * (x0_seed % MODoverMULT_I) - 
      MODmodMULT_I * (x0_seed / MODoverMULT_I);
    if (x0_seed <= 0)
      x0_seed += MODULUS_I;
    x0[i] = (x0_seed > (MODULUS_I/2)); /* convert to 0/1 for x0 */
  }
  /* set up rest of x0 */
  for (i = GFSR_P; i < 2*GFSR_P; i++) 
    x0[i] = x0[i - GFSR_P] ^ x0[i - GFSR_Q];
  /* finally set the state array */
  for (i = 0; i < GFSR_P; i++) 
    for (state[i] = j = 0; j < GFSR_L; j++)
      state[i] |= x0[i+16*(j+1)] << j;
  /* initialize the pointers */
  state_dst = state + (GFSR_P-GFSR_Q-1);
  state_src = state + (GFSR_P-1);
}
#endif

/* special for the TGFSR */
#if GenIs(TGFSR)
static
void
init_tgfsr(unsigned long seed) /* seed should not be 0 */
{
  int k;
  
  /* setting initial state using     */
  /* the generator Line 25 of Table 1 in          */
  /* [KNUTH 1981, The Art of Computer Programming */
  /*    Vol. 2 (2nd Ed.), pp102]                  */

  state[0]= seed & 0xffffffff;
  for (k=1; k<TGFSR_P; k++)
    state[k] = (69069 * state[k-1]) & 0xffffffff;
}
#endif

/*
 * Initialize state according to a given seed
 */
int
rng_init_state_seeded(int seed)
{
#if GenIs(RANDOM)
  initstate((unsigned) seed, (char *) state, sizeof(state));
  setstate((char *) state);

#elif GenIs(MINIMAL_REAL)
  state = seed; /* an assignment is enough */

#elif GenIs(MINIMAL_INT)
  state = seed; /* an assignment is enough */

#elif GenIs(DRAND48)
  srand48((long) seed);

#elif GenIs(GFSR521)
  init_gfsr(seed);

#elif GenIs(TGFSR)
  init_tgfsr((unsigned long) seed);

#else
#error "empty function"
#endif
  return(seed);
}

/* 
 * form a positive int to seed the random number generator 
 */

int 
rng_new_seed(void)
{
  unsigned int t = 0;
  FILE *fp;
  struct timeval tv;
  int ok = 0;
  int nread;

  /* try to read t from /dev/urandom */
  if (!ok) {
    fp = fopen("/dev/urandom", "r");
    if (fp != NULL) {
      nread = fread(&t, sizeof(t), (size_t) 1, fp);
      fclose(fp);
      if (nread == 1)
	ok = 1; // success
      else
	t = 0; // reset to uninitialized
    }
  }
  /* try to open a pipe to the system to get a seed */
  if (!ok) {
    // md5 is less portable, and no better for our purposes
    // plus, it returns a very long checksum
    fp = popen("ps auxww | cksum", "r");
    if (fp != NULL) {
      nread = fscanf(fp, "%u", &t);
      pclose(fp);
      if (nread == 1)
	ok = 1; // success
      else
	t = 0; // reset to uninitialized
    }
  }
  /* if not, use the old initialization method.
   * its seeds can repeat too often, although we have tried
   * to patch that up by using the "microseconds" from tv.
   */
  if (!ok) {
    gettimeofday(&tv, NULL);
    /* "logic" here:
     * tv.sec, tv.usec are time-of-day
     * clock() is time-since-invocation
     * Low-order bits of tv.sec change from run to run
     * Many bits of tv.usec change.
     * So, combine them so the changing bits do not overlap:
     * sec bits low, usec bits high.  Since 1e6 ~ 2^20, there
     * are 20 bits of tv.usec.
     * uid is not doing much.
     * For pids, make sure if the pid and ppid change by one each,
     * the changes do not cancel.
     */
    t = 
      ((tv.tv_sec  & 0xffffffff)  <<  0) ^
      ((tv.tv_usec & 0x000fffff)  << 11) ^
      ((clock()    & 0x00ffffff)  <<  7) ^
      ((getuid()   & 0x000fffff)  << 10) ^
      ((getpid()   & 0x000fffff)  <<  9) ^
      ((getppid()  & 0x0000ffff)  <<  1);
  }
  /* at this point, t is defined */
  t &= 0x7fffffff;   /* force t < 2^31 */
  if (t == 0) t = 1; /* force nonzero */
  return(t);
}


/*
 * Set the seed automatically.  Does not re-seed the rng if
 * it is already seeded.  This is useful if this is used as a 
 * dynamically loaded library: multiple uses of this routine will
 * not reset the seed.  Returns the seed used.
 *
 * assumes the generated seed is nonzero.
 */
int
rng_init_state(void)
{
  static int seed = 0; /* 0 is also the sentinel for unseeded */

  /* do not re-seed */
  if (seed == 0) {
    /* returned seed is guaranteed nonzero */
    rng_init_state_seeded(seed = rng_new_seed());
  }
  return seed;
}


/* == SAVING/RESTORING THE STATE ======================================== */

/*
 * number of bytes in the entire state of the generator
 */
int
rng_state_length(void)
{
  int bufused; /* sizeof the space needed to store the state */

#if GenIs(RANDOM)
  bufused = sizeof(state);

#elif GenIs(MINIMAL_REAL)
  bufused = sizeof(state);

#elif GenIs(MINIMAL_INT)
  bufused = sizeof(state);

#elif GenIs(DRAND48)
  bufused = 3*sizeof(unsigned short);

#elif GenIs(GFSR521)
  bufused = sizeof(state) + sizeof(state_dst) + sizeof(state_src);

#elif GenIs(TGFSR)
  bufused = sizeof(state);

#else
#error "empty function"
#endif
  return(bufused);
}

/*
 * get the entire state of the generator
 * inputs: buffer for the state, and its length
 * returns size of state, or 0 if failure (state too big for buf)
 */
int
rng_get_state(char *buf, int buflen)
{
  /* sizeof the space needed to store the seed */
  int bufused = rng_state_length();

#if GenIs(RANDOM)
  if (buflen >= bufused)
    bcopy((void *) state, (void *) buf, (size_t) bufused);
  else
    bufused = 0;

#elif GenIs(MINIMAL_REAL)
  if (buflen >= bufused)
    bcopy((void *) &state, (void *) buf, (size_t) bufused);
  else
    bufused = 0;

#elif GenIs(MINIMAL_INT)
  if (buflen >= bufused)
    bcopy((void *) &state, (void *) buf, (size_t) bufused);
  else
    bufused = 0;

#elif GenIs(DRAND48)
  if (buflen >= bufused) {
    unsigned short *buf_internal;

    /* get address of internal buffer where state is located */
    /* (as side-effect: sets state to contents of buf!) */
    buf_internal = seed48((unsigned short *) buf); 
    /* copy the internal buffer to our own buffer */
    bcopy((void *) buf_internal, (void *) buf, bufused);
    /* restore the old state, undoing side-effect above */
    seed48((unsigned short *) buf); 
  }
  else
    bufused = 0;

#elif GenIs(GFSR521)
  if (buflen >= bufused) {
    /*  save: state, state_dst, state_src */
    bcopy((void *) state, 
	  (void *) buf, 
	  (size_t) sizeof(state));
    bcopy((void *) &state_dst, 
	  (void *) (buf + sizeof(state)), 
	  (size_t) sizeof(state_dst));
    bcopy((void *) &state_src, 
	  (void *) (buf + sizeof(state) + sizeof(state_dst)), 
	  (size_t) sizeof(state_src));
  }
  else
    bufused = 0;

#elif GenIs(TGFSR)
  if (buflen >= bufused)
    bcopy((void *) state, (void *) buf, (size_t) bufused);
  else
    bufused = 0;

#else
#error "empty function"
#endif
  return(bufused);
}

/*
 * get the entire state to buf
 * inputs: buffer holding the state, and its length
 * returns size of state, or 0 if failure (buf too small for state)
 */
int
rng_set_state(char *buf, int buflen)
{
  int bufneed = rng_state_length();

#if GenIs(RANDOM)
  if (buflen >= bufneed) {
    bcopy((void *) buf, (void *) state, (size_t) bufneed);
    setstate((char *) state); /* tell random() to use this buffer */
  }
  else
    bufneed = 0;

#elif GenIs(MINIMAL_REAL)
  if (buflen >= bufneed)
    bcopy((void *) buf, (void *) &state, (size_t) bufneed);
  else
    bufneed = 0;

#elif GenIs(MINIMAL_INT)
  if (buflen >= bufneed)
    bcopy((void *) buf, (void *) &state, (size_t) bufneed);
  else
    bufneed = 0;

#elif GenIs(DRAND48)
  if (buflen >= bufneed) {
    /* copy our own buffer to the internal buffer */
    seed48((unsigned short *) buf); 
  }
  else
    bufneed = 0;

#elif GenIs(GFSR521)
  if (buflen >= bufneed) {
    /*  restore: state, state_dst, state_src */
    bcopy((void *) buf, (void *) state, (size_t) sizeof(state));
    bcopy((void *) (buf + sizeof(state)), 
	  (void *) &state_dst, 
	  (size_t) sizeof(state_dst));
    bcopy((void *) (buf + sizeof(state) + sizeof(state_dst)), 
	  (void *) &state_src, 
	  (size_t) sizeof(state_src));
  }
  else
    bufneed = 0;

#elif GenIs(TGFSR)
  if (buflen >= bufneed)
    bcopy((void *) buf, (void *) state, (size_t) bufneed);
  else
    bufneed = 0;

#else
#error "empty function"
#endif
  return(bufneed);
}

/* == GENERATING RANDOM NUMBERS ======================================== */

#if GenIs(RANDOM)
#define MODULUS 2147483648.0  /* 2^31 */
#define MODULUS_INV 4.656612873077393e-10 /* 1/MODULUS */

#elif GenIs(MINIMAL_REAL)
#define MODULUS 2147483647.0  /* 2^31 - 1*/
#define MODULUS_INV 4.656612875245797e-10 /* 1/MODULUS */
#define MULTIPLIER 16807.0 /* 7**5 */

#elif GenIs(MINIMAL_INT)
#define MODULUS 2147483647  /* 2^31 - 1*/
#define MODULUS_INV 4.656612875245797e-10 /* 1/MODULUS */
#define MULTIPLIER 16807 /* 7**5 */
#define MODoverMULT 127773 /* MODULUS / MULTIPLIER */
#define MODmodMULT 2836 /* MODULUS % MULTIPLIER */

#elif GenIs(GFSR521)
#define MODULUS 4294967296 /* 2^32 */
#define MODULUS_INV 2.3283064365387e-10 /* 1/MODULUS */

#elif GenIs(TGFSR)
/* Period parameters */  
#define N_twist 624
#define M_twist 397
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */
/* for tempering */   
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

#endif

/*
 * uniform number in (0,1)
 */
double 
rng_uniform (void)
{
#if GenIs(RANDOM)
  return(((double) random()) * MODULUS_INV);

#elif GenIs(MINIMAL_REAL)
  return(
	 (state = fmod(MULTIPLIER * state, MODULUS)) * MODULUS_INV
	 ); 
#elif GenIs(MINIMAL_INT)
  state = 
    MULTIPLIER * (state % MODoverMULT) - 
    MODmodMULT * (state / MODoverMULT);
  if (state <= 0)
    state += MODULUS;
  return(state * MODULUS_INV); 

#elif GenIs(DRAND48)
  return(drand48());

#elif GenIs(GFSR521)
  double rval; 

  rval = *state_dst * MODULUS_INV; /* use this value */
  *state_dst ^= *state_src; /* update the shift register */
  state_dst = (state_dst > state) ? state_dst-1 : state+(GFSR_P-1);
  state_src = (state_src > state) ? state_src-1 : state+(GFSR_P-1);
  return(rval);

#elif GenIs(TGFSR)
  unsigned long y;
  static int k = 1;
  static unsigned long mag01[2]={0x0, MATRIX_A};
  /* mag01[x] = x * MATRIX_A  for x=0,1 */
  
  if(k == N_twist){ /* generate N_twist words at one time */
    int kk;
    for (kk=0;kk<N_twist-M_twist;kk++) {
      y = (state[kk]&UPPER_MASK)|(state[kk+1]&LOWER_MASK);
      state[kk] = state[kk+M_twist] ^ (y >> 1) ^ mag01[y & 0x1];
    }
    for (;kk<N_twist-1;kk++) {
      y = (state[kk]&UPPER_MASK)|(state[kk+1]&LOWER_MASK);
      state[kk] = state[kk+(M_twist-N_twist)] ^ (y >> 1) ^ mag01[y & 0x1];
    }
    y = (state[N_twist-1]&UPPER_MASK)|(state[0]&LOWER_MASK);
    state[N_twist-1] = state[M_twist-1] ^ (y >> 1) ^ mag01[y & 0x1];
    
    k = 0;
  }
  
  y = state[k++];
  y ^= TEMPERING_SHIFT_U(y);
  y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
  y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
  y &= 0xffffffff; /* you may delete this line if word size = 32 */
  y ^= TEMPERING_SHIFT_L(y);

  return ( (double)y * (1.0/ (unsigned long)0xffffffff) );

#else
#error "empty function"
#endif
}

/******************************************************************
 *
 * Normal random numbers
 *
 ******************************************************************/

/*
 * normal, mean 0, unit variance
 * 
 * formerly just:
 *  return(sqrt(-2 * log(rng_uniform())) * sin(2 * M_PI * rng_uniform()));
 * as in rng_normal_old().
 * In this version, we retain the "orthogonal" random number as a double
 * called ace_in_the_hole, and indicate its presence (on alternate calls
 * to this routine) via have_ace.  
 * The method used eliminates the sin, cos calls of the "pure" bivariate 
 * transformation by starting from a uniform tuple (u1,u2) that lives 
 * on the unit disk.  The random angle that used to be generated by 
 * sin(2 pi random()) is now just encoded in the coordinates of 
 * the uniform tuple, scaled by the square root of its magnitude, 
 * i.e. u1/sqrt(norm) and u2/sqrt(norm).  
 * And,  conditioned to be < 1, norm := u1^2 + u2^2 is U(0,1).  So norm 
 * is used for the length of the bivariate normal, which is generated by:
 *   n1 = sqrt(-2 * log(norm)) * (u1/sqrt(norm))
 *   n2 = sqrt(-2 * log(norm)) * (u2/sqrt(norm))
 * Tucking the sqrt(norm) into the first factor yields the code below.
 * It runs about 2.5 times faster than the old version.
 */
double
rng_normal(void) {
  static double ace_in_the_hole; /* a N(0,1) side-effect of the prior call */
  static int have_ace = 0;       /* 1 if we have the above, 0 if not */
  double u1, u2;                 /* a pair of uniform variables */
  double norm;                   /* their norm taken as a tuple */

  if (have_ace) {
    have_ace = 0;                /* we're about to use it up */
    return ace_in_the_hole;
  } else {
    do {
      /* get a pair (u1,u2) that is uniform on the unit disk */
      u1 = 2.0*rng_uniform() - 1.0;  /* U([-1,1)) */
      u2 = 2.0*rng_uniform() - 1.0;  /* U([-1,1)) */
    } while (1.0 <= (norm = u1*u1 + u2*u2));
    norm = sqrt(-2.0 * log(norm)/norm); /* just a scale factor, see above */
    ace_in_the_hole = u1 * norm; /* first N(0,1): save for later */
    have_ace = 1;                /* remember we saved it */
    return u2 * norm;            /* second N(0,1): return it now */
  }
}

/*
 * normal, mean 0, unit variance
 */
double
rng_normal_old(void) {
    return (sqrt (-2.0 * log(rng_uniform())) *
	    sin(2.0 * M_PI * rng_uniform()));
}


/*
 * normal, mean mu, standard deviation sigma
 */
double
rng_normal_sdev(double mu, double sigma) {
  return rng_normal()*sigma + mu;
}

/*
 * normal, mean mu, variance sigma2
 * if rng_normal were more clever, we could eliminate the sqrt()
 */
double
rng_normal_var(double mu, double sigma2) {
  return rng_normal()*sqrt(sigma2) + mu;
}

/******************************************************************
 *
 * Non-normal real-valued random numbers
 *
 ******************************************************************/

/*
 * rng_gamma(alpha,theta): Returns a draw from gamma(alpha,theta).
 * 
 * The corresponding density is:
 * 
 *                       x^(alpha - 1) exp(-x/theta)
 * p(x; alpha, theta) =  ---------------------------
 *                       Gamma(alpha) theta^alpha
 * 
 * In this parameterization, theta is purely a scale parameter 
 * (similar to sigma in the Normal distribution) and alpha is 
 * a shape parameter.
 * 
 * For alpha <= 1, the distribution is one-sided and cusped like 
 * an exponential, and for alpha > 1, the distribution is (increasingly) 
 * two-sided and rounded like a gaussian.
 * 
 * If X ~ gamma(alpha, theta), then
 * 
 *   E X   = alpha theta
 *   var X = alpha theta^2
 * 
 * The algorithm below is primarily comprised of two rejection-based 
 * methods, recommended by Devroye, one for alpha < 1 and one for
 * alpha > 1.  (We use the exponential for alpha = 1.)  Both are 
 * crafted to have bounded (and small) rejection proportion for 
 * all values of alpha and of the variate x.  
 * I further separate out the alpha = 2 case because
 * it is easy and relatively common.
 * Regarding the alpha = 1 case, Devroye's book seems to indicate 
 * his version of the Best (1978) "XG" algorithm works for alpha = 1.
 * It has a problem there, as his web errata point out.  This is 
 * another reason to special-case the exponential.
 * 
 * Some of the acceptance crafting seems more useful for older hardware. 
 * For example, extra tests have been added for early rejection to avoid
 * computation of exponentials.  The exponentials give the exact rejection
 * criterion, but the extra tests are quicker to compute.  I have not 
 * checked the impact of tweaking the algorithms.
 * 
 * References: 
 *   Luc Devroye, Non-Uniform Random Variate Generation, Springer, 
 *     1986, chap. IX.3-IX.6.
 *   http://cg.scs.carleton.ca/~luc/rnbookindex.html
 *   http://cg.scs.carleton.ca/~luc/errors.pdf
 * 
 * which cites:
 * 
 *   Best 1978: D. J. Best, "A simple algorithm for the computer 
 *   generation of
 *   random samples from a Student's t or symmetric beta distribution," in
 *   Compstat 1978: Proceedings in Computational statistics, ed. L.C.A. 
 *   Corsten and J. Hermans, pp. 341-7, Physica Verlag, Vienna, 1978.
 *   => a.k.a. "XG" algorithm for alpha > 1.
 *   See Devroye, page 410.
 * 
 * and
 * 
 *   Best 1983: D. J. Best, "A note on gamma variate generators with shape
 *   parameter less than unity," Computing, vol. 30, pp. 185-8, 1983.
 *   => a.k.a. "XGS" algorithm for alpha < 1.
 *   See Devroye's exercise 6 on page 426.
 * 
 * I tested this algorithm for proper moments (1 thru 4), and 
 * distributional shape, for many values of alpha (excercising all 
 * cases below).
 *  
 */
double 
rng_gamma(double alpha, double theta)
{
  // follow DeVroye's notation
  double X, U, V, W, Y, Z;  // random variates
  double b, c, t;           // parameters

  if (alpha <= 0 || theta < 0)
    return rng_getnand();
  if (alpha == 1.0) {
    // Use the exponential, e.g. Devroye p. 405
    X = -log(rng_uniform());
  } else if (alpha == 2.0) {
    // Use the sum of 2 exponentials (Devroye p. 405)
    X = -log(rng_uniform() * rng_uniform());
  } else if (alpha < 1.0) {
    // Use RGS algorithm (Best, 1983; Devroye, p. 426)
    c = 1.0 / alpha; 
    t = 0.07 + 0.75 * sqrt(1.0 - alpha); // breakpoint: min reject prob
    b = 1.0 + exp(-t) * alpha / t;
    // generate a variate x
    while (1) {
      U = rng_uniform(); 
      W = rng_uniform(); 
      V = b * U;
      if (V <= 1.0) {
	X = t * pow(V, c);
	/* the first step is an early bailout because
	 * exp(-x) >= (2-x)/(2+x) for x>=0
	 */
	// test is equivalent to w <= (2-x)/(2+x) since x>=0 
	if (W * (2.0 + X) <= (2.0 - X)) break;
	if (W <= exp(-X)) break;
      } else {
	X = -log(c * t * (b - V));
	Y = X / t;
	if ((W * (alpha + Y - alpha*Y)) <= 1.0) break;
	if (W <= pow(Y, alpha - 1.0)) break;
      }
    }
  } else {
    // Use rejection algorithm XG (Best, 1978; Devroye, p. 410)
    //   (This algorithm is OK for alpha > 1.  If alpha = 1,
    //   b = 0, and the last rejection step with log(X/b) is
    //   illegal.  This is noted in Devroye's errata on the web.)
    b = alpha - 1.0;
    c = 3.0 * alpha - 0.75;
    // generate a variate x
    while (1) {
      U = rng_uniform();  
      V = rng_uniform();
      W = U * (1.0 - U);  
      Y = sqrt(c / W) * (U - 0.5);
      X = b + Y;
      if (X >= 0.0) {
	Z = 64.0 * W*W*W * V*V;
	/* the first step is an early bailout which
	 * Devroye remarks can be omitted at little cost.
	 */
	// test is same as Z <= 1 - 2*Y*Y/X since X >= 0 here
	if ((Z * X) <= (X - 2.0 * Y * Y)) break;
	if (log(Z) <= (2.0 * (b * log(X/b) - Y))) break;
      }
    }
  }
  return X * theta;
}

/*
 * beta(a,b) sample: Just use the gamma generator above.
 * Based on the remarks of Devroye chap. IX.4, pp 432-433:
 * "Roughly speaking, we will be able to improve
 * over this generator by at most 50%."
 */
double 
rng_beta(double a, double b)
{
  double Ga = rng_gamma(a, 1.0);
  double Gb = rng_gamma(b, 1.0);
  return Ga / (Ga + Gb);
  // note: Ga+Gb is independent of the returned value, and is
  // distributed as gamma(a+b).
}

/*
 * chi-square(r) sample: Use the gamma generator above.
 * See Devroye, chap. IX.3, pg 403.
 */
double 
rng_chi2(double r)
{
  return rng_gamma(r * 0.5, 2.0);
}


/*
 * exponential(1) sample: Simple inversion method.
 * If you want a different rate, just multiply.
 */
double 
rng_exponential()
{
  return -log(rng_uniform());
}



/******************************************************************
 *
 * Discrete-valued random numbers
 *
 ******************************************************************/

/******************************************************************
 * Poisson
 */

/*
 * return a Poisson r.v., using an inter-arrival time method
 * I.e., a Poisson process of unit intensity has exponential(1)
 * inter-arrival times, thus the number of arrivals between
 * time = 0 and time = lambda is Poisson(lambda), thus generate
 * Poisson R.V.'s by seeing how many exponential(1) r.v.'s it
 * takes until the sum exceeds lambda.  Since an exponential(1)
 * is just -log(uniform), this is the same as seeing when:
 *   exp(Sum_k [ -log(U(k)) ] ) > lambda, 
 * i.e. when 
 *   Prod_k(U(k)) < lambda.
 * See also Devroye, chapter X.3, Lemma 3.3.
 */
static
double
rng_poisson_count(double lambda)
{
  const double gate = exp(-lambda);
  double prod, x;

  if (gate == 1) {
    // e.g., lambda = 1e-20 => exp(-lambda) == 1 for IEEE doubles
    // The below line will always return 0 with most rng_uniform
    // implementations, because they don't output numbers less than
    // floating point epsilon.
    return (rng_uniform() < lambda) ? 1.0 : 0.0;
  }
  // it is guaranteed that rate < 1, so at exit, x is at least 1
  for (prod = 1.0, x = 0.0; prod > gate; prod *= rng_uniform(), x++)
    ;
  return x-1.0; // the previous "prod" was the last one > gate
}

/*
 * transfomed rejection algorithm, due to W. Hörmann.
 *
 * This is algorithm PTRS of the paper:
 *   W. Hörmann, 
 *   "The transformed rejection method for generating Poisson random 
 *   variables," Insurance: Mathematics and Economics 12, 39-45 (1993).
 *
 * See also:
 *   http://statmath.wu.ac.at/staff/hoermann/publications.html
 *
 * This generator is faster than that of Devroye, less complex,
 * and does not seem to suffer from unpleasant bugs!  Checked
 * 2/2010 against the pmf for several rates mu.
 *
 * Valid only for rate (mu) >= 10.
 */
static
double
rng_poisson_trans_reject(double mu)
{
  double b, a, v_r;
  double one_alpha;  // "(1/alpha)" in the paper
  double U, US, V;
  double k;
  
  // find parameters
  b = 0.931 + 2.53 * sqrt(mu);
  a = -0.059 + 0.02483 * b;
  v_r = 0.9277 - 3.6224 / (b - 2.0);
  // loop until accept
  while (1) {
    U = rng_uniform();
    V = rng_uniform();
    U = U - 0.5;
    US = 0.5 - fabs(U);
    k = floor((2.0*a/US + b) * U + mu + 0.43);
    if (k < 0)
      continue;
    if (US >= 0.07 && V <= v_r)
      break;
    if (US < 0.013 && V > US)
      continue;
    // further setup
    one_alpha = 1.1239 + 1.1328 / (b - 3.4);
    // acceptance-rejection test
    if (log(V * one_alpha / (a/(US*US) + b))
	<= (-mu + k*log(mu) - lgamma(k+1.0)))
      break;
  }
  return k;
}


double
rng_poisson(double lambda)
{
  // 30 is the approximate switchover point (intel core 2, 2010)
  const double mode_switch = 30;

  if (lambda <= 0.0)
    return rng_getnand(); // illegal
  else if (lambda < mode_switch)
    // fast for small lambda
    //   note, exp(-mode_switch) must comfortably fit in a double
    return rng_poisson_count(lambda);
  else
    // longer setup time, but better for large lambda
    return rng_poisson_trans_reject(lambda);
}

/******************************************************************
 * Geometric
 */

/*
 * geometric(p) sample: Simple inversion method.
 * See Devroye chap. X.2.
 */
double 
rng_geometric(double p)
{
  if (p < 0.25)
    // avoid loss of precision in computing log(1-p) 
    // when p is small by computing:
    //    log1p(-p) = log(1 + (-p)) = log(1-p)
    return ceil(log(rng_uniform()) / log1p(-p));
  else
    return ceil(log(rng_uniform()) / log(1.0 - p));
}


/******************************************************************
 * Binomial
 */

/*
 * Waiting time algorithm.  Works for small n*p.
 * Idea is that a geometric(p) r.v. is the number of 
 * Bernoulli trials up to and including the first success.
 * After the total of trials exceeds n, we count the number 
 * of summands.  This is one more (because it exceeded n) than
 * the number of 1's in the corresponding Bernoulli sequence.
 * The latter is the binomial(p) draw.
 * See Devroye, "First waiting time algorithm", p. 525.
 */

static
double
rng_binomial_wait(double n, double p)
{
  double sum, x;

  // take care of endpoints now, for robustness
  if (n == 0.0)
    return 0.0;
  else if (p == 0.0)
    return 0.0;
  else if (p == 1.0)
    return n;
  // wait for geometric
  for (sum = 0, x = 0; sum <= n; sum += rng_geometric(p), x++)
    ;
  return x - 1.0;
}

/*
 * This is algorithm BTRS from Hörmann.
 * It requires n*p > 10.
 * This generator is faster than that of Devroye, less complex,
 * and (once again) does not seem to suffer from unpleasant bugs!  
 * Checked extensively against binomial probabilities from
 * Matlab, 2/2010.
 *
 * W. Hörmann, :The generation of binomial random variates,"
 * Journal of Statistical Computation and Simulation 46, 
 * 101-110 (1993), online in TR form at:
 * http://statmath.wu.ac.at/staff/hoermann/publications.html
 *
 * The tech report version has a missing log() in step
 * 3.1 of the algorithm; noted below.  The published paper
 * may correct this.  One can also determine this rather blatant
 * error by comparison with the more elaborate BTRD algorithm,
 * which uses the same v <- log(v) construction in step 3.2, 
 * before the final acceptance test in step 3.4 in which the 
 * new "v" is compared against the log-probability.  Thus, the
 * log needs to be there.
 */
static
double
rng_binomial_trans_reject(double n, double p)
{
  const double spq = sqrt(n*p*(1-p));
  double b, a, c, v_r;
  double alpha, lpq, m, h;
  double U, US, V;
  double k;
  
  // find parameters
  b = 1.15 + 2.53 * spq;
  a = -0.0873 + 0.0248 * b + 0.01 * p;
  c = n * p + 0.5;
  v_r = 0.92 - 4.2 / b;
  // loop until accept
  while (1) {
    U = rng_uniform();
    V = rng_uniform();
    U = U - 0.5;
    US = 0.5 - fabs(U);
    k = floor((2.0 * a/US + b) * U + c);
    if (k < 0 || k > n)
      continue;
    if (US >= 0.07 && V <= v_r)
      break;
    // further setup
    alpha = (2.83 + 5.1/b) * spq;
    lpq = log(p/(1-p));
    m = floor(n*p+p);
    h = lgamma(m+1) + lgamma(n-m+1);
    // acceptance-rejection test
    V = log(V * alpha/(a / (US*US) + b)); // log() missing from TR
    if (V <= h - lgamma(k+1) -lgamma(n-k+1) + (k-m)*lpq)
      break;
  }
  return k;
}


/*
 * final binomial random variable entry point.
 * 
 * Calls one of two RNG's, a simple one for small n*p, and
 * a more sophisticated one with larger setup time for
 * larger n*p.
 */

double
rng_binomial(double n, double p)
{
  // 16 is the approximate switchover point (intel core 2, 2010)
  // (must be > 10 anyway to accommodate restrictions on fancy RNG)
  const double mode_switch = 16;
  int flip;
  double X;    // random draw

  // take care of errors now
  //   (tests written to yield error for nans and infinites)
  if (!(n >= 0.0) || !(floor(n) == n) || !finite(n))
    return rng_getnand();
  if (!(p >= 0.0 && p <= 1.0))
    return rng_getnand();
  // take care of endpoints to avoid trouble later
  if (n == 0.0)
    return 0.0;
  else if (p == 0.0)
    return 0.0;
  else if (p == 1.0)
    return n;
  // if p > 1/2, flip to: n - bin(n, 1-p) 
  if (p > 0.5) {
    p = 1.0 - p;
    flip = 1;
  } else {
    flip = 0;
  }
  // generate bin(n, p)
  if ((n*p) < mode_switch) {
    // fast for small rates
    X = rng_binomial_wait(n, p);
  } else {
    // np is large: transformed rejection method
    X = rng_binomial_trans_reject(n, p);
  }
  // account for flip if need be
  return flip ? (n - X) : X;
}

/******************************************************************
 *
 * Linear algebra tools
 *
 ******************************************************************/

/* 
 * Cholesky decomposition (in-place).
 * No pivoting.  Works only on strictly positive-definite matrices 
 * for this reason.  
 * Checked against Matlab chol and found agreement to floating-point
 * precision, 4/2006.
 * We write  R = G' G, where G is the Cholesky decomposition.
 * References only the diagonal and the "upper triangle" of the input R, 
 * and over-writes it with the cholesky factor G.  The upper half of 
 * the input matrix is left alone.
 * For 0 <= j <= i <= N-1, mat[i*N+j] contains the Cholesky factor.  If we 
 * identify G(i,j) with mat[i*d+j], the "lower triangle" corresponds to 
 * G(i,j) for j <= i, and the other entries in G are zero.
 * About shape:
 *   The C Cholesky decomposition below produces a triangular
 *   Cholesky factor.  It is laid out in memory the same way
 *   as a Cholesky factor from matlab's "chol(R)".  
 *
 * tested versus matlab chol, turmon, 6/2006; again 2/2010.
 * (randomized inputs, size 2x2, 5x5, 12x12, 128x128, 256x256)
 *
 * 0 is returned on failure (ie, non pd matrix), 1 for success.
 */

static
int
rng_cholesky(double *a, 
              int N)
{
  double fac;
  int i, j, k;

  for (k = 0; k < N; k++) {
    if (a[k*N+k] <= 0.0) 
      return 0;
    fac = a[k*N+k] = sqrt(a[k*N+k]);
    fac = 1/fac;
    for (j = k+1; j < N; j++)
      a[j*N+k] *= fac;
    for (j = k+1; j < N; j++)
      for (i = j; i < N; i++)
        a[i*N+j] -= a[i*N+k] * a[j*N+k];
  }
  return 1;
}

/*
 * Find the inverse of a nonsingular upper-triangular nXn matrix a.
 * The computation is done in place, so the inverse replaces a.
 * The lower triangle is not referenced or used.
 * If a does not have non-zero diagonals, the inverse does not
 * exist, and 0 is returned.  Otherwise, 1 is returned.
 *
 * This is included here because it would help with an
 * inverse Wishart, if we ever need it.
 *
 * tested versus matlab inv, turmon, 6/2006 
 * (randomized inputs, size 2x2, 5x5, 12x12, 128x128, 256x256)
 *
 */

// not currently used!
#ifdef NOT_DEFINED

static
int
rng_invert_upper(double *a,
		 int N)
{
  int i,j,k;
  double sum;

  for (j = 0; j < N; j++) {
    if (a[j*N+j] == 0.0)
      return 0; /* failure */
    a[j*N+j] = 1 / a[j*N+j];
    for (i = j+1; i < N; i++) {
      sum = 0.0;
      for (k = j; k < i; k++) {
        sum -= a[i*N+k] * a[k*N+j];
      }
      a[i*N+j] = sum / a[i*N+i];
    }
  }
  return 1;
}

#endif


/* multiply c = a * b, with all square upper-triangular matrices,
 * and all elements stored in standard matlab ordering
 */
static
void 
rng_mult_upper_triangles(double* c, double* a, double* b, int N)
{
  int i, j, k;

  for (i = 0; i < N; i++)
    for (j = i; j < N; j++) {
      // zero out c(j,i) (lower triangle, or diagonal)
      c[j+i*N] = 0;
      // add up c(i,j) (j>=i)
      c[i+j*N] = 0;
      for (k = i; k <= j; k++)
	// Computing:
        //   c(i,j) += a(i,k) * b(k,j)
	// Note:
	//   a(i,k) = 0 if k < i, giving lower loop bound
	//   b(k,j) = 0 if k > j, giving upper loop bound
        c[i+j*N] += a[i+k*N] * b[k+j*N];
    }
}


/* multiply b = a * a', with a square and upper triangular,
 * and all elements stored in standard matlab ordering.
 * The result is a full matrix.
 */
static
void 
rng_square_upper_triangle(double* b, double* a, int N)
{
  int i, j, k;

  for (i = 0; i < N; i++)
    for (j = i; j < N; j++) {
      // sum up b(i,j) = b(j,i) (note j >= i)
      b[i+j*N] = 0;
      for (k = 0; k <= i; k++)
	// Note:
	//   a(k,i) = 0 if k > i, and
	//   a(k,j) = 0 if k > j, giving k-loop lower bound: min(i,j)
	// since j >= i, we start at i.
        b[i+j*N] += a[k+i*N] * a[k+j*N];
      // reflect into b(j,i)
      b[j+i*N] = b[i+j*N];
    }
}

/*************************************************************
 *
 * Normal random vectors
 *
 *************************************************************/

/*
 * normal d-vector, mean mu (d), covariance sigma (dxd)
 *
 * If mu is supplied as NULL, assume mu is all zeros.
 * 
 * turmon 4/2006: tested for d=1 and d=6 with dense random covariances.
 * First, second (full covariance), and fourth (elementwise) moments check
 * out, and show correct convergence to their specified values as n grows.
 * 3D scatter plots (of the 6D vectors) look normal.  1-d histograms
 * also look normal.  The Kolmogorov-Smirnov distributional test along
 * each of 6 dimensions shows close adherence to normality (p-values
 * around 0.1-0.5 at n=1e6).
 *
 * Regarding strides.  stridex1 is the distance in doubles between
 * adjacent vector samples (typically equals d).  stridex2 is the
 * distance in doubles between adjacent elements within one vector
 * (typically equals 1).  stridemu is the same, between elements
 * of mu (typically equals 1).
 */
int 
rng_normal_vector(
		  double* x,     /* samples */
		  double* mu,    /* mean, taken to be 0.0 if NULL */
		  double* sigma, /* covariance */
		  int n,         /* number of samples */
		  int d,         /* dims of mean, samples */
		  int stridex1,  /* stride for each sample */
		  int stridex2,  /* stride within each sample */
		  int stridemu   /* stride for elements in mu */
		  )
{
  int i, j, k; /* counters */
  double *z;   /* temp space for N(0,1)'s */
  double *R;   /* temp space for cholesky factor */

  /* allocate temporary storage */
  R = calloc(d*d, sizeof(double));
  z = calloc(d  , sizeof(double));
  if (R==NULL || z==NULL) {
    (void)fprintf(stderr, "rng_normal_vector2: failed calloc\n");
    return 0;     
  }  
  /* copy sigma into R */
  memcpy((void *) R, (void *) sigma, sizeof(double)*d*d);
  /* find cholesky factor */
  if (!rng_cholesky(R, d)) {
    free(R); free(z); return(0);
  }
  /* now, R[i*d + j], j<=i, contains the Cholesky factor */
  /* sample each x(i), 1 <= i <= n */
  for (i = 0; i < n; i++) {
    /* make a fresh z for each x(i) */
    for (j = 0; j < d; j++) 
      z[j] = rng_normal(); /* z is d, iid, N(0,1) RV's */
    /* compute x(i) = mu + R * z  */
    /* note that R is lower triangular */
    for (j = 0; j < d; j++) {
      /* x(i,j) = mu(j) */
      x[i*stridex1+j*stridex2] = mu ? mu[j*stridemu] : 0.0;
      /* perform x(i,j) += R(j,1:j) <dot> z(1:j) */
      for (k = 0; k <= j; k++) 
	x[i*stridex1+j*stridex2] += R[j*d+k] * z[k];
    }
  }
  free((void *) R); free((void *) z); 
  return(1); /* OK */
}


/*************************************************************
 *
 * Wishart matrices
 *
 *************************************************************/

/*
  Wishart Random Matrices: Background
  
  References:
    James E. Gentle, "Random Number Generation and Monte Carlo Methods,"
    Springer, 1998.
    See the description surrounding Algorithm 5.8.  We use this idea.
    But, the notation in step 3 of that algorithm does not make clear
    that simple matrix multiplication is what's happening.

    Mark E. Johnson, "Multivariate statistical simulation," 
    Wiley, 1987.
    The algorithm on p. 204 has clearer notation which we follow.

    T. W. Anderson, "An introduction to multivariate statistical
    analysis," Wiley, 1984.
    The ideas are on pages 249-251.

    W. B. Smith and R. R. Hocking, R. R. "Algorithm AS 53: 
    Wishart Variate Generator". JRSS C, 1972, 21 (3): 341-345. 
    This paper apparently is the first modern reference for the 
    simplified algorithm.

 
  turmon 2/2010: checked against Matlab wishrnd for agreement
  with wishrnd from the statistics toolbox.  Found good agreement
  for randomized sigma, both in expected value of ensembles of
  ~10^4 draws of W, and in standard deviation of the ensembles of W's.
*/

/* 
 * Generate NxN Wishart matrix with dof degrees of freedom.
 * This is equivalent to an outer product of "dof" random vectors
 * from N(0,sigma), but done using only NxN matrix operations using
 * the "Bartlett decomposition" referred to above.
 * The idea is as follows.  Suppose:
 *    T = [upper triangular, as described below]
 * A Wishart matrix for sigma = I is, according to the Bartlett
 * decomposition,
 *    W = T T'
 * Suppose we have the upper-triangular cholesky factorization:
 *    sigma = G' G
 * A Wishart matrix for general sigma is
 *    W = G' T' T G = (T G)' * (T G) = Z' Z
 * where 
 *    Z = T G.
 * All operations are at at most O(N^3), which is better
 * than the outer product method of O(N^2 * dof) if dof >> N.
 */

static
int
rng_wishart_chol(double *sigma, double dof, int N, double *w)
{
  int i,j;
  double *G = NULL;  // cholesky factor
  double *T = NULL;  // triangular random matrix for Bartlett decomp.
  double *Z = NULL;  // workspace
  int ok = 0;        // suppose not OK

  /* create temporary storage */
  G = (double *) calloc(N*N, sizeof(double));
  T = (double *) calloc(N*N, sizeof(double));
  Z = (double *) calloc(N*N, sizeof(double));
  if (!G || !T || !Z)
    goto done; // ok=0, so will free and return failure

  /* get upper triangle and diagonal values:
   *   upper triangle: iid ~N(0,1)
   *   diagonal elements ~ Chi2(dof) ... Chi2(dof-N+1)
   */
  for (j = 0; j < N; j++)
    for (i = 0; i < N; i++)
      if (i == j)
	T[i+N*i] = sqrt(rng_chi2(dof - ((double) i)));
      else if (i < j)
	T[i+N*j] = rng_normal();
      else
	T[i+N*j] = 0.0;

  // G = chol(sigma)
  bcopy(sigma, G, N*N*sizeof(double)); // in-place cholesky
  if (!rng_cholesky(G, N))
    goto done;

  // Z = T * G.  Everything is upper-triangular.
  rng_mult_upper_triangles(Z, T, G, N);

  // w = Z * Z'
  rng_square_upper_triangle(w, Z, N);

  // mxt_put_matrix("S", -1, sigma, N, N);
  // mxt_put_matrix("G", -1, G, N, N);
  // mxt_put_matrix("T", -1, T, N, N);
  // mxt_put_matrix("Z", -1, Z, N, N);

  ok = 1; // successful exit
 done:
  free(G);
  free(T);
  free(Z);
  return ok;
}

/*
 * rng_wishart_outer: Wishart sample "w" using outer products
 * Generates the Wishart matrix directly, using outer products
 * of normal random vectors drawn from N(0,sigma).
 * Returns 0 on failure, 1 on success.
 * 
 * Best for dof not much more than N, but any dof is legal.
 */
static
int
rng_wishart_outer(double *sigma, double dof, int N, double *w)
{
  double *x, *x1;    /* indexes into random draw */
  int i, j, s;
  int ok;

  /* create temporary storage */
  if ((x = (double *) calloc(N*dof, sizeof(double))) == NULL)
    return 0;
  /* make normal RV's */
  ok = rng_normal_vector(x,     /* samples */
			 NULL,    /* mean, taken to be 0.0 if NULL */
			 sigma,   /* covariance */
			 dof,     /* number of samples */
			 N,       /* dims of mean, samples */
			 N,       /* stride for each sample */
			 1,       /* stride within each sample */
			 1);      /* stride for elements in mu */
  if (!ok) {
    free(x);
    return 0; /* sigma not > 0, or alloc fail in rng */
  }
  /* find unnormalized outer product of dof vectors */
  bzero(w, N*N*sizeof(*w));
  for (s = 0; s < dof; s++) {
    x1 = x + s*N; /* the s'th random vector */
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
	w[i*N+j] += x1[i] * x1[j];
  }
  /* done */
  free(x);
  return 1;
}

/*
 * Draw one n-by-n matrix "w" from Wishart distribution with 
 * "dof" degrees of freedom, parameterized by covariance matrix 
 * sigma (symmetric positive definite, n-by-n).
 * Returns 1 for OK, 0 for not.
 *   (not OK: bad sigma, or allocation failure)
 * Any integer dof > 0 is OK.
 * For convenience, multiplies the sample by a constant "scale"
 *   (Use scale = 1 for plain outer product matrix, which is the
 *    true Wishart r.v.  Use scale = 1/dof for a sample covariance 
 *    matrix that approaches sigma as dof grows.)
 */

int
rng_wishart(double* sigma, double dof, double scale, int N, double* w)
{
  int ok;
  int i;

  if (dof > N) {
    /* Wishart sample using Bartlett decomposition, good for large dof */
    ok = rng_wishart_chol(sigma, dof, N, w);
  } else {
    /* Wishart sample using naive outer products, good for small dof */
    ok = rng_wishart_outer(sigma, dof, N, w);
  }
  // rescale if needed
  if (ok && scale != 1.0)
    for (i = 0; i < N*N; i++)
      w[i] *= scale;

  return ok; /* 1 for OK, 0 for not */
}

