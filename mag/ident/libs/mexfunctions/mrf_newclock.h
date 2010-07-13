/***************************************************************

turmon nov 2006

Routines for generating a geometrically-distributed clock,
where we do not care if the clock is larger than M.
That is, if the return value is greater than M, the routine
can choose to not compute it exactly, and return any number
larger than M.

tested new_clock in all three regimes: low q, low p, and moderate
p, corresponding to its three cases.  
In the low-q case, I tried q=0.7, 0.5, 0.1, and 0.01 with good 
agreement between method 0 and method 2, versus the true density.
In the low p case, I used 
q = 1-p = 0.999, M=17, n=1e8, so M*p=0.017 < 1.75e-2, while giving
good contrast between Pr(clk == 1) and Pr(clk == M) as estimated
over n samples (ie, 10.0e4 versus 9.84e4 counts, a difference of 
1600 counts in a quantity with standard deviation about 314 counts.
I also did a rough check at q=0.999999, but there is no appreciable
variation of the density across clock value there, so I was not able
to validate the shape of the PDF.
In the moderate-p case, I used q=0.9.  This code is identical for
all "mode" inputs, so there's not much to check.

****************************************************************/

/*
 * this code needs:
 *   typedef for Clock, currently unsigned
 *   maximum value CLOCK_MAX (for type Clock), currently UINT_MAX
 */

/*
 * get a new clock number based on a geometric(1-q) random variable
 *
 * this is based on the dead-simple inversion method for the
 * exponential random variable, and then quantizing it
 *
 */
static
Clock
new_clock_easy(double q, int M)
{
  double rng_uniform();
  double clock; /* need the range */

  clock = log(rng_uniform())/log(q) + 1; /* initial choice */
  /* fix it up if needed */
  if (clock > CLOCK_MAX)
    clock = CLOCK_MAX; 
  return((Clock) floor(clock));
}

/*
 * get a new clock number based on a geometric(1-q) random variable
 *
 * this is the older version of the new_clock routine, which I 
 * replaced in 2006.  It special-cases the low and high q values
 * with deterministic values.  Because this is a crude way to do
 * the approximation, only very extreme values are special-cased.
 */
static
Clock
new_clock_orig(double q, int M)
{
  double rng_uniform();
  double clock; /* need the range */

  /* Cheat a little here by handling interval edges, which occur often,
   * as special cases, without drawing a new random variable.
   * This is somewhat speedier.  Size of edges is 10^-8 on each
   * side.
   * The deterministic approximation only fails on either endpoint
   * with probability about 10^-8, so even if all pixels are approximated,
   * only about 2 pixels in a 10K x 10K image will be in error.
   */
  if (q < 1e-8)
    clock = 1.0;
  else if (q > (1 - 1e-8)) /* -> log(p) < 1e-8 in magnitude by Taylor */
    clock = CLOCK_MAX;
  else {
    clock = log(rng_uniform())/log(q) + 1; /* initial choice */
    /* fix it up if needed */
    if (clock > CLOCK_MAX)
      clock = CLOCK_MAX; 
    /* NB, this approximation is compensated for by the "staleness"
     * checks within the main loop, that expire all clocks before 
     * CLOCK_MAX iterations.  There is no loss of precision from
     * this source.
     */
  }
  return((Clock) floor(clock));
}

/*
 * Get a new draw X of a random clock 
 * This is a geometric(p) random variable, where p = 1-q.
 * That is, the hold probability is q, and the change probability is p,
 * so X is formed by sequences like:
 *
 * X = 1: C    Pr = p
 * X = 2: HC   Pr = q p
 * X = 3: HHC  Pr = q^2 p
 * X = 4: HHHC Pr = q^3 p
 * etc.
 * 
 * In general, Pr(X=k) = p q^k for k=1,...
 * And, P(X>k) = q^k, because all such events start out with k "H" letters
 * 
 * We special-case the low and high q values to avoid the rather expensive
 * double log() required by the inversion method.  The important special
 * case is for q near 1, because it is most common.  This case requires
 * some careful reasoning.  The less-common special case is for small
 * q, where we just mean q less than 0.5 or so.  In this case the
 * geometric RV is easily simulated by implementing the Hold/Change
 * sequences illustrated above.
 */
static
Clock
new_clock(double q, int M)
{
  double rng_uniform();
  const double p = 1-q;
  double x; /* the random variable: need the range */

  if (q < 0.75) {
    /* CASE 1.  q small, or hold time is low.
     * Draw by the iteration method.
     * The loop below accumulates events of probability q.
     */
    /* The constant in the test above affects run time (by
     * selecting the fastest generation method) but not accuracy.  
     * At q = 0.75, on average we draw 1/(1-q) = 4 uniforms, 
     * about equal to the one uniform plus two logs of the inversion
     * method below in CASE 3 (G5, 11/2006).  At q=0.2, the code
     * below is about 2.2x faster than inversion.  At q=0.01, it
     * is 3x faster.  But note, this is not the dominant case.  */
    for (x = 1; rng_uniform() < q && x<=M; x++);
  } else if (M*p < 1.75e-2) {
    /* CASE 2.  q near 1, or p near 0.
     * Hold time is high relative to number of iterations
     * remaining; typically will return a value X > M.
     * This is the most frequent case in our application.
     * The code below is a factor of 2.5x faster than the default
     * inversion method on a G5.
     */
    /* Below, we take one of two branches with probability
     * Pr(X>M) = pow(1-p,M) versus Pr(X<=M).
     * In the former, overwhelmingly probable case, just return
     * CLOCK_MAX.  In the latter, generate X conditioned on X<=M.
     * 
     * Because we cannot quickly compute pow(1-p,M), we approximate it
     * with the first two terms of its binomial expansion:
     *   P(X>M) = (1-p)^M ~= 1 - M * p + M(M-1)/2 * p^2
     * For M*p < 1.75e-2, the relative accuracy 
     *   A(M,p) = [ 1 - M * p + M(M-1)/2 * p^2 ] / [ (1-p)^M ]
     * is better than 1 part in 10^6.  For M*p < 3.5e-3, the accuracy
     * is better than 1 part in 10^8.
     * The test below replaces 
     *   rng_uniform() < 1-rho 
     * with its equivalent,
     *   rng_uniform() > rho
     * where rho = M * p - M(M-1)/2 * p^2
     *           = M*p * (1 - 0.5 * (M-1) * p)
     *           = M*p * (1 - 0.5 * (M * p - p) )
     */
    if (rng_uniform() > M*p*(1 - 0.5*(M*p-p))) {
      /* CASE 2a.  Generate X | X>M, which is deterministic. */
      /* This branch is taken with probability approximately P(X>M),
       * and returns X conditioned on X > M.  We almost always
       * take this branch since M*p << 1.  */
      x = CLOCK_MAX;
    } else {
      /* CASE 2b.  Generate X | X<=M */
      /* This branch is taken with probability approximately P(X<=M),
       * and returns X conditioned on X <= M.  Censoring the 
       * realizations with X>M would be one way.  But it's easier to 
       * note that, if X is geometric, then (X-1) % M + 1 is geometric
       * conditioned on X<=M.  We produce X-1 by taking floor
       * and not ceil in the log-of-uniform below.
       * We almost never take this branch, so its efficiency is moot. */
      x = fmod(floor(log(rng_uniform())/log(q)), (double) M) + 1;
    }
  } else {
    /* CASE 3. Hold probability is in between: draw by inversion method. */
    x = floor(log(rng_uniform())/log(q)) + 1; /* initial choice */
    /* don't overflow the clock size; note, CLOCK_MAX >= M anyway */
    if (x > CLOCK_MAX)
      x = CLOCK_MAX; 
  }
  // mexPrintf("\tp = %g\tclock=%u\n", p, (Clock) x);
  // note, x is already an integer > 0
  return((Clock) x);
}

