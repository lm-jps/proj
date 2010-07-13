
/*
 * THIS IS A GENERATED FILE; DO NOT EDIT.
 *
 * Declarations for calling `mrf_segment_wts' a.k.a. `msw' as a C library.
 *
 * Made by intermediate binary `mrf_segment_wts.out' on:
 * 	Mon Jun  7 15:37:56 2010
 *
 * Code for include-generation driver `../Gen-include.c' last modified on:
 * 	Mon Jun  7 15:11:28 2010
 *
 */

// Original documentation block
/*
 * mrf_segment_wts: segment with a discrete Markov random field
 * 
 *  [yp,post] = mrf_segment_wts(iter,T,beta,alpha,dist,y,lprob1,...,lprobK)
 *  * Performs iter full sweeps (using temperature schedule T) of
 *  Gibbs sampling on an input labeling y, with entries in 1..K, to produce
 *  an output labeling yp.  Posterior probability is optionally returned.
 *  Posterior is actually an energy function, which is correct with respect
 *  to changes in y, but does not have the correct scale factors which
 *  vary with beta, alpha, and dist.
 *  * Conditional distributions of pixel data x given class y are
 *  calculated externally, and given by lprob1...lprobK, as log-
 *  probabilities.
 *  * Per-class biases are given by alpha, which can be empty, indicating
 *  no bias.  Otherwise, the interpretation is alpha(k) is the prior
 *  log-probability of seeing class k.
 *  * The smoothness parameter of the MRF (Potts model) is beta.  If beta
 *  is a matrix, beta(k,l) is the smoothness "reward" given to a site of
 *  class k for having a neighbor of class l.  The scalar beta thus
 *  corresponds to a diagonal matrix with repeated beta entries.
 *  (beta and alpha agree with definitions in the Besag paper below.)
 *  * Pixel-pixel distances are given by dist, where dist(nu,m,n) gives
 *  the distance between pixel (i,j) and its neighbor number nu, looking
 *  up or left.  For the 3x3 neighborhood, pixel s=(m,n) has 8 neighbors s',
 *  and distances to 4 of them, where s'<s, at offsets:
 *    (-1,-1),(0,-1),(1,-1),(-1,0),
 *  are given in that order, in dist(i,j,:).  The other neighbors
 *  have s'>s, and the corresponding distances are listed in the
 *  symmetric entries of dist.  If dist=[], it is taken to be
 *  everywhere 1, thus removing the direction-sensitive smoothing
 *  (but still smoothing).
 *  * Classes are 1..K, but labels NaN and 0 are not updated or counted
 *  as neighbors.  Any NaN in a log-probability forces a NaN in the
 *  output class.
 *  * A `clock' at each pixel may speed computation.  This recognizes that,
 *  if the neighbors of a pixel do not change, the Gibbs sampler is rolling
 *  a stationary die once per iteration to determine the pixel label.
 *  The waiting time until another label change is then geometric,
 *  and this waiting time can be sampled once to short-circuit a series
 *  of die rolls.  See Ripley, below.
 *  * The speedup of this method is greatest when few labels change.  If
 *  supplied, iter(2) is the threshold (in [0,1]) of #changed/#labels
 *  for a switchover from the ordinary method to the clock method.
 *  (Only one switch is permitted in the iteration sequence.)
 *  If iter(2) = 0, the ordinary coin-flip method is used throughout;
 *  if =1, the clock method is used throughout.
 *  * Annealing is permitted through the `temperature' parameter T.
 *  Initial temperature is T(1), reduced by a factor of T(2) at each
 *  iteration. The default T(2) = 1 suppresses annealing.
 *  Clocks are reset whenever the cumulative drop in temperature
 *  reaches a factor of T(3).  (Roughly, then, clocks are reset every
 *  log(T(3))/log(T(2)) iterations.)  If T(4) is present and equals
 *  zero, additional iterations are done at zero temperature, at the
 *  end of the normal annealing schedule, until the labels reach a
 *  fixed point.  (This is equivalent to Besag's ICM.)
 *  * If iter[1] has a fractional part, that part is multiplied by 2**31
 *  and rounded to set the desired random number seed; if not, a
 *  pseudorandom seed is generated.  This is provided to allow repeatability.
 *  * This is implemented as a MEX file.
 * 
 *  Inputs:
 *    int iter[1] or [2] = [0 0.05];
 *    real T[0] or [1] or [2] or [3] or [4] = [1 1 0.8 0];
 *    real beta[1] or [K,K];
 *    real alpha[K] or [0] = [];
 *    real dist[nbr,m,n] or [0];
 *    int y[m,n];  -- 1 <= y <= K, or 0 or NaN
 *    real lprob1[m,n];
 *    ...
 *    real lprobK[m,n];
 * 
 *  Outputs:
 *    int yp[m,n];  -- 1 <= y <= K, or 0 or NaN
 *    opt real post;
 * 
 *  See Also:  mrf_segment, makemrfdiscwts
 *  Geman and Geman, Stochastic relaxation, Gibbs distributions, and
 *  the Bayesian restoration of images, IEEE PAMI Nov. 1984
 *  J. Besag, On the statistical analysis of dirty pictures, JRSSB, 1986
 *  B. D. Ripley, Statistical Inference for Spatial Processes,
 *  Cambridge U., 1988, p. 99
 * 
 *  turmon sep/oct 2006, weighted distance
 *  turmon 5 march 1999, modified to streamline for batch operations
 *  MJT 18 Mar 1996
 * 
 * 
 */

#ifndef _mexfn_mrf_segment_wts_h_
#define _mexfn_mrf_segment_wts_h_

// function entry point
mexfn_lib_t main_mrf_segment_wts;

// argument counts
#define MXT_msw_NARGIN_MIN 	7
#define MXT_msw_NARGIN_MAX 	14
#define MXT_msw_NARGOUT_MIN	0
#define MXT_msw_NARGOUT_MAX	2

// input argument numbers
#define MXT_msw_ARG_iters	0
#define MXT_msw_ARG_T	1
#define MXT_msw_ARG_beta	2
#define MXT_msw_ARG_alpha	3
#define MXT_msw_ARG_dist	4
#define MXT_msw_ARG_y	5
#define MXT_msw_ARG_lprob1	6
#define MXT_msw_ARG_lprob2	7
#define MXT_msw_ARG_lprob3	8
#define MXT_msw_ARG_lprob4	9
#define MXT_msw_ARG_lprob5	10
#define MXT_msw_ARG_lprob6	11
#define MXT_msw_ARG_lprob7	12
#define MXT_msw_ARG_lprob8	13

// output argument numbers
#define MXT_msw_ARG_yp	0
#define MXT_msw_ARG_post	1


#endif // _mexfn_mrf_segment_wts_h_

// (file ends)
