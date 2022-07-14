struct fill_hash_struct {
  int nbad; // Number of bad points for this entry
  int nhit; // Number of times used
  int *wbad; // List of bad pixels for this entry
  double *coeff; // Saved coefficients
  float cnorm;
  float ierror;
  struct fill_hash_struct *next; // Next item in list
};

struct fill_struct {
  int order; // Interpolation order
  double *a0,*a00,*rh0,*a0t,*rh0t,acort00;
  double *a,*rh,*a1b,*a1r;
  int *wgood;
  int *wbad;
  int hashmod; // Modules for hash calculation
  int *hashcount; // How many of this hash value
  int ndiff; // Number of different masks seen
  int ncollision; // Number of times with same hash and nbad
  struct fill_hash_struct **hashtable;
};

  int init_fill(
  int method, // Interpolation method
  double pnoise, // Level of photon noise for trade-off
  int order, // Interpolation order. Generally odd.
  int targetx, // Target point in x (normally (order-1)/2)
  int targety, // Target point in y (normally (order-1)/2)
  struct fill_struct *pars // Structure to save setup information etc.
);

int free_fill(
  struct fill_struct *pars
);

int fgap_fill(
  struct fill_struct *pars, // Parameters from call of init_fill
  float *im, // Input image
  int nx, // Size of array in dimension adjacent in memory
  int ny, // Size of array in dimension not adjacent in memory
  int nlead, // Leading dimension of array nlead>=nx
  unsigned char *mask, // Mask is 0 data point is there
             // Mask is 1 if data point is not there and should be filled
             // Mask is 2 if data point is not there and should not be filled
  float *cnorm, // White noise magnification
  float *ierror // Estimated interpolation error
);
