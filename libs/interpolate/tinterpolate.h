#define tavg_boxcar 1
#define tavg_cosine 2
#define tavg_fourth 3
#define tavg_hathaway 4

int tinterpolate(
  int nsample, // Number of input times
  double *tsample, // Input times
  double tint, // Target time
  int nconst, // Number of polynomial terms exactly reproduced
  float **images, // Pointer array to input images
  unsigned char **masks, // Pointer array to input masks. 0=good, 1= missing
  float *image_out, // Interpolated image
  int nx, // Number of points in dimension adjacent in memory
  int ny, // Number of points in dimension not adjacent in memory
  int nlead, // Leading dimension of arrays. nlead>=nx
  int method, // Interpolation method
  char **filenamep // Pointer to name of file to read covariance from.
                   // Set to actual file used if method > 0.
);

int taverage(
  int nsample, // Number of input times
  double *tsample, // Input times
  double tint, // Target time
  int nconst, // Number of polynomial terms exactly reproduced
  float **images, // Pointer array to input images
  unsigned char **masks, // Pointer array to input masks. 0=good, 1= missing
  float *image_out, // Interpolated image
  int nx, // Number of points in dimension adjacent in memory
  int ny, // Number of points in dimension not adjacent in memory
  int nlead, // Leading dimension of arrays. nlead>=nx
  int method, // Interpolation method
  char **filenamep, // Pointer to name of file to read covariance from.
                   // Set to actual file used if method > 0.
  int avmethod, // averaging method
  int order, // Interpolation order
  double cadence, // Cadence to interpolate to
  int hwidth, // In units of cadence. Total width is 2*hwidth+1
  double par1 // In units of cadence.
);

