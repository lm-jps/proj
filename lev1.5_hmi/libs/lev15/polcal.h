struct polcal_struct {
  int method;
  int nin;
  double *xin,*yin;
  double tsela,tfronta;
  double *fqq_0,*fqq_1,*fqq_2,*fqu_0,*fqv_0,*fuu_0,*fuu_1,*fuu_2,*fuv_0,*fvv_0,*fvv_1,*fvv_2;
  double *ret1_0,*ret1_1,*ret2_0,*ret2_1,*ret3_0,*ret3_1,*phi1_0,*phi2_0,*phi3_0;
};

int init_polcal(struct polcal_struct *pars, int method, const char *paramFile);

int free_polcal(
  struct polcal_struct *pars
);

int polcal(
  struct polcal_struct *pars,
  int nframe, // Number of input polarization states
  int mode, // 1 for IQUV, 2 for LCP+RCP from full, 3 for LCP+RCP from 2 pol
  float **input, // Pointers to input filtergrams
  float **output, // Pointers to output polarization images
  int *ps1, // PS1 positions
  int *ps2, // PS2 positions
  int *ps3, // PS3 positions
  float tsel, // Polarization selector temperature
  float tfront, // Front window temperatures
  int nx, // Number of colums
  int ny, // Number of rows
  int nlead // Declared number of rows. nlead>=nx
);

char *polcal_version(); // Returns CVS version of polcal.c
