//#include <fresize.h>

//char *X0_MP_key="X0_MP";
//char *Y0_MP_key="Y0_MP";
//char *RSUN_OBS_key="RSUN_OBS";
//char *IMSCL_MP_key="IMSCL_MP";

//char *HCAMID_key="HCAMID";
//char *HCFTID_key="HCFTID";

//const int light_val1=2;
//const int light_val2=3;

  const int malign=32;
  const float fwhm=1.0; //full width half maximum for gaussian filter 
  const float sigmamin=0.0;  //lower limit for standard deviation
  const float sigmamax=700.0;  //upper limit for standard deviation
  const float kernel_size=4.0;

const int cent_frac=8; //center portion of image that is used to calculate std of image
const float limit=10.5; //detection limit (limit * sigma)
const float maxval=119815.0; //90 % of max value
