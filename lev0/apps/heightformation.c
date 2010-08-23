//CORRECTION OF HEIGHT FORMATION (from Sebastien 10Aug2010)
//returns 0 if corrections were successful, 1 otherwise
int heightformation(int FID, double OBSVR, float *CDELT1, float *RSUN, float *CRPIX1, float *CRPIX2, float CROTA2)
{
  int wl=0;
  int status=0;
  float correction=0.0,correction2=0.0;
 
  wl = (FID/10)%20;  //temp is now the filter index
 
  if( (wl >= 0) && (wl < 20) )
    {
      correction  = 0.445*exp(-(wl-10.-(float)OBSVR/(0.690/6173.*3.e8/20.)-0.25)*(wl-10.-(float)OBSVR/(0.690/6173.*3.e8/20.)-0.25)/7.1);
      correction2 = 0.39*(-2.0*(wl-10.- (float)OBSVR/(0.690/6173.*3.e8/20.)-0.35)/6.15)*exp(-(wl-10.-(float)OBSVR/(0.690/6173.*3.e8/20.)-0.35)*(wl-10.-(float)OBSVR/(0.690/6173.*3.e8/20.)-0.35)/6.15);
     
      *CDELT1 = *CDELT1*(*RSUN)/((*RSUN)-correction);
      *RSUN   = *RSUN-correction;
      *CRPIX1 = *CRPIX1-cos(M_PI-CROTA2*M_PI/180.)*correction2;
      *CRPIX2 = *CRPIX2-sin(M_PI-CROTA2*M_PI/180.)*correction2;
    }
  else status=1;
 
  return status;
}
