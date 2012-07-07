// Place for help tips for add-on code.
// first param in Tip is HTML id string for element where tip should be linked.
// second param is tip message,
// final param is tipstyle.
// This function is called from exportdata.html.

function AddOnTips(tipstyle, Tip)
  {
  new Tip("HgPatchHelp","Extract patch in heliographic coordinates, trac across the disk, " +
      "provide Carrington location and size of patch in degrees or ref time and location.",tipstyle);
  new Tip("HgNoaaTip","Enter 4- or 5-digit NOAA AR number to preload location type, reference time, " +
      "and location for the NOAA AR region.  Will find entry in NOAA region table that is closest to CM.",tipstyle);
  new Tip("ImagesHelp","Select details for images and movies.",tipstyle);
  new Tip("AiaScaleHelp","When normalize scale processing is selected, the images will be scaled to AIA normal scale of 0.6 arcsec per pixel.",tipstyle);
  new Tip("RebinHelp","Reduce or Enlarge image size.  If scale < 1 reduce size by scale, reciprocal of scale rounded to integer.  Reduction can be either by boxcar binning or Gaussian smoothing to specified FWHM(pixels) averaged over a vector of specified length. If size enlarged, scale will be rounded to integer and pixel values are replicated.  ",tipstyle);
  new Tip("ResizeHelp","Perform sub-pixel registration with optional rescaling.  Each image will be co-registered to either Sun center moved to center of array or to center offset of the first image present.  Setting scale to 0.6 is same as aia_scale processing for aia.lev1p5.  The interpolation method is same as used in ANA and IDL for co-registration.",tipstyle);
  }

