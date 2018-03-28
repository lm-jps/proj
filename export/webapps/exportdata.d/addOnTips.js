// Place for help tips for add-on code.
// first param in Tip is HTML id string for element where tip should be linked.
// second param is tip message,
// final param is tipstyle.
// This function is called from exportdata.html.

function AddOnTips(tipstyle, Tip)
  {
  new Tip("MaprojHelp","Extract patch in heliographic coordinates, trac across the disk, " +
      "provide Carrington location and size of patch in degrees or ref time and location.",tipstyle);
  new Tip("MaprojXHelp","Enter the Carrington longitude (degrees) for the center of the desired patch cutout.  If the disk center " +
      "box is checked, the longitude at disk center at the time of the first image processed will be used.",tipstyle);
  new Tip("MaprojYHelp","Enter the Carrington latitude (degrees) for the center of the desired patch cutout.",tipstyle);
  new Tip("MaprojProjHelp","Select desired projection from the drop-down list.",tipstyle);
  new Tip("MaprojScaleHelp","Select target plate scale in degrees per pixel. Note should be <= given limit or undersampling may occur. " +
      "This can be pushed up a factor of 1.5 or 2 if you are not planning on making time a sequence of images.",tipstyle);
  new Tip("MaprojWideHelp","Select target image width in pixels of width defined by scale above.",tipstyle);
  new Tip("MaprojHighHelp","Select target image height in pixels of height defined by scale above.",tipstyle);
  new Tip("MaprojVerifyHelp","Click to run internal consistency check on selections above.",tipstyle);
  new Tip("MaprojNoaaTip","Enter 4- or 5-digit NOAA AR number to preload location type, reference time, " +
      "and location for the NOAA AR region.  Will find entry in NOAA region table that is closest to CM.",tipstyle);
  new Tip("ImPatchHelp","Extract patch in heliographic coordinates, trac across the disk, " +
      "provide Carrington location and size of patch in degrees or ref time and location.",tipstyle);
  new Tip("ImNoaaTip","Enter 4- or 5-digit NOAA AR number to preload location type, reference time, " +
      "and location for the NOAA AR region.  Will find entry in NOAA region table that is closest to CM.",tipstyle);
  new Tip("ImPatchResetHelp","After change in series use these to reset all impatch params or just the T_START and T_STOP times.", tipstyle);
  new Tip("ImTStartHelp", "T_START, T_STOP, and CADENCE take default values from the RecordSet specification. " +
      "The final query is formed from T_START and T_STOP.", tipstyle);
  new Tip("ImTStopHelp", "T_START, T_STOP, and CADENCE take default values from the RecordSet specification. " +
      "The final query is formed from T_START and T_STOP.", tipstyle);
  new Tip("ImCadenceHelp", "Cadence for image patch extraction.  Taken from RecordSet spec if present.  Value here " +
      "will override the RecordSet spec if different.  Use <number><s|m|h|d> format.", tipstyle);
  new Tip("ImLocUnitsHelp", "The patch center locationa(x,y) can be specified in Carrington Rot, Lat, and Longitude or with " +
      "a reference time and location in arcsec, or degrees from CM and equator or in pixels from the lower left of " +
      "the image (1,1) before any implied rotation or center offset.", tipstyle);
  new Tip("ImTRefHelp", "One of T_REF or CAR_ROT is required to locate the patch in time. " +
      "If T_REF is absent, T_START will be used for T_REF.  Data must exist within 4 hours of T_REF. ", tipstyle);
  new Tip("ImCARROTHelp", "One of T_REF or CAR_ROT is required to locate the patch in time. " +
      "If CAR_ROT is absent for Carrington lat-lon locunits then it will be computed as the " +
      "rotation with CM at T_REF. ", tipstyle);
  new Tip("ImXHelp", "X will be in arcsec, degrees, or pixels from lower-left corner as (1,1) depending on LocUnits.", tipstyle);
  new Tip("ImYHelp", "Y will be in arcsec, degrees, or pixels from lower-left corner as (1,1) depending on LocUnits.", tipstyle);
  new Tip("ImBoxUnitsHelp", "the image patch extract size can be specified in pixels, arcsecs, or degrees.  In the case of " +
      "degrees, the pixel size is computed for the time that the patch is centered at central meridian.", tipstyle);
  new Tip("ImWideHelp", "Width of extract patch.  pixels, or arcsecs converted to pixels with CDELT1, or degrees converted " +
      " to pixels with img2sphere projection based on DSUN_OBS, RSUN_REF, T_Ref, X, and Y.and CDELT1.", tipstyle);
  new Tip("ImHighHelp", "Height of extract patch.  pixels, or arcsecs converted to pixels with CDELT1, or degrees converted " +
      " to pixels with img2sphere projection based on DSUN_OBS, RSUN_REF, T_Ref, X, and Y.and CDELT1.", tipstyle);
  new Tip("ImVerifyHelp", "Check all parameters for presence if required and some checks for consistency.", tipstyle);
  new Tip("ImagesHelp","Select details for images and movies.",tipstyle);
  new Tip("HmiB2ptrHelp","HmiB2ptr processing reprojects HMI Vector field B to spherical coordinates (phi,theta,r).",tipstyle);
  new Tip("AiaScaleHelp",'When normalize scale processing is selected, the images will be scaled to AIA normal scale ' +
      'of 0.6 arcsec per pixel.',tipstyle);
  new Tip("AiaScaleMptHelp",'Not sure what to say here ', tipstyle);
  new Tip("AiaScaleCutoutXcHelp",'Not sure what to say here ', tipstyle);
  new Tip("AiaScaleCutoutYcHelp",'Not sure what to say here ', tipstyle);
  new Tip("AiaScaleCutoutWideHelp",'Not sure what to say here ', tipstyle);
  new Tip("AiaScaleCutoutHighHelp",'Not sure what to say here ', tipstyle);
  new Tip("RebinHelp",'Reduce or Enlarge image size.  If scale < 1 reduce size by scale, reciprocal of scale ' +
      'rounded to integer.  Reduction can be either by boxcar binning or Gaussian smoothing to specified FWHM(pixels) ' +
      'averaged over a vector of specified length. If size enlarged, scale will be rounded to integer and pixel values ' +
      ' are replicated.  ',tipstyle);
  new Tip("ResizeHelp",'Perform sub-pixel registration with optional rescaling. ' +
      'Each image will be co-registered to either Sun center moved to center of array or to center offset of the first ' +
      'image present.  Setting scale to 0.6 is same as aia_scale processing for aia.lev1p5.  ' +
      'The interpolation method is same as used in ANA and IDL for co-registration.',tipstyle);
  }

