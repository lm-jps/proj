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
  }

