/*
 *  Bugs:
 *    Reads from and writes to segment #0 of input and output records: only
 *	safe if there is but one segment each in the input and ouput series
 *    Unconditionally reads and writes the data arrays as shorts, possibly
 *	invalid if BSCALE != 1 or BZERO != 0.
 *    Assumes that output record segments are 1024*1024, no checks.
 *    Assumes fixed path for crop lists in /home/soi/CM/tables/bin_list
 */

int get_vlist_crop (char *table_list, int *list, int *img_row_length) {
  DRMS_Array_t *table;
  HContainer_t *keylist;
  int status;
  short *tblvals;

  int addr, addr_lo, addr_hi, row_length, ct, n, pixct;

  keylist = NULL;
  if (!(table = drms_fitsrw_read (drms_env, table_list, 0, &keylist,
      &status))) {
    fprintf (stderr, "Error unable to read file %s as FITS\n", table_list);
    return -1;
  }
  tblvals = (short *)table->data;
  *img_row_length = tblvals[1];
  ct = tblvals[2];
  if (table->axis[0] < (3 * ct + 16)) {
						   /*  table list incomplete  */
    fprintf (stderr,"Error: pixel list %s is incomplete\n", table_list);
    return -1;
  }
	/*  whether table list contains extraneous entries seems unimportant  */
/*
  if (table->axis[0] > (3 * ct + 16)) {
    fprintf (stderr, "Warning: pixel list %s\n", table_list);
    fprintf (stderr,"         contains extraneous entries\n");
    fprintf (stderr, "%d vs %d (= 3*%d + 16)expected\n", table->axis[0],
	3 * ct + 16, ct);
  }
*/
  tblvals += 16;
  pixct = 0;
  while (ct--) {
    addr_lo = *tblvals++;
    addr_hi = *tblvals++;
    row_length = *tblvals++;
    if (addr_lo < 0) addr_lo += 65536;
    if (addr_hi < 0) addr_hi += 65536;
    if (row_length < 0) row_length += 65536;
    n = addr = addr_lo + 65536 * addr_hi;
    while (n <= addr + row_length) {
      list[pixct++] = n++;
    }
  }
  drms_free_array (table);

  return (pixct);
}

char *get_limbcrop_from_keylist (HContainer_t *keylist) {
  DRMS_Keyword_t *keyword = NULL;
  HIterator_t hit;
  unsigned long dpc;
  static char limbcrop[5];

  hiter_new_sort (&hit, keylist, drms_keyword_ranksort);
  while (keyword = hiter_getnext (&hit)) {
    if (!strcmp (keyword->info->name, "DPC")) {
			   /*  assume keyword->info->type = DRMS_TYPE_STRING  */
      dpc = strtol (keyword->value.string_val, NULL, 16);
      sprintf (limbcrop, "%04x", (int)(dpc & 0xffff));
      hiter_free (&hit);
      return limbcrop;
    }
  }
  strcpy (limbcrop, "0000");
  hiter_free (&hit);
  return limbcrop;
}

char *get_limbcrop_from_keyword (DRMS_Record_t *irec) {
  static char limbcrop[5];

  int status;

  char *dpcstr = drms_getkey_string (irec, "DPC", &status);
  unsigned long dpc = strtol (dpcstr, NULL, 16);

  if (status) strcpy (limbcrop, "0000");
  else {
    dpc = strtol (dpcstr, NULL, 16);
    sprintf (limbcrop, "%04x", (int)(dpc & 0xffff));
  }
  return limbcrop;
}

int embed_limbpixels (const char *ifile, DRMS_Record_t *orec, float *cropmin,
    float *cropmax) {
  static DRMS_Array_t *img, *limbdat;
  static float cropmn, cropmx;
  static int *vlist;
  static int img_vals, ntot, initial = 1;
  static unsigned short binparam, limbparam;
  static short *blank, *recon;
  static char lastcrop[] = {"0000"};
  static char filedir[] = {"/home/soi/CM/tables/bin_list"};

  DRMS_Segment_t *iseg, *oseg;
  HContainer_t *keylist;
  int axes[2];
  int cols, rows, n;
  short *orig;
  char *limbcrop;
  char bin_list[DRMS_MAXPATHLEN];

  int status = 0;

  if (initial) {
    axes[1] = rows = axes[0] = cols = 1024;
    ntot = cols * rows;
    vlist = (int *)malloc (ntot * sizeof (int));
    blank = (short *)malloc (ntot * sizeof (short));
    recon = (short *)malloc (ntot * sizeof (short));
    for (n = 0; n < ntot; n++) blank[n] = -32768;
						  /*  create the image array  */
    img = drms_array_create (DRMS_TYPE_SHORT, 2, axes, (void *)recon, &status);
    initial = 0;
  }
				/*  assume only one segment in output series  */
  oseg = drms_segment_lookupnum (orec, 0);
  keylist = NULL;
  if (!(limbdat = drms_fitsrw_read (drms_env, ifile, 0, &keylist,
      &status))) {
    fprintf (stderr, "Error unable to read file %s as FITS \n", ifile);
    return status;
  }
  orig = (short *)limbdat->data;
  limbcrop = get_limbcrop_from_keylist (keylist);
  if (strcmp (limbcrop, lastcrop)) {
				   /*  a new limb crop DPC: update the vlist  */
    sprintf (bin_list, "%s/limb%s.fits", filedir, limbcrop);
    img_vals = get_vlist_crop (bin_list, vlist, &cols);
    sscanf (limbcrop, "%04x", &limbparam);
    binparam = limbparam & 0x0fff;
    cropmn = 0.125 * binparam;
    cropmx = cropmn + 0.5 * (limbparam >> 12);
    strcpy (lastcrop, limbcrop);
  }
						   /*  blank the image array  */
  recon = (short *)(img->data);
  memcpy (recon, blank, ntot * sizeof (short));
  for (n = 0; n < img_vals; n++) recon[vlist[n]] = orig[n];
						     /*  write output record  */
  if (drms_segment_write (oseg, img, 0)) {
    fprintf (stderr, "Warning: unable to write to output segment; skipped\n");
    return 1;
  }

  *cropmin = cropmn;
  *cropmax = cropmx;

  return status;
}

