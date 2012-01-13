/*
 * mharp_match.c
 *
 * Utilities for matching active regions (HARP <-> NOAA association)
 *
 * reads matching AR info from su_rsb.NOAA_ActiveRegions
 * puts the info into our own data structures
 * develops correspondence between NOAA regions and HARPs
 */


/*
 * Summarize a match ar by printing it -- diagnostic aid
 *
 * (s is a prefix string)
 */
static
void
match_ar_print(FILE *fp, char *s, match_info_t *m)
{
  const int human_friendly = 0;
  char trec_str[24];

  if (human_friendly) {
    if (m->valid != Patch_Normal) {
      fprintf(fp, "%sMatch-Patch not valid\n", s);
    } else {
      sprint_time(trec_str, m->t, "TAI", 0);
      fprintf(fp, 
	      "%sID = %d (%s) @ (%.2f,%.2f) A = %.1f\t[%s/%s]\n", 
	      s, m->id, trec_str, m->lat, m->lon, m->area, m->class_mag, m->class_zur);
    }
  } else {
    // machine-friendly output format
    // (time is #seconds since epoch)
    fprintf(fp, 
	    "%d\t%d\t%.1f\t%.2f\t%.2f\t%.1f\t%.1f\n", 
	    m->id, m->valid, (double) m->t, m->lat, m->lon, m->lon_wid, m->area);
  }
}


/*
 * interpolate match AR info at time t, by filling in *m 
 * We assume that time t is between m0 (before) and m1 (after)
 */
static
void
match_ar_interpolate(TIME t, match_info_t *m0, match_info_t *m1, match_info_t *m)
{
  double dt0 = t - m0->t; // > 0
  double dt1 = m1->t - t; // > 0
  double wt0, wt1;

  // check for perfect matches at either endpoint 
  //  (weight computation would yield 0/0)
  if (m0 == m1) {
    wt0 = 1; // m0 == m1: code for a close-enough perfect match to both m0 and m1
  } else if (dt0 == 0) {
    wt0 = 1; 
  } else if (dt1 == 0) {
    wt0 = 0;
  } else {
    // weights for imperfect matches
    wt0 = dt1 / (dt0 + dt1);
  }
  wt1 = 1.0 - wt0;
  // nearest-neighbor interpolation by default (discrete quantities)
  if (wt0 >= 0.5) {
    memcpy(m, m0, sizeof(*m));
  } else {
    memcpy(m, m1, sizeof(*m));
  }
  // linear interpolation overwrites the continuous quantities
  m->lat      = wt0 * m0->lat + wt1 * m1->lat;
  m->lon      = wt0 * m0->lon + wt1 * m1->lon;
  m->lon_carr = wt0 * m0->lon_carr + wt1 * m1->lon_carr;
  m->lon_wid  = wt0 * m0->lon_wid + wt1 * m1->lon_wid;
  m->area     = wt0 * m0->area + wt1 * m1->area;
  // set time explicitly because if m0 == m1, t can be != m0->t 
  // (i.e. we allow inexact matches)
  m->t = t;
}

/*
 * Generate a sample from the raw match-ar stream MI_raw,
 * for AR numbered `id', at time `t'.
 * This routine just finds the bracketing times, before and after, and 
 * calls the interpolator within that bracket.
 * Returns nonzero status if there no way to interpolate at time t,
 * which is not an error condition in our setup.
 */
static
int
match_ar_resample(int id, TIME t, match_info_t *MI_raw, int Nraw, match_info_t *MI)
{
  int i, iLO, iHI;
  double dt, dt0, dt1;
  const double TOO_FAR = 3600*24*5;  // N days away is too far for interpolation
  const double CLOSE_ENOUGH = 3600*1;// ok to extrapolate N hours away

  // loop over all trec's, finding bracketing times
  iLO = iHI = -1; // nonsense flag value
  dt0 = -HUGE_VAL;
  dt1 =  HUGE_VAL;
  for (i = 0; i < Nraw; i++) {
    if (MI_raw[i].id != id) continue;
    dt = MI_raw[i].t - t;
    if ((dt < -TOO_FAR) || (dt > TOO_FAR))
      continue;
    // "earlier than t" bracket
    if (dt <= 0) {
      // if this AR was before t...
      if (dt > dt0) {
	// ...and, the gap is better than the best existing gap,
	// then make this one the new best
	dt0 = dt; // dt0 <= 0
	iLO = i;
      }
    } 
    // "later than t" bracket
    if (dt >= 0) {
      // if this AR was after t...
      if (dt < dt1) {
	// ...and, the gap is better than the best existing gap,
	// then make this one the new best
	dt1 = dt;  // dt1 >= 0
	iHI = i;
      }
    } 
  }
  // allow a match even if t is not strictly in [t_first, t_last]
  // (a) no later bracket found, and t is just earlier than the early bracket
  if ((iHI == -1) && (iLO >= 0 && fabs(dt0) < CLOSE_ENOUGH))
    iHI = iLO;
  // (b) no earlier bracket found, and t is just later than the late bracket
  if ((iLO == -1) && (iHI >= 0 && fabs(dt1) < CLOSE_ENOUGH))
    iLO = iHI;
  // this condition means no bracketing times were found
  // note: this is not a hard error; the AR could be inactive at
  // this time
  if (iLO < 0 || iHI < 0)
    return 1;
  // OK to interpolate
  match_ar_interpolate(t, MI_raw+iLO, MI_raw+iHI, MI);
  return 0;
}

/*
 * Get info for a single match AR by reading keywords.
 * returns nonzero for trouble
 */
static
int
match_get_single_info(int *id, match_info_t *mI, DRMS_Record_t *inRec)
{
  char *str;
  int s1 = 0; // status from single command
  int s  = 0; // overall status, OR of all s1's

  *id = drms_getkey_int(inRec, "RegionNumber", &s1); s |= s1;
  if (mI == NULL)
    return s;
  // get keys in order of structure
  // note, LongitudeCM is the location referenced to the current central meridian
  // LongitudeHG is the carrington longitude
  // LongitudeCM changes from day to day, LongitudeHG changes little
  mI->valid    = Patch_Normal;
  mI->id       = *id;
  mI->t        = drms_getkey_time(inRec, "ObservationTime",    &s1); s |= s1;
  mI->lat      = drms_getkey_int (inRec, "LatitudeHG",         &s1); s |= s1;
  mI->lon      = drms_getkey_int (inRec, "LongitudeCM",        &s1); s |= s1;
  mI->lon_carr = drms_getkey_int (inRec, "LongitudeHG",        &s1); s |= s1;
  mI->lon_wid  = drms_getkey_int (inRec, "LongitudinalExtent", &s1); s |= s1;
  mI->area     = drms_getkey_int (inRec, "Area",               &s1); s |= s1;
  mI->spot     = drms_getkey_int (inRec, "SpotCount",          &s1); s |= s1;
  // strings: be careful about lengths
  // note: strncpy does not ensure null termination, so use snprintf
  // 1: magnetic class
  str = drms_getkey_string(inRec, "MagneticType", &s1); s |= s1;
  if (str && (strcmp(str, "?") == 0))
    str = "Plage"; // use this instead of ?
  snprintf(mI->class_mag, sizeof(mI->class_mag), "%s", str ? str : "");
  // 2: zurich class
  str = drms_getkey_string(inRec, "ZurichClass", &s1);  s |= s1;
  snprintf(mI->class_zur, sizeof(mI->class_zur), "%s", str ? str : "");
  return s;
}

/* 
 * Get all Marp metadata.  Returns the number of Marps (nMarp), a pointer to
 * the individual Marp IDs, and sets up the contents of a Marp-info list containing
 * all Marp metadata.  MInfo contains just one AR appearance per day for NOAA ARs.
 * We allocate *Marps_p, but not *MInfo.  We write in to both.
 * Returns 0 for success (even if some attributes could not be found), 1 for error
 */
static
int
match_get_all_info(DRMS_RecordSet_t *rs, 
		   int *nMarp, 
		   marp_info_t **Marps_p,
		   match_info_t *MInfo)
{
  DRMS_Record_t *rec;   // DRMS record for one Marp at one time
  marp_info_t *Marps;   // all the Marps
  int MAXid = 32;       // number of Marps allocated so far in above
  int Nid = 0;          // number of Marps used so far
  int Nrec, i, j;       // record iteration
  int id1;              // match-ar ID code
  int already_seen;

  // sensible defaults in case of premature return
  *nMarp = 0;
  *Marps_p = NULL;
  if ((Marps = malloc(MAXid*sizeof(*Marps))) == NULL) {
    V_printf(-1, "", "Failed malloc for %d Marps\n", MAXid);
    return 1;
  }
  /* Loop over all records => all appearances of all Marps
   * Loop runs across all records, which are returned sorted by the
   * first prime key (ObservationTime), and within that, by ID.
   */
  Nrec = rs->n;
  for (i = 0; i < Nrec; i++) {
    rec = rs->records[i];
    if (match_get_single_info(&id1, MInfo+i, rec) != 0) {
      V_printf(-1, "", "Could not get some attributes for Match AR %d\n", id1);
      continue;
    }
    for (already_seen = 0, j = 0; j < MAXid && Marps[j].id > 0; j++)
      if (Marps[j].id == id1) {
        already_seen = 1;
	break;
      }
    // printf("[% 3d] ID = % 5d %s\n", i+1, id1, already_seen ? "" : "!");
    if (!already_seen) {
      // insert the new AR
      if (Nid + 1 == MAXid) {
	// oops, must enlarge the Marp list
        MAXid *= 4;
        Marps = realloc(Marps, MAXid*sizeof(*Marps));
	if (Marps == NULL) {
	  V_printf(-1, "", "Failed re-alloc for %d Marps\n", MAXid);
	  return 1;
	}
      }
      Marps[Nid++].id = id1;
    }
  }
  // set up outputs, return OK
  *nMarp = Nid;
  *Marps_p = Marps;
  return 0;
}

/*
 * Load matching ARs from DRMS, and resample these sightings to return
 * a 2d matrix of match-patches (match_info_t's) that closely parallels
 * the harp-patches.  The matrix contains nMarp x nRec_all entries, and
 * is allocated here.  Additionally, we return a nMarp-length list of
 * matching AR id's, which is analogous to the harp_info list.  Both
 * the match-patch list and the match-info list are allocated here and
 * should be freed by the caller.
 *
 * Returns 0 on success, nonzero for failure.
 */
static
int
match_load_resample(int nRec_all,
		    trec_info_t *tRec, 
		    int *nMarp_p, 
		    marp_info_t **marp_p, 
		    match_info_t **match_p)
{
  const char *matchQueryRoot = "su_rsb.NOAA_ActiveRegions";
  const int SECDAY = 3600*24;  // seconds per day
  DRMS_RecordSet_t *matchRS;   // matching record set
  int status;                  // error status
  char t0_S[24], t1_S[24];     // time range for query
  char matchQuery[STR_MAX];    // DRMS query
  marp_info_t *marpInfo;
  match_info_t *matchInfo_raw; // daily match-patches (from query)
  match_info_t *matchInfo;     // resampled match-patches (to return)
  int Nmatch_raw;              // total length of daily match-patches
  int nMarp;                   // matchInfo is nRec_all x nMarp
  int i, m, inx, id;           // counters

  // initialize now in case of error exit
  *nMarp_p = 0;
  *marp_p = NULL;
  *match_p = NULL;
  /*
   * 1: Load the list of daily ARs from DRMS into our own struct's.
   */
  // set up query for matching ARs (add one day in both directions)
  if (nRec_all > 0) {
    TIME t0 = tRec[0         ].t - SECDAY;
    TIME t1 = tRec[nRec_all-1].t + SECDAY;
    sprint_time(t0_S, t0, "TAI", 0);
    sprint_time(t1_S, t1, "TAI", 0);
  } else {
    snprintf(t0_S, sizeof(t0_S), "%s", "2011.02.14_00:00:00_TAI");
    snprintf(t1_S, sizeof(t1_S), "%s", "2011.02.14_00:00:00_TAI");
    V_printf(verbflag > 0, "", "No HARPs, using nominal time range for matches.\n");
  }
  V_printf(verbflag > 0, "", "Match start = %s\n", t0_S);
  V_printf(verbflag > 0, "", "Match end   = %s\n", t1_S);
  snprintf(matchQuery, sizeof(matchQuery), "%s[%s-%s]", matchQueryRoot, t0_S, t1_S);
  // get the list of ARs from DRMS
  status = 0;
  matchRS = drms_open_records(drms_env, matchQuery, &status);
  if (status != DRMS_SUCCESS) {
    V_printf(-1, "", "Could not get match AR records from DRMS (%s)\n", matchQuery);
    return 1;
  }
  Nmatch_raw = matchRS->n;
  V_printf(verbflag > 0, "", "Loading %s (%d records)\n", matchQuery, Nmatch_raw);
  // load match metadata in original (daily) form into matchInfo structures
  if ((matchInfo_raw = calloc(Nmatch_raw, sizeof(*matchInfo_raw))) == NULL) {
    V_printf(-1, "", "Could not allocate space for %d match records.\n", Nmatch_raw);
    drms_close_records(matchRS, DRMS_FREE_RECORD); // release records
    return 1;
  }
  // this does not currently error out if AR attributes were not found
  status = match_get_all_info(matchRS, &nMarp, &marpInfo, matchInfo_raw);
  drms_close_records(matchRS, DRMS_FREE_RECORD); // can release these now
  if (status != 0) {
    V_printf(-1, "", "Could not extract info from match AR records\n");
    return 1;
  }
  V_printf(verbflag > 1, "\t", "Found %d distinct match ARs\n", nMarp);
  /*
   * 2: Resample the list of daily ARs at the times given by the tRec array
   */
  /* Note: The nested loops below take O(nMarp * nRec_all * Nmatch_raw) time.
   * Nmatch_raw = total # daily AR records ~ O(nMarp*#days) ~ O(nMarp * nRec_all/120)
   * Therefore, this is "sort of" quadratic in both nMarp and nRec_all.
   * It could be linear in nRec_all by doing the scan more intelligently, but it 
   * does not matter for now because the numbers are low and the inner loop is fast.
   */
  // allocate one record for each time, for each Marp
  if ((matchInfo = calloc(nMarp*nRec_all, sizeof(*matchInfo))) == NULL) {
    free(matchInfo_raw);
    return 1;
  }
  for (m = 0; m < nMarp; m++) {
    id = marpInfo[m].id;
    // resample AR numbered "id"
    for (i = 0; i < nRec_all; i++) {
      inx = nRec_all * m + i; // offset into match-patch array (marp #m, time #i)
      // generate a new AR sample of a given id at a new time t
      //   (scans the whole matchInfo_raw list)
      status = match_ar_resample(id, tRec[i].t, matchInfo_raw, Nmatch_raw, matchInfo+inx);
      if (status != 0) {
	// not an error, just informational
	// printf("No match for %d at %s\n", id, tRec[i].str);
      } else {
	// match_ar_print(stdout, "\t", matchInfo+inx);
      }
    }
  }
  free(matchInfo_raw);
  // set up output args
  *nMarp_p = nMarp;
  *marp_p = marpInfo;
  *match_p = matchInfo;
  return 0;
}

