/*
 *  getorbitinfo.
 */

#include "jsoc_main.h"
#include "astro.h"

#define kORBSERIES "in"
#define kINTERPALG "alg"
#define kTGTTIMESARR "tgtarr"
#define kTGTTIMES "tgt"
#define kGRIDTIMES "grid"
#define kTESTINTERP "i"
#define kNOVELOCITY "n"
#define kSAAHLZ "s"


ModuleArgs_t module_args[] =
{
    {ARG_STRING, kORBSERIES, "sdo.fds_orbit_vectors", "Series containing orbit position and velocity vectors; it testing, then contains record-set query."},
    {ARG_STRING, kINTERPALG, "linear", "Supported algorithm to use when interpolating between grid times."},
    {ARG_DOUBLES, kTGTTIMESARR, "0.0", "Array of internal times (doubles) for which orbit information is to be returned. "},
    {ARG_STRING, kTGTTIMES, "problem", "File containing array of internal times (doubles) for which orbit information is to be returned. "},
    {ARG_STRING, kGRIDTIMES, "unspecified", "Record-set query that specifies grid-vector internal times."},
    {ARG_FLAG, kTESTINTERP, "", "Test interpolation routine."},
    {ARG_FLAG, kNOVELOCITY, "", "The orbit series has no velocity keywords (like iris.orbitvectors)"},
    {ARG_FLAG, kSAAHLZ, "", "Test the SAA-HLZ stuff out."},
    {ARG_END}
};

char *module_name = "getorbitinfo";

enum sdoorb_stat_enum
{
   kSDOORB_success = 0,
   kSDOORB_failure
};

typedef enum sdoorb_stat_enum sdoorb_stat_t;

int DoIt(void)
{
   sdoorb_stat_t status = kSDOORB_success;
   IORBIT_Info_t *info = NULL;
   IORBIT_Alg_t interpalg;
   const char *orbseries = NULL;
   const char *alg = NULL;
   const char *tgttimesstr = NULL;
   double *tgttimes = NULL;
   const char *gridtimes = NULL;
   int ntimes = 0;
   int testinterp;
    int novelocity;
    int saahlz;
   TIMER_t *tmr = NULL;

   orbseries = cmdparams_get_str(&cmdparams, kORBSERIES, NULL);
   testinterp = cmdparams_isflagset(&cmdparams, kTESTINTERP);
    novelocity = cmdparams_isflagset(&cmdparams, kNOVELOCITY);
    saahlz = cmdparams_isflagset(&cmdparams, kSAAHLZ);

   alg = cmdparams_get_str(&cmdparams, kINTERPALG, NULL);

   /* cmdparams or hash_insert has some kind of bug in it - if the string value of a cmdparams is too
    * long, cmdparams will fail. */
   tgttimesstr = cmdparams_get_str(&cmdparams, kTGTTIMES, NULL);
   if (tgttimesstr[0] == '@')
   {
      /* read in the tgttimes from a file */
      struct stat stBuf;
      FILE *atfile = NULL;
      const char *filestr = NULL;
      char lineBuf[LINE_MAX];
      char *fullline = NULL;
      int stgt = 48;

      if (drms_env->verbose)
      {
         /* time file read */
         tmr = CreateTimer();
      }

      filestr = tgttimesstr + 1;

      if (!stat(filestr, &stBuf))
      {
         if (S_ISREG(stBuf.st_mode))
         {
            /* read a line */
            if ((atfile = fopen(filestr, "r")) == NULL)
            {
               fprintf(stderr, "Cannot open @file %s for reading, skipping.\n", filestr);
            }
            else
            {
               int len = 0;
               int itgt = 0;
               tgttimes = malloc(sizeof(double) * stgt);

               while (!(fgets(lineBuf, LINE_MAX, atfile) == NULL))
               {
                  /* strip \n from end of lineBuf */
                  len = strlen(lineBuf);

                  fullline = strdup(lineBuf);

                  if (len == LINE_MAX - 1)
                  {
                     /* may be more on this line */
                     while (!(fgets(lineBuf, LINE_MAX, atfile) == NULL))
                     {
                        fullline = realloc(fullline, strlen(fullline) + strlen(lineBuf) + 1);
                        snprintf(fullline + strlen(fullline),
                                 strlen(lineBuf) + 1,
                                 "%s",
                                 lineBuf);
                        if (strlen(lineBuf) > 1 && lineBuf[strlen(lineBuf) - 1] == '\n')
                        {
                           break;
                        }
                     }
                  }

                  len = strlen(fullline);

                  if (fullline[len - 1] == '\n')
                  {
                     fullline[len - 1] = '\0';
                  }

                  /* fullline has a single time (string) */
                  if (itgt == stgt)
                  {
                     stgt *= 2;
                     tgttimes = realloc(tgttimes, sizeof(double) * stgt);
                  }

                  sscanf(fullline, "%lf", &(tgttimes[itgt]));

                  if (fullline)
                  {
                     free(fullline);
                     fullline = NULL;
                  }

                  itgt++;
               }

               ntimes = itgt;
            }
         }
      }

      if (drms_env->verbose)
      {
         fprintf(stdout, "file read seconds elapsed: %f\n", GetElapsedTime(tmr));
         DestroyTimer(&tmr);
      }
   }
   else
   {
      ntimes = cmdparams_get_dblarr(&cmdparams, kTGTTIMESARR, &tgttimes, NULL);
   }

   gridtimes = cmdparams_get_str(&cmdparams, kGRIDTIMES, NULL);

   if (strcasecmp(gridtimes, "unspecified") == 0)
   {
      gridtimes = NULL;
   }

   if (strcasecmp(alg, "linear") == 0)
   {
      interpalg = IORBIT_Alg_Linear;
   }
   else if (strncasecmp(alg, "quad", strlen("quad")) == 0)
   {
      interpalg = IORBIT_Alg_Quadratic;
   }
   else
   {
      fprintf(stderr, "Unsupported interpolation algorithm '%s'.\n", alg);
      status = kSDOORB_failure;
   }

   if (status == kSDOORB_success)
   {
       int rv = 0;

      if (drms_env->verbose)
      {
         /* time solar velocity calculation */
         tmr = CreateTimer();
      }

       if (novelocity)
       {
           HContainer_t *keymap = NULL;

           keymap = hcon_create(DRMS_MAXKEYNAMELEN, DRMS_MAXKEYNAMELEN, NULL, NULL, NULL, NULL, 0);

           if (keymap)
           {
               hcon_insert(keymap, IORBIT_STRINGIFY(IORBIT_KEYWORD_GAE_X), "geixobs");
               hcon_insert(keymap, IORBIT_STRINGIFY(IORBIT_KEYWORD_GAE_Y), "geiyobs");
               hcon_insert(keymap, IORBIT_STRINGIFY(IORBIT_KEYWORD_GAE_Z), "geizobs");
               hcon_insert(keymap, IORBIT_STRINGIFY(IORBIT_KEYWORD_HAE_X), "heixobs");
               hcon_insert(keymap, IORBIT_STRINGIFY(IORBIT_KEYWORD_HAE_Y), "heiyobs");
               hcon_insert(keymap, IORBIT_STRINGIFY(IORBIT_KEYWORD_HAE_Z), "heizobs");
               hcon_insert(keymap, IORBIT_STRINGIFY(IORBIT_KEYWORD_RSUN_OBS), "rsunobs");
               hcon_insert(keymap, IORBIT_STRINGIFY(IORBIT_KEYWORD_OBS_VR), "obsvr");
               hcon_insert(keymap, IORBIT_STRINGIFY(IORBIT_KEYWORD_DSUN_OBS), "dsunobs");
               hcon_insert(keymap, IORBIT_STRINGIFY(IORBIT_KEYWORD_OBS_DATE), "obsdate");

               rv = iorbit_getinfo_ext(drms_env,
                                       orbseries,
                                       gridtimes,
                                       interpalg,
                                       tgttimes,
                                       ntimes,
                                       kIORBIT_CacheAction_DontCache,
                                       &info,
                                       keymap);

               hcon_destroy(&keymap);
           }
           else
           {
               rv = kLIBASTRO_OutOfMemory;
           }
       }
       else if (saahlz)
       {
           /* Test out the SAA-HLZ stuff. */
           IORBIT_SaaHlzInfo_t *saahlzinfo = NULL;
           int itime;
           LinkedList_t *list = NULL;
           LinkedList_t **pList = NULL;
           ListNode_t *node = NULL;
           char *eventType = NULL;
           int nhlz;
           int shlz;
           int saa;
           HContainer_t *colToList = NULL;
           HContainer_t **pColToList = NULL;
           char timeStr[IORBIT_SAAHLZINFO_TIME_KEY_LEN];

           rv = iorbit_getSaaHlzInfo(drms_env, "iris.saa_hlz", tgttimes, ntimes, &saahlzinfo);

           if (rv == kLIBASTRO_Success)
           {
               for (itime = 0; itime < ntimes; itime++)
               {
                   nhlz = 0;
                   shlz = 0;
                   saa = 0;

                   snprintf(timeStr, sizeof(timeStr), "%lf", tgttimes[itime]);

                   if ((pColToList = hcon_lookup(saahlzinfo, timeStr)) != NULL)
                   {
                       colToList = *pColToList;

                       if ((pList = hcon_lookup(colToList, IORBIT_SAAHLZINFO_KW_EVENT_TYPE)) != NULL)
                       {
                           list = *pList;
                           list_llreset(list);

                           while((node = list_llnext(list)) != NULL)
                           {
                               eventType = (char *)node->data;

                               if (strcasecmp(eventType, "NHLZ") == 0)
                               {
                                   nhlz = 1;
                               }
                               else if (strcasecmp(eventType, "SHLZ") == 0)
                               {
                                   shlz = 1;
                               }
                               else if (strcasecmp(eventType, "SAA") == 0)
                               {
                                   saa = 1;
                               }
                           }

                           printf("For time %lf:\n", tgttimes[itime]);

                           if (nhlz && shlz)
                           {
                               /* ERROR - set HLZ to 0. */
                               printf("  ERROR - setting HLZ to 0.\n");
                           }
                           else if (nhlz)
                           {
                               printf("  Setting HLZ to 1.\n");
                           }
                           else if (shlz)
                           {
                               printf("  Setting HLZ to 2.\n");
                           }
                           else
                           {
                               printf("  Setting HLZ to 0.\n");
                           }

                           if (saa)
                           {
                               printf("  Setting SAA to 1.\n");
                           }
                           else
                           {
                               printf("  Setting SAA to 0.\n");
                           }
                       }
                   }
               }

               iorbit_cleanSaaHlzInfo(&saahlzinfo);
           }
           else
           {
               printf("Couldn't query db properly.\n");
               status = kSDOORB_failure;
           }

           return status;
       }
       else
       {
           rv = iorbit_getinfo(drms_env,
                          orbseries,
                          gridtimes,
                          interpalg,
                          tgttimes,
                          ntimes,
                          kIORBIT_CacheAction_DontCache,
                          &info);
       }

      if (rv != kLIBASTRO_Success)
      {
         status = kSDOORB_failure;
      }
      else
      {
          if (drms_env->verbose)
          {
              fprintf(stdout, "iorbit_getinfo() seconds elapsed: %f\n", GetElapsedTime(tmr));
              DestroyTimer(&tmr);
          }

          /* This is a demonstration module - just print */
          IORBIT_Info_t *infoitem = NULL;
          char timestr[128];
          int iinfo;

          if (novelocity)
          {
              fprintf(stdout, "%-24s%-22s%-22s%-22s%-22s%-22s%-22s%-22s%-22s%-22s\n",
                      "obstime", "hae_x-interp", "hae_y-interp", "hae_z-interp", "gciX-interp", "gciY-interp", "gciZ-interp", "rsunobs-interp", "obsvr-interp", "dsunobs-interp");

              for (iinfo = 0; iinfo < ntimes; iinfo++)
              {
                  /* interpolated values */
                  infoitem = &(info[iinfo]);
                  sprint_time(timestr, infoitem->obstime, "UTC", 0);

                  if (!testinterp)
                  {
                      fprintf(stdout, "%-24s%-22.4f%-22.4f%-22.4f%-22.4f%-22.4f%-22.4f%-22.4f%-22.4f%-22.4f\n",
                              timestr, infoitem->hae_x, infoitem->hae_y, infoitem->hae_z,
                              infoitem->gae_x, infoitem->gae_y, infoitem->gae_z, infoitem->rsun_obs, infoitem->obs_vr, infoitem->dsun_obs);
                  }
              }
          }
          else
          {
              fprintf(stdout, "%-32s%-20s%-15s%-15s%-15s%-12s%-12s%-12s%-10s%-32s\n",
                      "obstime", "dsun_obs", "obs_vr", "obs_vw", "obs_vn", "rsun_obs", "crln_obs", "crlt_obs", "car_rot", "orb_rec");

              for (iinfo = 0; iinfo < ntimes; iinfo++)
              {
                  infoitem = &(info[iinfo]);
                  sprint_time(timestr, infoitem->obstime, "UTC", 0);

                  fprintf(stdout, "%-32s%-20.3f%-15.3f%-15.3f%-15.3f%-12.3f%-12.3f%-12.3f%-10d%-32s\n",
                          timestr, infoitem->dsun_obs, infoitem->obs_vr, infoitem->obs_vw, infoitem->obs_vn,
                          infoitem->rsun_obs, infoitem->crln_obs, infoitem->crlt_obs, infoitem->car_rot, infoitem->orb_rec);
              }

              /* Now, let's print out interpolated heliocentric vectors for kicks */
              fprintf(stdout, "\n\n");

              if (!testinterp)
              {
                  fprintf(stdout, "%-24s%-14s%-22s%-22s%-22s%-22s%-22s%-22s\n",
                          "obstime", "secs", "hae_x-interp", "hae_y-interp", "hae_z-interp", "hae_vx-interp", "hae_vy-interp", "hae_vz-interp");
              }
              else
              {
                  fprintf(stdout, "%-24s%-12s%-22s%-22s%-22s%-22s%-22s%-22s%-20s%-20s%-20s%-14s%-14s%-14s\n",
                          "obstime", "secs",
                          "hae_x-interp", "hae_x-actual",
                          "hae_y-interp", "hae_y-actual",
                          "hae_z-interp", "hae_z-actual",
                          "hae_vx-interp", "hae_vx-actual",
                          "hae_vy-interp", "hae_vy-actual",
                          "hae_vz-interp", "hae_vz-actual");
              }

              for (iinfo = 0; iinfo < ntimes; iinfo++)
              {
                  /* interpolated values */
                  infoitem = &(info[iinfo]);
                  sprint_time(timestr, infoitem->obstime, "UTC", 0);

                  if (!testinterp)
                  {
                      fprintf(stdout, "%-24s%-14.1f%-22.4f%-22.4f%-22.4f%-22.4f%-22.4f%-22.4f\n",
                              timestr, infoitem->obstime, infoitem->hae_x, infoitem->hae_y, infoitem->hae_z,
                              infoitem->hae_vx, infoitem->hae_vy, infoitem->hae_vz);
                  }
                  else
                  {
                      /* actual values - this will only work if the tgttimes are whole seconds, and they
                       * are in the range of the grid vectors. If so, print out the actual values
                       * because the user is obviously testing this module.
                       * Use psql to find the actual value given timestr - not the most efficient way to do this
                       * but this is just a test anyway.
                       */
                      char query[DRMS_MAXQUERYLEN];
                      int drmsstatus = DRMS_SUCCESS;
                      DRMS_RecordSet_t *rs = NULL;
                      TIME obstime;
                      double actual_hae_x = 0;
                      double actual_hae_y = 0;
                      double actual_hae_z = 0;
                      double actual_hae_vx = 0;
                      double actual_hae_vy = 0;
                      double actual_hae_vz = 0;

                      snprintf(query, sizeof(query), "%s[%s]", orbseries, timestr);

                      rs = drms_open_records(drms_env, query, &drmsstatus);
                      if (rs && rs->n == 1)
                      {
                          obstime = drms_getkey_time(rs->records[0], IORBIT_KEYWORD_OBS_DATE, &drmsstatus);

                          if (obstime == infoitem->obstime)
                          {
                              actual_hae_x = drms_getkey_double(rs->records[0], IORBIT_KEYWORD_HAE_X, &drmsstatus);
                              actual_hae_y = drms_getkey_double(rs->records[0], IORBIT_KEYWORD_HAE_Y, &drmsstatus);
                              actual_hae_z = drms_getkey_double(rs->records[0], IORBIT_KEYWORD_HAE_Z, &drmsstatus);
                              actual_hae_vx = drms_getkey_double(rs->records[0], IORBIT_KEYWORD_HAE_VX, &drmsstatus);
                              actual_hae_vy = drms_getkey_double(rs->records[0], IORBIT_KEYWORD_HAE_VY, &drmsstatus);
                              actual_hae_vz = drms_getkey_double(rs->records[0], IORBIT_KEYWORD_HAE_VZ, &drmsstatus);
                          }
                      }

                      fprintf(stdout, "%-24s%-12.1f%-22.8f%-22.8f%-22.8f%-22.8f%-22.8f%-22.8f%-22.8f%-22.8f%-22.8f%-14.8f%-14.8f%-14.8f\n",
                              timestr, infoitem->obstime,
                              infoitem->hae_x, actual_hae_x,
                              infoitem->hae_y, actual_hae_y,
                              infoitem->hae_z, actual_hae_z,
                              infoitem->hae_vx, actual_hae_vx,
                              infoitem->hae_vy, actual_hae_vy,
                              infoitem->hae_vz, actual_hae_vz);
                  }
              }
          }

          /* clean up info */
          if (info)
          {
              free(info);
          }
      }

      iorbit_cleanup();
   }

   return status;
}
