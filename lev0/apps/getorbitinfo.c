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

#define kOBSDATE "OBS_DATE"
#define kHCIX "HCIEC_X"
#define kHCIY "HCIEC_Y"
#define kHCIZ "HCIEC_Z"
#define kHCIVX "HCIEC_VX"
#define kHCIVY "HCIEC_VY"
#define kHCIVZ "HCIEC_VZ"

ModuleArgs_t module_args[] =
{
  {ARG_STRING, kORBSERIES, "sdo.fds_orbit_vectors", "Series containing orbit position and velocity vectors; it testing, then contains record-set query."},
  {ARG_STRING, kINTERPALG, "linear", "Supported algorithm to use when interpolating between grid times."},
  {ARG_DOUBLES, kTGTTIMESARR, "0.0", "Array of internal times (doubles) for which orbit information is to be returned. "},
  {ARG_STRING, kTGTTIMES, "problem", "File containing array of internal times (doubles) for which orbit information is to be returned. "},
  {ARG_STRING, kGRIDTIMES, "unspecified", "Record-set query that specifies grid-vector internal times."},
  {ARG_FLAG, kTESTINTERP, "", "Test interpolation routine."},
  {ARG_FLAG, kNOVELOCITY, "", "The orbit series has no velocity keywores (like iris.orbitvectors)"},
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
   TIMER_t *tmr = NULL;

   orbseries = cmdparams_get_str(&cmdparams, kORBSERIES, NULL);
   testinterp = cmdparams_isflagset(&cmdparams, kTESTINTERP);
    novelocity = cmdparams_isflagset(&cmdparams, kNOVELOCITY);

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
               hcon_insert(keymap, "kXGCI", "geixobs");
               hcon_insert(keymap, "kYGCI", "geiyobs");
               hcon_insert(keymap, "kZGCI", "geizobs");
               hcon_insert(keymap, "kXHCI", "heixobs");
               hcon_insert(keymap, "kYHCI", "heiyobs");
               hcon_insert(keymap, "kZHCI", "heizobs");
               hcon_insert(keymap, "kRSUNOBS", "rsunobs");
               hcon_insert(keymap, "kOBSVR", "obsvr");
               hcon_insert(keymap, "kDSUNOBS", "dsunobs");
               hcon_insert(keymap, "kOBSDATE", "obsdate");
               
               rv = iorbit_getinfo_ext(drms_env, 
                                       orbseries,
                                       gridtimes, 
                                       interpalg, 
                                       tgttimes, 
                                       ntimes, 
                                       kIORBIT_CacheAction_DontCache,
                                       &info,
                                       keymap);
           }
           else
           {
               rv = kLIBASTRO_OutOfMemory;
           }
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
                      "obstime", "hciX-interp", "hciY-interp", "hciZ-interp", "gciX-interp", "gciY-interp", "gciZ-interp", "rsunobs-interp", "obsvr-interp", "dsunobs-interp");
              
              for (iinfo = 0; iinfo < ntimes; iinfo++)
              {
                  /* interpolated values */
                  infoitem = &(info[iinfo]);
                  sprint_time(timestr, infoitem->obstime, "UTC", 0);
                  
                  if (!testinterp)
                  {
                      fprintf(stdout, "%-24s%-22.4f%-22.4f%-22.4f%-22.4f%-22.4f%-22.4f%-22.4f%-22.4f%-22.4f\n", 
                              timestr, infoitem->hciX, infoitem->hciY, infoitem->hciZ, 
                              infoitem->gciX, infoitem->gciY, infoitem->gciZ, infoitem->rsun_obs, infoitem->obs_vr, infoitem->dsun_obs);
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
                          "obstime", "secs", "hciX-interp", "hciY-interp", "hciZ-interp", "hciVX-interp", "hciVY-interp", "hciVZ-interp");
              }
              else
              {
                  fprintf(stdout, "%-24s%-12s%-22s%-22s%-22s%-22s%-22s%-22s%-20s%-20s%-20s%-14s%-14s%-14s\n", 
                          "obstime", "secs", 
                          "hciX-interp", "hciX-actual",
                          "hciY-interp", "hciY-actual",
                          "hciZ-interp", "hciZ-actual",
                          "hciVX-interp", "hciVX-actual",
                          "hciVY-interp", "hciVY-actual",
                          "hciVZ-interp", "hciVZ-actual");
              }
              
              for (iinfo = 0; iinfo < ntimes; iinfo++)
              {
                  /* interpolated values */
                  infoitem = &(info[iinfo]);
                  sprint_time(timestr, infoitem->obstime, "UTC", 0);
                  
                  if (!testinterp)
                  {
                      fprintf(stdout, "%-24s%-14.1f%-22.4f%-22.4f%-22.4f%-22.4f%-22.4f%-22.4f\n", 
                              timestr, infoitem->obstime, infoitem->hciX, infoitem->hciY, infoitem->hciZ, 
                              infoitem->hciVX, infoitem->hciVY, infoitem->hciVZ);
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
                      double actual_hciX = 0;
                      double actual_hciY = 0;
                      double actual_hciZ = 0;
                      double actual_hciVX = 0;
                      double actual_hciVY = 0;
                      double actual_hciVZ = 0;
                      
                      snprintf(query, sizeof(query), "%s[%s]", orbseries, timestr);
                      
                      rs = drms_open_records(drms_env, query, &drmsstatus);
                      if (rs && rs->n == 1)
                      {
                          obstime = drms_getkey_time(rs->records[0], kOBSDATE, &drmsstatus);
                          
                          if (obstime == infoitem->obstime)
                          {
                              actual_hciX = drms_getkey_double(rs->records[0], kHCIX, &drmsstatus);
                              actual_hciY = drms_getkey_double(rs->records[0], kHCIY, &drmsstatus);
                              actual_hciZ = drms_getkey_double(rs->records[0], kHCIZ, &drmsstatus);
                              actual_hciVX = drms_getkey_double(rs->records[0], kHCIVX, &drmsstatus);
                              actual_hciVY = drms_getkey_double(rs->records[0], kHCIVY, &drmsstatus);
                              actual_hciVZ = drms_getkey_double(rs->records[0], kHCIVZ, &drmsstatus);
                          }
                      }
                      
                      fprintf(stdout, "%-24s%-12.1f%-22.8f%-22.8f%-22.8f%-22.8f%-22.8f%-22.8f%-22.8f%-22.8f%-22.8f%-14.8f%-14.8f%-14.8f\n", 
                              timestr, infoitem->obstime, 
                              infoitem->hciX, actual_hciX, 
                              infoitem->hciY, actual_hciY,
                              infoitem->hciZ, actual_hciZ,
                              infoitem->hciVX, actual_hciVX,
                              infoitem->hciVY, actual_hciVY,
                              infoitem->hciVZ, actual_hciVZ);
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
