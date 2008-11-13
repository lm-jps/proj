/* 
 *  getorbitinfo.
 */

#include "jsoc_main.h"
#include "astro.h"

#define kORBSERIES "in"
#define kINTERPALG "alg"
#define kTGTTIMES "tgt"
#define kGRIDTIMES "grid"
#define kTESTHCI "t"

ModuleArgs_t module_args[] =
{
  {ARG_STRING, kORBSERIES, "sdo.fds_orbit_vectors", "Series containing orbit position and velocity vectors; it testing, then contains record-set query."},
  {ARG_STRING, kINTERPALG, "linear", "Supported algorithm to use when interpolating between grid times."},
  {ARG_DOUBLES, kTGTTIMES, "problem", "Array of internal times (doubles) for which orbit information is to be returned. "},
  {ARG_STRING, kGRIDTIMES, "unspecified", "Record-set query that specifies grid-vector internal times."},
  {ARG_FLAG, kTESTHCI, "", "Test MOC HCI values."},
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
   LinkedList_t *info = NULL;
   IORBIT_Alg_t interpalg;
   char *orbseries = NULL;
   char *alg = NULL;
   char *tgttimesstr = NULL;
   double *tgttimes = NULL;
   const char *gridtimes = NULL;
   int ntimes = 0;
   int doHCItest;

   doHCItest = cmdparams_isflagset(&cmdparams, kTESTHCI);
   orbseries = cmdparams_get_str(&cmdparams, kORBSERIES, NULL);
   
   if (!doHCItest)
   {
      alg = cmdparams_get_str(&cmdparams, kINTERPALG, NULL);

      /* cmdparams or hash_insert has some kind of bug in it - if the string value of a cmdparams is too
       * long, cmdparams will fail. */
      tgttimesstr = cmdparams_get_str(&cmdparams, kTGTTIMES, NULL);
      if (tgttimesstr[0] == '@')
      {
         /* read in the tgttimes from a file */
         struct stat stBuf;
         FILE *atfile = NULL;
         char *filestr = NULL;
         char lineBuf[LINE_MAX];
         char *fullline = NULL;
         int stgt = 48;

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
      }
      else
      {
         ntimes = cmdparams_get_dblarr(&cmdparams, kTGTTIMES, &tgttimes, NULL);
      }

      gridtimes = cmdparams_get_str(&cmdparams, kGRIDTIMES, NULL);

      if (strcasecmp(gridtimes, "unspecified") == 0)
      {
         gridtimes = NULL;
      }
   }

   if (doHCItest)
   {
      if (iorbit_test(drms_env, orbseries) != kLIBASTRO_Success)
      {
         status = kSDOORB_failure;
      }
   }
   else
   {
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
         if (iorbit_getinfo(drms_env, 
                            orbseries, 
                            interpalg, 
                            tgttimes, 
                            ntimes, 
                            gridtimes,
                            &info) != kLIBASTRO_Success)
         {
            status = kSDOORB_failure;
         }
         else
         {
            /* This is a demonstration module - just print */
            ListNode_t *node = NULL;
            IORBIT_Info_t *infoitem = NULL;
            char timestr[128];

            fprintf(stdout, "%-32s%-20s%-20s%-20s%-20s\n", 
                    "obstime", "dsun_obs", "obs_vr", "obs_vw", "obs_vn");

            list_llreset(info);
            while ((node = list_llnext(info)) != NULL)
            {
               infoitem = (IORBIT_Info_t *)(node->data);
               sprint_time(timestr, infoitem->obstime, "UTC", 0);

              
               fprintf(stdout, "%-32s%-20.8f%-20.8f%-20.8f%-20.8f\n", 
                       timestr, infoitem->dsun_obs, infoitem->obs_vr, infoitem->obs_vw, infoitem->obs_vn);
            }

            /* Now, let's print out interpolated heliocentric vectors for kicks */
            fprintf(stdout, "\n\n");
            fprintf(stdout, "%-32s%-12s%-20s%-20s%-20s%-20s%-20s%-20s\n", 
                    "obstime", "secs", "hciX", "hciY", "hciZ", "hciVX", "hciVY", "hciVZ");
            list_llreset(info);
            while ((node = list_llnext(info)) != NULL)
            {
               infoitem = (IORBIT_Info_t *)(node->data);
               sprint_time(timestr, infoitem->obstime, "UTC", 0);

              
               fprintf(stdout, "%-32s%-12.1f%-20.8f%-20.8f%-20.8f%-20.8f%-20.8f%-20.8f\n", 
                       timestr, infoitem->obstime, infoitem->hciX, infoitem->hciY, infoitem->hciZ, 
                       infoitem->hciVX, infoitem->hciVY, infoitem->hciVZ);
            }
         }

         iorbit_cleanup();
      }
   }

   if (info)
   {
      list_llfree(&info);
   }

   return status;
}
