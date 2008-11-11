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
   double *tgttimes = NULL;
   const char *gridtimes = NULL;
   int ntimes = 0;
   int doHCItest;

   doHCItest = cmdparams_isflagset(&cmdparams, kTESTHCI);
   orbseries = cmdparams_get_str(&cmdparams, kORBSERIES, NULL);
   
   if (!doHCItest)
   {
      alg = cmdparams_get_str(&cmdparams, kINTERPALG, NULL);
      ntimes = cmdparams_get_dblarr(&cmdparams, kTGTTIMES, &tgttimes, NULL);
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
