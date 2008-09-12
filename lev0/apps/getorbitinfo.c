/* 
 *  getorbitinfo.
 */

#include "jsoc_main.h"
#include "astro.h"

#define kORBSERIES "in"
#define kINTERPALG "alg"
#define kTGTTIMES "tgt"
#define kTESTHCI "t"

ModuleArgs_t module_args[] =
{
  {ARG_STRING, kORBSERIES, "sdo.fds_orbit_vectors", "Series containing orbit position and velocity vectors; it testing, then contains record-set query."},
  {ARG_STRING, kINTERPALG, "linear", "Supported algorithm to use when interpolating between grid times."},
  {ARG_DOUBLES, kTGTTIMES, "problem", "Array of internal times (doubles) for which orbit information is to be returned. "},
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
   int ntimes = 0;
   int doHCItest;

   doHCItest = cmdparams_isflagset(&cmdparams, kTESTHCI);
   orbseries = cmdparams_get_str(&cmdparams, kORBSERIES, NULL);
   
   if (!doHCItest)
   {
      alg = cmdparams_get_str(&cmdparams, kINTERPALG, NULL);
      ntimes = cmdparams_get_dblarr(&cmdparams, kTGTTIMES, &tgttimes, NULL);
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
      else
      {
         fprintf(stderr, "Unsupported interpolation algorithm '%s'.\n", alg);
         status = kSDOORB_failure;
      }

      if (status == kSDOORB_success)
      {
         if (iorbit_getinfo(drms_env, orbseries, interpalg, tgttimes, ntimes, &info) != kLIBASTRO_Success)
         {
            status = kSDOORB_failure;
         }

         iorbit_cleanup();
      }
   }

   return status;
}
