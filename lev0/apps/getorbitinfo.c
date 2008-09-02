/* 
 *  getorbitinfo.
 */

#include "jsoc_main.h"
#include "astro.h"

#define kRECSETIN "in"
#define kTESTHCI "t"

ModuleArgs_t module_args[] =
{
  {ARG_STRING, kRECSETIN, "", "Set of records containing orbit position and velocity vectors."},
  {ARG_FLAG, kTESTHCI, "", "Test MOC HCI values."},
  {ARG_END}
};

char *module_name = "getorbitinfo";

enum sdoorb_stat_enum
{
   kSDOORB_success = 0,
   kSDOORB_libastrofailure
};

typedef enum sdoorb_stat_enum sdoorb_stat_t;

int DoIt(void)
{
   sdoorb_stat_t status = kSDOORB_success;
   LinkedList_t *info = NULL;
   int doHCItest = 0;

   char *recin = cmdparams_get_str(&cmdparams, kRECSETIN, NULL);
   doHCItest = cmdparams_isflagset(&cmdparams, kTESTHCI);

   if (doHCItest)
   {
      if (testiorbit(drms_env, recin, &info) != kLIBASTRO_Success)
      {
         status = kSDOORB_libastrofailure;
      }
   }
   else
   {
      if (iorbit(drms_env, recin, &info) != kLIBASTRO_Success)
      {
         status = kSDOORB_libastrofailure;
      }
   }

   return status;
}
