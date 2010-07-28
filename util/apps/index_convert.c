/*
 *
 * index_convert - Converts a slot value to an index, or vise versa, for a given series.
 *
 */

/**
\defgroup index_convert - Convert index to key value or the reverse
@ingroup drms_util

\par Synopsis:
index_convert [-h]
index_convert ds=<seriesname> [<pkey>=<value>] | [<pkey>_index=<indexvalue>] 
index_convert {same options as above}
\endcode

\details

\b Index_convert converts a prime key index value to a prime key value or the reverse.  The series must be
slotted.  

The operation is to get the template record then if the <pkey>_index aregument is present convert
that value to the equivalent prime key value.  If the <pkey> argument is present, i.e. no "_index" in
the keyname, then the equivalent index value is computed and printed.

\par Options:

\par Flags:

\li \c -h: help - print usage information and exit

\par Parameters:

\li \c ds=<seriesname> - required parameter specified the series to examine.
\li \c <pkey>=<value> - optional prime key name and its test value.
\li \c <pkey>_index=<indexvalue> - optional prime key slotted index name and its test index value.

\par JSOC flags:
\ref jsoc_main

\par Usage:

The program may be used to convert a single prime key value into its equivalent index value or the reverse.
The keyname provided must be for a slotted prime key in the specified series.  If the seriesname is
not a valid series or if the keyname is not a prime key in that series an error code is returned.
The direction of the conversion is determined by the presence or absence of the "_index" suffix on
the single specified keyword argument.

\par Examples:

\b Example 1:
\code
  index_convert ds=mdi.fd_V_lev18 T_REC=1996.05.01
\endcode
reports: 1751040, the minute number for the given time using the proper epoch.

\b Example 2:
\code
  index_convert ds=mdi.fd_V_lev18 T_REC_index=1751040
\endcode
reports: 1996.05.01_00:00:00_TAI, the time of the center of the specified slot.

\bug

\sa
time_convert time_index(SOI)

*/

#include "jsoc_main.h"
#include "drms.h"
#include "atoinc.h"

#define NOT_SPECIFIED	"NOT SPECIFIED"

ModuleArgs_t module_args[] =
{ 
    {ARG_STRING, "ds", NOT_SPECIFIED,  "Input data series."},
    {ARG_FLAG, "h", "0", "Help - Print usage and exit"},
    {ARG_END}
};

#define DIE(msg) {fprintf(stderr,"%s\n",msg);exit(1);}

char *module_name = "index_convert";

int nice_intro ()
  {
  int usage = cmdparams_get_int (&cmdparams, "h", NULL);
  if (usage)
    {
    printf ("Usage:\nndex_convert [-h] "
        "ds=<seriesname> {<pkey>=<value>} | {<pkey>_index=<indexvalue>} \n"
        "  -h: help - show this message then exit\n"
        "ds=<seriesname> - required\n"
        "<pkey> - prime key to use\n"
        "<pkey>_index - prime index key to use\n"
        "Only one of <pkey> or <pkey>_index is allowed\n"
        );
    return(1);
    }
  return (0);
  }

/* Module main function. */
int DoIt(void)
  {
  int status = 0;
  DRMS_Record_t *template;
  DRMS_Keyword_t *skey, *pkey;
  DRMS_Type_t ptype;
  int npkeys;
  char *pname=NULL;
  char *piname=NULL;
  char *seriesname;
  char *inbracket;
  const char *ds = cmdparams_get_str (&cmdparams, "ds", NULL);
  int ikey;

  if (nice_intro()) return(0);

  /* check for minimum inputs */
  if (strcmp(ds, NOT_SPECIFIED) == 0 )
    DIE("No dataseries: at least ds must be specified");

  /* get series info and prime key type, etc. */
  inbracket = index(ds, '[');
  if (inbracket)
	  *inbracket = '\0';
  template = drms_template_record (drms_env, ds, &status);
  if (!template || status)
	DIE("Series not found");

  npkeys = template->seriesinfo->pidx_num;
  if (npkeys < 1)
    DIE("Series has no prime keys");

  // for each prime key check command line for key or key_index args
  // ignore non-slotted prime keys
  // The first appropriate key is used.
  for (ikey=0; ikey < npkeys; ikey++)
	{
	pkey = template->seriesinfo->pidx_keywords[ikey];
	if (pkey->info->recscope > 1)
		{
		if (cmdparams_exists(&cmdparams, pkey->info->name))
			{
			// pkey->value.longlong_val = params_get_int64(&cmdparams, pkey->info->name);
			skey = drms_keyword_slotfromindex(pkey);
			long long indexval = params_get_int64(&cmdparams, pkey->info->name);
			double epoch = drms_keyword_getslotepoch(skey, &status);
			double step = drms_keyword_getvalkeystep(skey, &status);
			skey->value.double_val = indexval * step + epoch;
			drms_keyword_printval(skey);
			printf("\n");
			return(0);
			}
		skey = drms_keyword_slotfromindex(pkey);
		if (cmdparams_exists(&cmdparams, skey->info->name))
			{
		        DRMS_Value_t indexval;
		        DRMS_Value_t inval;
		        if (skey->info->type == DRMS_TYPE_TIME)
			    inval.value.time_val = params_get_time(&cmdparams, skey->info->name);
		        else 
			    inval.value.double_val = params_get_double(&cmdparams, skey->info->name);
		        inval.type = skey->info->type;
		        drms_keyword_slotval2indexval(skey, &inval, &indexval, NULL);
		        pkey->value.longlong_val = indexval.value.longlong_val;
		        drms_keyword_printval(pkey);
			printf("\n");
			return(0);
			}
		}
	}
  return(1);
  }
