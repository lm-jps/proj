/* this block of code sets the mech table info for HMI
   it assumes the mech tables are located at $JSOCROOT/TABLEPATH
   and are named in the POL_TABLE, TUNE_TABLE, FOCUS_TABLE defines.
   declaration below.

   this code creates the kwywords in the lines:
	static char *pol_keys[] =     {"HPL1POS", "HPL2POS", "HPL3POS"};
	static char *tuning_keys[] = {"HWL1POS", "HWL2POS", "HWL3POS", "HWL4POS"};
	static char *focus_keys[] =  {"HCF1POS", "HCF2POS"};

 Call: status = set_HMI_mech_values(DRMS_Record_t *rec)

 rec should be the lev0 record after the ISP values are inserted.

 status != 0 probably should cause a severe warning that the mechanism
 position keywords are wrong.
*/

#include <drms_keyword.h>
#include <printk.h>

/* Function to set HMI mechanism position keywords */

/* The "mech" tables are found in $JSOCROOT at these locations */
/* #define TABLE_PATH "proj/lev0/apps/data/" */
/* path now absolute as users don't have a Development directory */
#define TABLE_PATH "/home/jsoc/cvs/Development/JSOC/proj/tables/hmi_mech/"
/* #define POL_TABLE   "in_air_cal3.p" */
/* #define TUNE_TABLE  "in_air_cal3.w" */
/* #define FOCUS_TABLE "in_air_cal.c" */
#define POL_TABLE   "std_flight.p"
#define TUNE_TABLE  "std_flight.w"
#define FOCUS_TABLE "std_flight.c"


#define MECHMAXROWS	2000
#define MECHINIT_OK     0
#define MECHINIT_NOFILE	1
#define MECHINIT_ROWS   2
#define MECH_KEY_MISSING 3
#define MECH_INDEX_MISSING 4

typedef struct mech_tabinfo
  {
  char *filename;
  char *index;
  char **keys;
  int *table;
  int cols;
  } MECH_TABINFO_t;


/* Initialize mechanism lookup tables.  Is used once on forst call
   of  set_HMI_mech_values and once if out of range index is found
   to reset if tables have been updated.
*/

static int init_HMI_mech_tables(MECH_TABINFO_t *tabinfo, int tab)
{
for (tab=0; tab<3; tab++)
    {
    char tablepath[1024];
    char *tableroot;
    int idx, *res, val, vals;
    FILE *fp;
    char line[1024];
    /* get file for this table */
    //if(!(tableroot = (char *)getenv("JSOCROOT"))) return(MECHINIT_NOFILE);
    //strcpy(tablepath, tableroot);
    //strcat(tablepath, "/");
    //strcat(tablepath, TABLE_PATH);
    strcpy(tablepath, TABLE_PATH);
    strcat(tablepath, tabinfo[tab].filename);
    fp = fopen(tablepath, "r");
    if (!fp)
      {
      printk("Failed to open mech table %s, die.\n",tablepath);
      return(MECHINIT_NOFILE);
      }
    res = tabinfo[tab].table;
    vals = tabinfo[tab].cols;
    /* fill table with missing entry flags */
    for (idx=0; idx<MECHMAXROWS; idx++)
     for (val=0; val<vals+1; val++)
       res[val + (vals+1)*idx] = -1;
    /* fill table from file */
    for (idx=0; idx<MECHMAXROWS && fgets(line,1024,fp); )
      {
      if (*line != '#')
        {
        char *e, *p=line;
        int d;
        for (val=0; val<vals+1; val++)
          {
          d = strtod(p, &e);
          if (e == p)
            break;
          else
            {
            p = e;
            res[val + (vals+1)*idx] = d;
            }
          }
        if (res[(vals+1)*idx] >= 0)
	  idx++;
        }
      }
    fclose(fp);
    if (idx >= MECHMAXROWS)
      {
      printk("Failed in mech table init, too many rows in %s, fix MECHMAXROWS die.\n",tablepath);
      return(MECHINIT_ROWS);
      }
    }
return(MECHINIT_OK);
}

/* set mechanism state lookup table info */
/* returns 0 for OK or other status if program should terminate */

#define HMI_MECH_TABLES 3

int set_HMI_mech_values(DRMS_Record_t *rec)
{
static int called = 0;

static int pol[MECHMAXROWS*4];
static int tuning[MECHMAXROWS*5];
static int focus[MECHMAXROWS*3];

/* DRMS keyword names for mech table values - "short" names will be used. */
static char *pol_keys[] =     {"HPL1POS", "HPL2POS", "HPL3POS"};
static char *tuning_keys[] = {"HWL1POS", "HWL2POS", "HWL3POS", "HWL4POS"};
static char *focus_keys[] =  {"HCF1POS", "HCF2POS"};

static char *camkey = "HCAMID";

static MECH_TABINFO_t  tabinfo[] = {
	{POL_TABLE, "HPLTID", pol_keys, pol, 3},
	{TUNE_TABLE, "HWLTID", tuning_keys, tuning, 4},
	{FOCUS_TABLE, "HCFTID", focus_keys, focus, 2}
};
int tab;

/* on first time only, fetch and read each table. */
if (!called)
  {
  int status = 0;
  called = 1;
  
  /* init each table */
  if (init_HMI_mech_tables(tabinfo, 4) != 0)
    {
    printk("Failed it initialize mechanism tables, code=%d\n",status);
    return(status);
    }

  }

for (tab=0; tab < HMI_MECH_TABLES; tab++)
  {
  int retried = 0;
  int found_index;
  int row, index;
  int status;
  int val, vals = tabinfo[tab].cols;
  int *res = tabinfo[tab].table;
  char **keys = tabinfo[tab].keys;
  found_index = 0;
  index = drms_getkey_int(rec, tabinfo[tab].index, &status);
  if (status)
    {
    printk("Mech Index %s not found.\n",tabinfo[tab].index);
    return(MECH_KEY_MISSING);
    }
  /* search for row index in table */
  for (row=0; row<=MECHMAXROWS; row++)
    {
    if (row == MECHMAXROWS)
      { /* index not found, reset table once else giveup */
      if (!retried)
	{
        if (init_HMI_mech_tables(tabinfo, 4) != 0)
          {
          printk("Failed it initialize mechanism tables, code=%d\n",status);
          return(status);
          }
        row = 0; /* start over */
        retried = 1;
        }
      else
        {
        printk("Failed to find given index %d in mech table %d\n",index,tab);
        return(MECH_INDEX_MISSING);
        }
      }
    if (index == res[(vals+1)*row])
      { /* found proper row for this image */
      for (val=0; val<vals; val++)
	  {
          drms_setkey_int(rec, keys[val], res[val+1+(vals+1)*row]);
	  }
      break;
      }
    }
  }
return(0);
}

