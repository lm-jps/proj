#include "jsoc_main.h"
#include "drms.h"
#include "dsdsmigr.h"
// #include "printk.h"

/*
set_gaps_missing

Set_gaps_missing uses the same scanning method as show_coverage
to find contiguous sections of absent records in a slotted or integer
prime-key dataseries.  For the range specified set_gaps_missing will
generate records with all default values except if there is a QUALITY
keyword it will be set to 0x8000000 and if there is a DATAVALS it will
be set to 0.

QUALITY must be an int, else DATAVALS will be used.

Command line
  ds={seriesname}
  low=<first prime value to examine>
  high=<last prime value to examine>
  {key=<primekeyname>}
  {<other_primekey=value>}...
  -i  print index values vs prime key values in table
  -q  omit header information

The program operates as follows:

1.  get prime key list and determine name and type of index to use.
2.  get any other prime keys to use to filter the results.
3.  Read existing file or create new empty list.
4.  Scan target series from low to high limits and categorize each record
as OK, Missing, or unknown.
5.  Create records of MISSING for slots with no records matching the
    prime_keys specified.
5.  Write report coverage map, as header and table
first rows of header contain:
series=<seriesname>
key=<primekeyname> - user name, not slotted index name
type=<slotted type, or variable type, e.g. int or short or longlong>
step=<value of keyname_step>
epoch=<epoch>
low=<table first value
high=table high value.
end -- marke end of header

next print coverage table, 3 words per line
kind start length
where kind is:
    OK - data is present for this range
    MISS - data is all known to be permanently missing for this range.
    SETMISSING - data was absent but is now set to MISSING

start is in user units
length is count of recordscontiguous 

List terminates with line containing

*/

#define NOT_SPECIFIED "NOT_SPECIFIED"

#define DATA_OK ('\0')
#define DATA_MISS ('\1')
#define DATA_UNK ('\2')


char primestr[100];

#define RSUNM   (6.96e5)
#define AUM     (1.49597870e8)
#define SECRAD  (206264.8062)

int set_soho_ephemeris_keys (DRMS_Record_t *rec, TIME t) 
{
   TIME tbl_mod;
   double dist, lat, lon, vr, vn, vw, rsun;
   int cr, status;

   status = soho_ephemeris (t, &dist, &lat, &lon, &vr, &vn, &vw, &tbl_mod);
   if (!status) 
   {
      drms_setkey_double (rec, "CRLN_OBS", lon);
      drms_setkey_double (rec, "CRLT_OBS", lat);
      drms_setkey_double (rec, "OBS_DIST", dist);
      drms_setkey_double (rec, "OBS_VR", vr);
      drms_setkey_double (rec, "OBS_VN", vn);
      drms_setkey_double (rec, "OBS_VW", vw);
   }
   return status;
}

char *primevalstr(TIME prime, DRMS_Type_t type, char *unit, char *format)
  {
  if (type==DRMS_TYPE_TIME)
    sprint_time(primestr, prime, unit, atoi(format));
  else
    sprintf(primestr, (type <= DRMS_TYPE_LONGLONG ? "%.0f" : "%f"), prime);
  return(primestr);
  }

void printprime(FILE *fp, TIME prime, DRMS_Type_t type, char *unit, char *format)
  {
  fprintf(fp, primevalstr(prime, type, unit, format));
  }

ModuleArgs_t module_args[] =
{ 
    {ARG_STRING, "ds", NOT_SPECIFIED,  "Input data series."},
    {ARG_STRING, "low", NOT_SPECIFIED, "Low limit for coverage map."},
    {ARG_STRING, "high", NOT_SPECIFIED, "High limit for coverage map."},
    {ARG_STRING, "key", NOT_SPECIFIED, "Prime key name to use, default is first prime"},
    {ARG_FLAG, "i", "0", "Index - Print index values instead of prime slot values"},
    {ARG_FLAG, "q", "0", "Quiet - omit series header info"},
    {ARG_END}
};

#define DIE(msg) {fprintf(stderr,"%s\n",msg);exit(1);}


char *module_name = "show_coverage";

/* Module main function. */
int DoIt(void)
  {
  FILE *out;
  int status = 0;
  int slotted;
  long long lowslot, highslot, serieslowslot, serieshighslot;
  DRMS_RecordSet_t *rs;
  DRMS_Record_t *rec, *template;
  DRMS_Keyword_t *skey, *pkey;
  DRMS_Type_t ptype;
  char name[DRMS_MAXNAMELEN];
  int npkeys;
  char *pname;
  char *piname;
  char *punit;
  char *pformat;
  char *seriesname;
  TIME step, epoch;
  TIME series_low, series_high, low, high;
  char in[DRMS_MAXQUERYLEN];
  char *inbracket;
  char otherkeys[20*DRMS_MAXQUERYLEN];
  const char *ds = cmdparams_get_str (&cmdparams, "ds", NULL);
  const char *lowstr = cmdparams_get_str (&cmdparams, "low", NULL);
  const char *highstr = cmdparams_get_str (&cmdparams, "high", NULL);
  const char *skeyname = cmdparams_get_str (&cmdparams, "key", NULL);
  int quiet = cmdparams_get_int (&cmdparams, "q", NULL) != 0;
  int useindex = cmdparams_get_int (&cmdparams, "i", NULL) != 0;
  char *map;
  long long islot, jslot, nslots;
  char *qualkey;
  int qualkind;
  int ikey;
  int nOtherPrimes;
  struct OtherKeyStruct
    {
    char *name;
    char *value;
    } *OtherPrimes;

  /* check for minimum inputs */
  if (strcmp(ds, NOT_SPECIFIED) == 0 )
    DIE("No files: at least ds must be specified");
  out = stdout;


  /* get series info, low and high, and prime key type, etc. */
  strcpy(in,ds);
  inbracket = index(in, '[');
  if (inbracket)
	  *inbracket = '\0';
  else
	  inbracket = in + strlen(in);
  template = drms_template_record (drms_env, in, &status);
  if (!template || status)
	DIE("Series not found or empty");

  npkeys = template->seriesinfo->pidx_num;
  if (npkeys < 1)
    DIE("Series has no prime keys");
  if (strcmp(skeyname, NOT_SPECIFIED) != 0)
	{
	for (ikey=0; ikey < npkeys; ikey++)
		{
		pkey = template->seriesinfo->pidx_keywords[ikey];
		if (pkey->info->recscope > 1)
			skey = drms_keyword_slotfromindex(pkey);
		if (strcmp(skeyname, skey->info->name) == 0)
			break;
		}
	if (ikey == template->seriesinfo->pidx_num)
		DIE("name in key command line arg is not a prime key of this series");
	}
  else
	{
	skey = pkey = template->seriesinfo->pidx_keywords[0];
	}
  if (pkey->info->recscope > 1)
	{ // pkey is actual "_index" internal primekey, skey is user name for keyword
	skey = drms_keyword_slotfromindex(pkey); // Get the indexed keyword
	slotted = 1;
	}
  else
	slotted = 0;

  // now skey contains DRMS_Keyword_t for users prime key, pkey contains DRMS_Keyword_t for index
  ptype = skey->info->type;
  pname = strdup(skey->info->name);
  piname = strdup(pkey->info->name);
  punit = strdup(skey->info->unit);
  pformat = strdup(skey->info->format);
  seriesname = strdup(template->seriesinfo->seriesname);
  // get optional other primekeys
  otherkeys[0] = '\0';
  OtherPrimes = (struct OtherKeyStruct *)malloc(npkeys * sizeof(struct OtherKeyStruct));
  nOtherPrimes = 0;
  for (ikey=0; ikey < npkeys; ikey++)
	{
	DRMS_Keyword_t *tmppkey = template->seriesinfo->pidx_keywords[ikey];
	if (tmppkey->info->recscope > 1)
		tmppkey = drms_keyword_slotfromindex(pkey);
	if (cmdparams_exists(&cmdparams, tmppkey->info->name))
		{
                const char *value;
		char tmp[DRMS_MAXQUERYLEN];
		if (strcmp(tmppkey->info->name, pname) == 0)
			DIE("Can not have main prime key listed explicitly");
		value = cmdparams_get_str(&cmdparams, tmppkey->info->name, NULL);
		sprintf(tmp, "[%s=%s]", tmppkey->info->name, value);
		strcat(otherkeys, tmp);
		OtherPrimes[nOtherPrimes].name = strdup(tmppkey->info->name);
		OtherPrimes[nOtherPrimes].value = strdup(value);
                nOtherPrimes += 1;
		}
	}

  // get method to determine if all data in a record is missing
  DRMS_Keyword_t *qualitykeyword = drms_keyword_lookup(template, "QUALITY", 1);
  if (qualitykeyword && qualitykeyword->info->type == DRMS_TYPE_INT)
	{
	qualkey = "QUALITY";
	qualkind = 1;
	}
  else if (drms_keyword_lookup(template, "DATAVALS", 1))
	{
	qualkey = "DATAVALS";
	qualkind = 2;
	}
  else
	{
	qualkey = NULL;
	qualkind = 0;
	}
  // check prime key type for valid type for this program
  if (slotted == 0 && ( ptype != DRMS_TYPE_SHORT && ptype != DRMS_TYPE_INT && ptype != DRMS_TYPE_LONGLONG))
	DIE("Must be slotted or integer type first prime key");
  if (ptype == DRMS_TYPE_TIME)
	{
	strcpy(name, pname);
	strcat(name, "_epoch");
	epoch = drms_getkey_time(template, name, &status);
	strcpy(name, pname);
	strcat(name, "_step");
	step = drms_getkey_double(template, name, &status);
	}
  else if (slotted)
	{
	strcpy(name, pname);
	strcat(name, "_base");
	epoch = (TIME)drms_getkey_double(template, name, &status);
	strcpy(name, pname);
	strcat(name, "_step");
	step = (TIME)drms_getkey_double(template, name, &status);
	}
  else
	{
	epoch = (TIME)0.0;
	step = (TIME)1.0;
	}
  // Get series low info
  sprintf(in, "%s[%s=^]", seriesname, pname);
  // sprintf(in, "%s[%s=^]%s", seriesname, pname, otherkeys);
  rs = drms_open_records (drms_env, in, &status); // first record
  if (status || !rs || rs->n == 0)
	DIE("Series is empty");
  rec = rs->records[0];
  if (ptype == DRMS_TYPE_TIME)
	series_low = drms_getkey_time(rec, pname, &status);
  else if (slotted)
	series_low = (TIME)drms_getkey_double(rec, pname, &status);
  else
	series_low = (TIME)drms_getkey_longlong(rec, pname, &status);
  if (slotted)
	serieslowslot = drms_getkey_longlong(rec, piname, &status);
  else
	serieslowslot = series_low;
  drms_close_records(rs, DRMS_FREE_RECORD);
  if (strcmp(lowstr, NOT_SPECIFIED) == 0)
		low = series_low;
  else
	{
	if (ptype == DRMS_TYPE_TIME)
		low = sscan_time((char *)lowstr);
	else
		low = (TIME)atof(lowstr);
	}

  sprintf(in, "%s[%s=$]", seriesname, pname);
  // sprintf(in, "%s[%s=$]%s", seriesname, pname, otherkeys);
  rs = drms_open_records (drms_env, in, &status); // last record
  rec = rs->records[0];
  if (ptype == DRMS_TYPE_TIME)
	series_high = drms_getkey_time(rec, pname, &status);
  else if (slotted)
	series_high = (TIME)drms_getkey_double(rec, pname, &status);
  else
	series_high = (TIME)drms_getkey_longlong(rec, pname, &status);
  if (slotted)
	serieshighslot = drms_getkey_longlong(rec, piname, &status);
  else
	serieshighslot = series_high;
  drms_close_records(rs, DRMS_FREE_RECORD);
  if (strcmp(highstr, NOT_SPECIFIED) == 0)
	high = series_high;
  else
	{
	if (ptype == DRMS_TYPE_TIME)
		high = sscan_time((char *)highstr);
	else
		high = atof(highstr);
	}

  // Now get lowslot and highslot using same code as drms library calls.
  if (slotted)
	{
	DRMS_Value_t indexval;
	DRMS_Value_t inval;
	inval.value.double_val = low;
	inval.type = skey->info->type;
	drms_keyword_slotval2indexval(skey, &inval, &indexval, NULL);
	lowslot = indexval.value.longlong_val;

	inval.value.double_val = high;
	inval.type = pkey->info->type;
	drms_keyword_slotval2indexval(skey, &inval, &indexval, NULL);
	highslot = indexval.value.longlong_val;
	}
  else
	{
	lowslot = low;
	highslot = high;
	}

  // NOW get the record coverage info
  nslots = highslot - lowslot + 1;
  map = (char *)malloc(sizeof(char) * nslots);
  for (islot=0; islot<nslots; islot++)
	map[islot] = DATA_UNK;
  islot = 0;
  while (islot < nslots)
	{
	DRMS_Array_t *data;
	int nrecs, irec; 
	char query[DRMS_MAXQUERYLEN];
	char keylist[DRMS_MAXQUERYLEN];
	int qualindex=0;
	jslot = islot + 1000000;
	if (jslot >= nslots) jslot = nslots - 1;
	sprintf(query, "%s[%s=#%lld-#%lld]%s", seriesname, pname, lowslot+islot, lowslot+jslot,otherkeys);
	strcpy(keylist, piname);
	if (qualkind)
		{
		strcat(keylist, ",");
		strcat(keylist, qualkey);
		qualindex = 1;
		}
	data = drms_record_getvector(drms_env, query, keylist, DRMS_TYPE_LONGLONG, 0, &status);
	if (!data || status)
		{
		fprintf(stderr, "getkey_vector failed status=%d\n", status);
		DIE("getkey_vector failure");
		}
	nrecs = data->axis[1];
	for (irec = 0; irec < nrecs; irec++)
		{
		long long thisslot = *((long long *)data->data + irec);
		long long qualval;
		char val = DATA_OK;
		if (qualkind)
			{
			qualval = *((long long *)data->data + qualindex*nrecs + irec);
			if ((qualkind == 1 && qualval < 0) || (qualkind == 2 && qualval == 0))
			     val = DATA_MISS;
			}
		map[thisslot - lowslot] = val;
		}
	islot = jslot + 1;
	drms_free_array(data);
	}

  // now have low, high and series_low, series_high, epoch and step, and lowslot and highslot and serieshighslot and serieslowslot.

if (!quiet)
  {
  fprintf(out, "series=%s\n", seriesname);
  fprintf(out, "key=%s\n", pname);
  fprintf(out, "type=%s\n", drms_type2str(ptype)); // ptype as string 
  fprintf(out, "slotted=%s\n", (slotted ? "T" : "F"));
  fprintf(out, "epoch="); printprime(out, epoch, ptype, punit, pformat); fprintf(out, "\n");
  fprintf(out, "step=%f\n", step); // print step as proper format
  fprintf(out, "low="); printprime(out, low, ptype, punit, pformat); fprintf(out, "\n");
  fprintf(out, "high="); printprime(out, high, ptype, punit, pformat); fprintf(out, "\n");
  fprintf(out, "series_low="); printprime(out, series_low, ptype, punit, pformat); fprintf(out, "\n");
  fprintf(out, "series_high="); printprime(out, series_high, ptype, punit, pformat); fprintf(out, "\n");
  fprintf(out, "qualkey=%s\n", qualkey);
  }
// fprintf(out, "lowslot=%ld, highslot=%ld, serieslowslot=%ld serieshighslot=%ld\n",lowslot, highslot, serieslowslot, serieshighslot);

  islot = 0;
  while (islot < nslots)
	{
	long long startslot = islot;
	char pval[DRMS_MAXQUERYLEN];
	char primeval[DRMS_MAXQUERYLEN];
	int nsame = 1;
	if (useindex)
		 sprintf(pval, "%lld", lowslot + islot);
	else
		 sprintf(pval, "%s",
			 primevalstr(epoch + (lowslot + islot) * step, ptype, punit, pformat));
	char thisval = map[islot], nval = 0;
	for (islot += 1; islot < nslots && map[islot] == thisval; islot++)
	  nsame += 1;
        // Now have a contiguous set of one kind.  If the kind is UNK, create
        // a block of all missing records.
        if (thisval == DATA_UNK)
		{
		int irec;
                DRMS_RecordSet_t *rs = drms_create_records(drms_env, nsame, seriesname, DRMS_PERMANENT, &status);
                TIME primekeytime;
		if (status || !rs)
			DIE("Can not create records.");
fprintf(stderr,"setting %d records missing.\n",nsame);
		for (irec = 0; irec<nsame; irec++)
			{
			int iother;
			DRMS_Record_t *rec = rs->records[irec];
			sprintf(primeval, "%s",
                         primevalstr(epoch + (lowslot + startslot + irec) * step, ptype, punit, pformat));
			drms_setkey_string(rec, pname, primeval);
			for (iother=0; iother<nOtherPrimes; iother++)
				drms_setkey_string(rec, OtherPrimes[iother].name, OtherPrimes[iother].value);
			if (qualkind == 1)
				drms_setkey_int(rec, qualkey, 0X80000000);
			else if (qualkind == 2)
				drms_setkey_int(rec, qualkey, 0);
			if (qualitykeyword && qualitykeyword->info->type == DRMS_TYPE_STRING)
				drms_setkey_string(rec, "QUALITY", "0X80000000");

                        primekeytime = sscan_time(primeval);

                        set_soho_ephemeris_keys(rec, primekeytime);
			}
		if (drms_close_records(rs, DRMS_INSERT_RECORD))
			DIE("Failed at close new records\n");
		}
	fprintf(out, "%12s %s %d\n",
		(thisval == DATA_OK ? "OK" :
		(thisval == DATA_MISS ? "MISS" : "NOWMISSING")),
		pval, nsame );
	}
  return(0);
  }


