/*
 *  ingest_mharp_log.c
 *  
 *
 *  Purpose:
 *    Ingest HARP log file, error file, and checkpoint file.  Record
 *    per-run information, such as the first and last image records 
 *    processed by the tracker.
 *
 *  Input:
 *    Files are taken from specified locations underneath a standard
 *    root directory, and ingested directly:
 *      log file (text):          root/track-latest.log
 *	error file (text):        root/track-latest.err
 *      checkpoint file (binary): root/track-post.mat
 *    Metadata is read from a key: value file, root/track-param.txt
 *
 *  Output:
 *    Files and metadata are inserted into a tailor-made data series.
 *
 *  Usage:
 *
 *  ingest_mharp_log root=/tmp22/HARPs/Tracks/jsoc out=hmi.Mharp_log_720s \
 *       "trec=2011.03.07_12:00_TAI" 
 *
 *
 *  Written by X. Sun, M. Turmon
 *
 *	History:
 *	v0.0	Jun 23 2011
 *	v1.0	Aug 2  2011
 *	v1.1	Feb 13 2012
 *	v1.2	Mar 8  2012
 *
 *	Notes:
 *
 *	v0.0
 *	v1.0 -- more error checking, default file args, name change
 *      v1.1 -- ingest error log as well
 *      v1.2 -- more metadata read from track-params file
 *
 */

#include "jsoc_main.h"
#include "drms_types.h"
#include <time.h>
#include <math.h>
#include <timeio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

// exit with error
//  (standard trick to swallow the semicolon)
//  (ensure nonzero exit even if zero status)
//  (for use *only* in DoIt() routine)
#define DIE(msg) do { \
	fflush(stdout); \
        fprintf(stderr, "%s: FATAL: %s. (status=%d)\n", module_name, msg, status); \
        return (status ? status : 1); \
        } while (0)
// report non-fatal error
//  (standard trick to swallow the semicolon)
#define WARN(msg) do { \
	fflush(stdout); \
        fprintf(stderr, "%s: WARNING: %s. Continuing.\n", module_name, msg); \
	fflush(stderr); \
        } while (0)
// V_printf: facilitate verbflag output
//   if flag is > 0, output is to stdout, if < 0, to stderr
//   if flag is 0, no output is made at all
// The message is printed in the form:
// <module_name>: <first><message>"
// Usage:
// V_printf(VERB > 0, "\t", "Mask(%d) = %d\n", 2048, mask[2048]);
void
V_printf(int flag, char *first, char *format, ...) {
  va_list args;
  extern char *module_name;
  FILE *fp = (flag > 0) ? stdout : stderr;

  va_start(args, format);
  if (flag != 0) {
    // first is a string, even "" -- print the module name too
    // otherwise, omit it
    if (first)
      fprintf(fp, "%s: %s", module_name, first);
    vfprintf(fp, format, args);
    fflush(fp);
  }
  va_end(args);
}

#define kRootDir	"root"
#define kRecSetOut	"out"
#define kTRec		"trec"
#define kLogFile	"log"
#define kErrFile	"err"
#define kCkpFile	"checkpoint"
#define kVerbosity      "verb"
#define kNOT_SPEC	"Not Specified"

#define STR_MAX 512
// per-run parameter info
#define FN_TRACK_PARAM "track-param.txt"

/*
 * Whole-run info
 * (run start/end, tracker parameters)
 */
typedef struct {
  int n_par;              // number of parameters
  char **par_names;       // parameter names
  char **par_vals;        // parameter values
} run_info_t;


char *module_name = "ingest_mharp_log";
static int verbflag;

// maps per-run parameter namess in run_info into key names
static const char *tracker_run_info_keys[] = {
  // per-run parameters -- exported from matlab
  "TKP_RUNN", "run_name",
  "TKP_RUNC", "run_number",
  "TKP_RUNT", "run_time",
  "T_FRST",   "first_frame_trec",
  "T_LAST",   "last_frame_trec",
  // end sentinel -- required
  NULL, NULL,
};


/*
 * load tracker parameters
 */
int
load_tracker_params(const char *rootDir, run_info_t *ri)
{
  FILE *fp;
  int p, np, llen;
  char fn[STR_MAX];
  char line[STR_MAX];
  char *val, *name_end;

  snprintf(fn, sizeof(fn), "%s/%s", rootDir, FN_TRACK_PARAM);
  if ((fp = fopen(fn, "r")) == NULL) {
    V_printf(-1, "", "Failed to open `%s'.\n", fn);
    return 1;
  }
  // find the number of parameters to be stored
  np = 0;
  while (fgets(line, sizeof(line), fp) != NULL) {
    llen = strlen(line);
    // we must have read all the way to the newline
    if (line[llen-1] != '\n') {
      V_printf(-1, "", "Line: `%.20s ...' in param file `%s' is too long.\n", line, fn);
      return 1;
    }
    if (line[0] == '#' || line[0] == '\n') continue;
    np++;
  }
  ri->n_par = np;
  if (np == 0) {
    fclose(fp);
    return 0;
  }
  // allocate space
  ri->par_names = calloc(np, sizeof(*(ri->par_names)));
  ri->par_vals  = calloc(np, sizeof(*(ri->par_vals )));
  if (ri->par_names == NULL || ri->par_vals == NULL) {
    V_printf(-1, "", "Calloc fail in load_tracker_params for np = %d\n", np);
    fclose(fp);
    return 1;
  }
  // read params
  rewind(fp);
  p = 0;
  while (fgets(line, sizeof(line), fp) != NULL) {
    llen = strlen(line);
    if (line[0] == '#' || line[0] == '\n') continue;
    // we know it's newline terminated, so drop the newline
    line[--llen] = '\0';
    // drop any other trailing whitespace
    while (line[llen-1] == ' ' || line[llen-1] == '\t')
      line[--llen] = '\0';
    // find the :
    val = strchr(line, (int) ':');
    if (val == NULL) {
      V_printf(-1, "", "load_tracker_params: no colon in line `%s' in file `%s'\n", 
	       line, fn);
      fclose(fp);
      return 1;
    }
    // trim trailing whitespace before : (if any)
    for (name_end = val-1; *name_end == ' ' || *name_end == '\t'; name_end--)
      *name_end = '\0';
    *val++ = '\0'; // null the colon, thus terminating *name
    // null any following whitespace
    while (*val == ' ' || *val == '\t')
      *val++ == '\0';
    // copy params
    ri->par_names[p] = strdup(line);
    ri->par_vals [p] = strdup( val);
    p += 1; // next param
  }
  fclose(fp);
  V_printf(verbflag > 0, "", "Got %d tracker params from `%s'\n", np, fn);
  return 0; // OK
}

/*
 * set_keys_runinfo: set up drms keys using per-tracker-run keys
 * which are stored in runInfo
 */
static
int
set_keys_runinfo(DRMS_Record_t *rec, 
		 run_info_t *runInfo)
{
  const char **keys = tracker_run_info_keys;
  const char *drms_key;
  const char *tracker_key;
  int not_ok, iKey;    // outer loop
  int p, found, ok;    // inner loop

  for (not_ok = iKey = 0; keys[iKey] != NULL; iKey += 2) {
    // map drms_key -> tracker_key
    drms_key    = keys[iKey];
    tracker_key = keys[iKey+1];
    // find a parameter p matching tracker_key
    for (found = p = 0; p < runInfo->n_par; p++) {
      if (strcasecmp(tracker_key, runInfo->par_names[p]) == 0) {
	found = 1;
	ok = drms_setkey_string(rec, drms_key, runInfo->par_vals[p]);
	if (ok != DRMS_SUCCESS) not_ok++;
      }
    }
    if (!found) not_ok++; // could not find a tracker param for the key
  }
  return not_ok;
}

/*
 * Count log lines indicating errors, return count or (-1)
 *
 * file has line groups like:
 *   T_REC = 2011.03.04_12:24:00_TAI
 *   Reason = Bad quality 0x00000200.  Skipping the image at 2011.03.04_12:24:00_TAI
 */
static
int
count_errlog_times(char *fn)
{
  FILE *fp;
  char line[STR_MAX];
  int count;
  const char *lookfor = "T_REC = ";

  if ((fp = fopen(fn, "r")) == NULL)
    return -1;
  count = 0;
  while (!feof(fp)) {
    fgets(line, sizeof(line), fp);
    if (strncmp(line, lookfor, sizeof(lookfor)) == 0)
      count++;
  }
  fclose(fp);
  return count;
}


/*
 * Main routine DoIt() and arguments
 *
 */

ModuleArgs_t module_args[] =
{
  {ARG_STRING, kRootDir,   kNOT_SPEC, "Input root directory name (typ. ends in Tracks/jsoc)"},
  {ARG_STRING, kRecSetOut, kNOT_SPEC, "Output data series"},
  {ARG_STRING, kTRec,      kNOT_SPEC, "T_REC as primary key"},
  {ARG_STRING, kLogFile,   "track-latest.log", "Log file name (typ. .log extension)"},
  {ARG_STRING, kErrFile,   "track-latest.err", "Error file name (typ. .err extension)"},
  {ARG_STRING, kCkpFile,   "track-post.mat",   "Checkpoint file name (.mat extension)"},
  {ARG_INT,    kVerbosity, "1",                "Verbosity: 0=errs only; 1, 2 = more"},
  {ARG_END}
};

int DoIt(void)

{
  int status = DRMS_SUCCESS;
  const char *outQuery, *trec;
  const char *rootDir;
  const char *log;             // basename
  const char *errlog;          // basename
  const char *ckp;             // basename
  char logFile[STR_MAX];       // full filename
  char errlogFile[STR_MAX];    // full filename
  char ckpFile[STR_MAX];       // full filename
  int errlog_present;          // did we find the error log?
  int errlog_count;            // count of #fails in error log
  DRMS_Record_t *outRec;
  DRMS_Segment_t *outSegLog, *outSegErr, *outSegChkpt;
  TIME trec_num;
  char trec_buf[32];
  run_info_t runInfo;
  struct stat buf; // for stat call
	
  /* Command line parameters */
  rootDir = cmdparams_get_str(&cmdparams, kRootDir, NULL);
  if (strcmp(rootDir, kNOT_SPEC) == 0) DIE("Root directory must be specified\n");
  outQuery = cmdparams_get_str(&cmdparams, kRecSetOut, NULL);
  if (strcmp(outQuery, kNOT_SPEC) == 0) DIE("Output series must be specified\n");	
  trec = cmdparams_get_str(&cmdparams, kTRec, NULL);
  if (strcmp(trec, kNOT_SPEC) == 0) DIE("T_REC must be specified\n");	
  // optional args
  log      = cmdparams_get_str(&cmdparams, kLogFile,   NULL);
  errlog   = cmdparams_get_str(&cmdparams, kErrFile,   NULL);
  ckp      = cmdparams_get_str(&cmdparams, kCkpFile,   NULL);
  verbflag = cmdparams_get_int(&cmdparams, kVerbosity, NULL); // an int, not a flag
	
  /* Get and check filenames */
  if (snprintf(logFile, sizeof(logFile), "%s/%s", rootDir, log) >= sizeof(logFile))
    DIE("Log filename was too long.\n");
  if (snprintf(errlogFile, sizeof(errlogFile), "%s/%s", rootDir, errlog) >= sizeof(errlogFile))
    DIE("Error log filename was too long.\n");
  if (snprintf(ckpFile, sizeof(ckpFile), "%s/%s", rootDir, ckp) >= sizeof(ckpFile))
    DIE("Checkpoint filename was too long.\n");
  printf("%s: Log filename: `%s'\n",        module_name, logFile);
  printf("%s: Error log filename: `%s'\n",  module_name, errlogFile);
  printf("%s: Checkpoint filename: `%s'\n", module_name, ckpFile);
  if (stat(logFile, &buf) < 0 || !S_ISREG(buf.st_mode))
    DIE("Log file does not exist or is not a regular file.\n");
  if (stat(ckpFile, &buf) < 0 || !S_ISREG(buf.st_mode))
    DIE("Checkpoint file does not exist or is not a regular file.\n");
  printf("%s: Found log and checkpoint files\n", module_name);
  // error log: OK to be absent
  if (stat(errlogFile, &buf) < 0 || !S_ISREG(buf.st_mode)) {
    errlog_present = 0;
    printf("%s: Error log file does not exist: this is OK.\n", module_name);
  } else {
    errlog_present = 1;
    printf("%s: Found an error log file `%s', some T_RECs probably were skipped.\n", 
	   module_name, errlog);
  }
  // tracker parameters
  if (load_tracker_params(rootDir, &runInfo))
    WARN("Failed to load tracker parameter file.  Some keys will not be set");

  /* Get and check the time */
  trec_num = sscan_time((char *) trec);
  if (time_is_invalid(trec_num))
    DIE("Could not understand T_REC input time\n");
  sprint_at(trec_buf, trec_num); // round-trip it to TAI
  printf("%s: T_REC = %s\n", module_name, trec_buf);

  /* Open output record */
  outRec = drms_create_record(drms_env, (char *) outQuery, DRMS_PERMANENT, &status);
  if (status) {
    drms_close_record(outRec, DRMS_FREE_RECORD);
    DIE("Output record not found\n");
  }
	
  /* Ingest files */
  // regular log
  outSegLog = drms_segment_lookup(outRec, "Log");
  if (!outSegLog)
    DIE("Log segment lookup failed.\n");
  if (drms_segment_write_from_file(outSegLog, logFile))
    DIE("Log file ingestion failed.\n");
  printf("%s: ingested `%s'\n", module_name, logFile);

  // error log
  outSegLog = drms_segment_lookup(outRec, "ErrLog");
  if (!outSegLog)
    DIE("Error log segment lookup failed.\n");
  if (errlog_present) {
    if (drms_segment_write_from_file(outSegLog, errlogFile))
      DIE("Error log file ingestion failed (even though error log was present).\n");
    printf("%s: ingested `%s'\n", module_name, errlogFile);
    errlog_count = count_errlog_times(errlogFile);
  } else {
    if (drms_segment_write_from_file(outSegLog, "/dev/null"))
      DIE("Error log file ingestion failed (no error log present -> used /dev/null).\n");
    printf("%s: No error log, set errors to empty.\n", module_name);
    errlog_count = 0;
  }
  // checkpoint
  outSegChkpt = drms_segment_lookup(outRec, "Checkpoint");
  if (!outSegChkpt)
    DIE("Checkpoint segment lookup failed.\n");
  if (drms_segment_write_from_file(outSegChkpt, ckpFile))
    DIE("Checkpoint file ingestion failed.\n");
  printf("%s: ingested `%s'\n", module_name, ckpFile);

  /* Keywords */
  // errors here are not fatal, just count them
  int errs = 0;
  errs += (drms_setkey_time  (outRec, "T_REC",    trec_num)           != DRMS_SUCCESS);
  errs += (drms_setkey_string(outRec, "BLD_VERS", jsoc_version)       != DRMS_SUCCESS);
  errs += (drms_setkey_time  (outRec, "DATE",     CURRENT_SYSTEM_TIME)!= DRMS_SUCCESS);
  errs += (drms_setkey_int   (outRec, "NFAILED",  errlog_count)       != DRMS_SUCCESS);
  // per tracker run keys
  errs += set_keys_runinfo(outRec, &runInfo);
  printf("%s: finished keywords, %d keyword errors.\n", module_name, errs);
	
  /* Finish up */
  drms_close_record(outRec, DRMS_INSERT_RECORD);
  printf("%s: done.\n", module_name);
  return DRMS_SUCCESS;
}

