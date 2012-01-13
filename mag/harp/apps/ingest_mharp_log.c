/*
 *  ingest_mharp_log.c
 *  
 *
 *  Purpose:
 *	Ingesting HARP log file and checkpoint binary file
 *
 *	Input:
 *	  text log file (from a named file)
 *        binary checkpoint file (another named file)
 *
 *	Output:
 *	  ingest files into a tailor-made data series
 *
 *	Usage:
 *
 *  ingest_mharp_log "root=/tmp22/xudong/jsoc" "out=su_xudong.test_mharp_log" \
 *       "trec=2011.03.07_12:00_TAI" [ "run_name=2011.03.07_12:00_TAI/10d" ]
 *
 *  The "run_name" argument is optional, and gives a run name which can be helpful
 *  in determining which run ended in the given checkpoint.  It can be any
 *  string.
 *  
 *
 *  Written by X. Sun, M. Turmon
 *
 *	History:
 *	v0.0	Jun 23 2011
 *	v1.0	Aug 2  2011
 *
 *	Notes:
 *
 *	v0.0
 *	v1.0 -- more error checking, default file args, name change
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

#define DIE(msg)	return(fprintf(stderr, "FATAL: %s", msg), 1)
#define DIE_status(msg)	return(fprintf(stderr, "FATAL: %s, status=%d", msg, status), 1)
#define WARN(msg)	{printf("WARNING: %s", msg); fflush(stdout); continue;}

#define kRootDir	"root"
#define kRecSetOut	"out"
#define kTRec		"trec"
#define kLogFile	"log"
#define kCkpFile	"checkpoint"
#define kRunName        "run_name"
#define kNOT_SPEC	"Not Specified"

#define STR_MAX 240

char *module_name = "ingest_mharp_log";

ModuleArgs_t module_args[] =
{
  {ARG_STRING, kRootDir,   kNOT_SPEC, "Input root directory name (typ. ends in Tracks/jsoc)"},
  {ARG_STRING, kRecSetOut, kNOT_SPEC, "Output data series"},
  {ARG_STRING, kTRec,      kNOT_SPEC, "T_REC as primary key"},
  {ARG_STRING, kRunName,   "Un-named run",     "Indentifying name for run, optional string"},
  {ARG_STRING, kLogFile,   "track-latest.log", "Log file name (typ. .log extension)"},
  {ARG_STRING, kCkpFile,   "track-prior.mat",  "Checkpoint file name (.mat extension)"},
  {ARG_END}
};

int DoIt(void)

{
	int status = DRMS_SUCCESS;
	const char *outQuery, *trec, *run_name;
	const char *rootDir, *log, *ckp;             // portions of filenames
	char logFile[STR_MAX], ckpFile[STR_MAX];	     // full filenames
	DRMS_Record_t *outRec;
	DRMS_Segment_t *outSegLog, *outSegChkpt;
	TIME trec_num;
	char trec_buf[32];
	struct stat buf; // for stat call
	
	
	/* Command line parameters */
	rootDir = cmdparams_get_str(&cmdparams, kRootDir, NULL);
	if (strcmp(rootDir, kNOT_SPEC) == 0) DIE("Root directory must be specified\n");
	outQuery = cmdparams_get_str(&cmdparams, kRecSetOut, NULL);
	if (strcmp(outQuery, kNOT_SPEC) == 0) DIE("Output series must be specified\n");	
	trec = cmdparams_get_str(&cmdparams, kTRec, NULL);
	if (strcmp(trec, kNOT_SPEC) == 0) DIE("T_REC must be specified\n");	
	// optional
	log = cmdparams_get_str(&cmdparams, kLogFile, NULL);
	ckp = cmdparams_get_str(&cmdparams, kCkpFile, NULL);
	run_name = cmdparams_get_str(&cmdparams, kRunName, NULL);
	
	/* Get and check filenames */
	if (snprintf(logFile, sizeof(logFile), "%s/%s", rootDir, log) >= sizeof(logFile))
		DIE("Log filename was too long.\n");
	if (snprintf(ckpFile, sizeof(ckpFile), "%s/%s", rootDir, ckp) >= sizeof(ckpFile))
		DIE("Checkpoint filename was too long.\n");
	printf("%s: Log filename: `%s'\n", module_name, logFile);
	printf("%s: Checkpoint filename: `%s'\n", module_name, ckpFile);
        if (stat(logFile, &buf) < 0 || !S_ISREG(buf.st_mode))
		DIE("Log file does not exist or is not a regular file.\n");
        if (stat(ckpFile, &buf) < 0 || !S_ISREG(buf.st_mode))
		DIE("Checkpoint file does not exist or is not a regular file.\n");
	printf("%s: Found log and checkpoint files\n", module_name);


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
	outSegLog = drms_segment_lookup(outRec, "Log");
	if (drms_segment_write_from_file(outSegLog, logFile))
		DIE("Log file ingestion failed.\n");
	printf("%s: ingested `%s'\n", module_name, logFile);

	outSegChkpt = drms_segment_lookup(outRec, "Checkpoint");
	if (drms_segment_write_from_file(outSegChkpt, ckpFile))
		DIE("Checkpoint file ingestion failed.\n");
	printf("%s: ingested `%s'\n", module_name, ckpFile);

	/* Keywords */
	int errs = 0;
	errs += (drms_setkey_time  (outRec, "T_REC",    trec_num)           != DRMS_SUCCESS);
	errs += (drms_setkey_string(outRec, "BLD_VERS", jsoc_version)       != DRMS_SUCCESS);
	errs += (drms_setkey_string(outRec, "RUN_NAME", run_name)           != DRMS_SUCCESS);
	errs += (drms_setkey_time  (outRec, "DATE",     CURRENT_SYSTEM_TIME)!= DRMS_SUCCESS);
	printf("%s: finished keywords, %d keyword errors\n", module_name, errs);
	
	/* Finish up */
	drms_close_record(outRec, DRMS_INSERT_RECORD);
	printf("%s: done.\n", module_name);
	return(DRMS_SUCCESS);
}	// end of main
