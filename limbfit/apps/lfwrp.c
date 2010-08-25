/* 
	I.Scholl 

		High level wrapper for the limbfit code
		Can be used as a stand-alone module (includes a 'DoIt' function)
			calls: do_one_limbfit() that calls:
										limbfit() and reformat output data to be rewritten in the DB
			do_one_limbfit() can be plugged into lev0lev1 module
			
	
	#define CODE_NAME 		"limbfit"
	#define CODE_VERSION 	"V1.4r0" 
	#define CODE_DATE 		"Tue Aug 24 17:46:37 PDT 2010" 
*/

#include "limbfit.h"

char *module_name = "lfwrp";

ModuleArgs_t module_args[] = {
  {ARG_STRING, 	"dsin", "hmi.lev1c_nrt[]", "input data set"},
  {ARG_STRING, 	"dsout", "su_scholl.limbfit", "output data set"},
  {ARG_STRING, 	"logdir", "./", "logs directory"},
  {ARG_STRING, 	"tmpdir", "./", "tmp directory"},
  {ARG_INTS, 	"bfsn", "0", "first lev1 fsn# to process"},
  {ARG_INTS, 	"efsn", "0", "last  lev1 fsn# to process."},
  {ARG_INT, 	"debug", "0", "debug level (default: 0)"},
  {ARG_END}

};

void get_sdate(char *sdate)
{

  time_t t = time(NULL);
  struct tm *timeptr;

  timeptr = localtime(&t);
  sprintf(sdate, "%d.%02d.%02d_%02d:%02d:%02d", 
	  (timeptr->tm_year+1900), (timeptr->tm_mon+1),
	  timeptr->tm_mday, timeptr->tm_hour, timeptr->tm_min, timeptr->tm_sec);
}


void lf_logmsg(char * type1, char * type2, int return_code, int status, char *message, char *code_name, FILE *opf)
{
	static char sdate[32];
	get_sdate(sdate);
	if ( !(strcmp(type1,"INFO")) || (status==0 &&return_code==0))
		fprintf(opf,"%s/%20s: %s %5s %5s msg: '%s' \n",LOGMSG1, code_name, sdate,type1,type2, message);
	else
		fprintf(opf,"%s/%20s: %s %5s %5s msg: '%s' return code: %d exit: %d\n",LOGMSG1, code_name, sdate,type1,type2, message, return_code,status);
}

void close_on_error(DRMS_Record_t *record_in,DRMS_Record_t *record_out,DRMS_Array_t *data_array) //,FILE *opf)
{
	drms_free_array (data_array);
	drms_close_record (record_out, DRMS_FREE_RECORD);
	drms_close_record (record_in, DRMS_FREE_RECORD);
	//fclose(opf);
}

//************************************************************************
// Get N records, decide or not to process them
//************************************************************************

int process_n_records(char * open_dsname, char *dsout, char *tmp_dir, FILE *opf, int debug, int *status)    
{
	static char *log_msg_code="process_n_records";
	char log_msg[120];
	sprintf(log_msg,"doing process for %s",open_dsname);
	lf_logmsg("INFO", "APP", 0,0, log_msg, log_msg_code, opf);			
	static char errcode[20]=PROCSTAT_NOK;

	//************************************************************************
	//  Open DRMS connexion, get the range of data, decide what to do, 
	//			and either call the next step or skip
	//************************************************************************
	
	DRMS_RecordSet_t *drs_in,*drs_out;
	DRMS_Record_t *record_in,*record_out;
	int rstatus, ncnt,r;

    drs_in = drms_open_records(drms_env, open_dsname, &rstatus);
	if (!drs_in) {
		sprintf(log_msg,"unable to open record set %s",open_dsname);
		lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ, rstatus, log_msg, log_msg_code, opf);			
		fprintf (stderr, log_msg);
		*status=ERR_DRMS_READ;
		return(ERR_DRMS_READ);
	}

	ncnt = drs_in->n;
	if (ncnt < 1) {
		sprintf(log_msg,"no records in selected set %s",open_dsname);
		lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ, rstatus, log_msg, log_msg_code, opf);			
		fprintf (stderr, log_msg);
		*status=ERR_DRMS_READ;
		return(ERR_DRMS_READ);
	}

    drs_out = drms_create_records(drms_env, ncnt, dsout, DRMS_PERMANENT,&rstatus);
	if (!drs_out) {
		sprintf(log_msg,"unable to create record set %s",dsout);
		lf_logmsg("ERROR", "DRMS", ERR_DRMS_WRITE, rstatus, log_msg, log_msg_code, opf);			
		fprintf (stderr, log_msg);
		*status=ERR_DRMS_WRITE;
		return(ERR_DRMS_WRITE);
	}
    
	unsigned int fsn = 0;
	char *imgtype;
    char *hwltnset;
	int missvals;
	
    for(r=0; r < ncnt; r++) 
    {
		record_in = drs_in->records[r];
		record_out = drs_out->records[r]; 

		fsn = drms_getkey_int(record_in, "FSN", &rstatus);       
		if(rstatus) {
			lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_KW, rstatus, "drms_getkey_string(FSN)", log_msg_code, opf);			
			write_mini_output(PROCSTAT_NO_LF_DB_READ_PB,record_in,record_out,opf,debug); 
			*status=ERR_DRMS_READ_MISSING_KW;   
			return(0);   
		}
		sprintf(log_msg,"FSN:%u processing",fsn);
		lf_logmsg("INFO", "APP", 0,0, log_msg, log_msg_code, opf);			

		hwltnset = drms_getkey_string(record_in, "HWLTNSET", &rstatus);
		if(rstatus) {
			lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_KW, rstatus, "drms_getkey_string(HWLTNSET)", log_msg_code, opf);			
			write_mini_output(PROCSTAT_NO_LF_DB_READ_PB,record_in,record_out,opf,debug); 
			*status=ERR_DRMS_READ_MISSING_KW;   
			return(0);   
		}
		missvals = drms_getkey_int(record_in, "MISSVALS", &rstatus);
		if(rstatus) {
			lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_KW, rstatus, "drms_getkey_int(MISSVALS)", log_msg_code, opf);			
			write_mini_output(PROCSTAT_NO_LF_DB_READ_PB,record_in,record_out,opf,debug); 
			*status=ERR_DRMS_READ_MISSING_KW;   
			return(0);   
		}
		imgtype = drms_getkey_string(record_in, "IMG_TYPE", &rstatus);
		if(rstatus) {
			lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_KW, rstatus, "drms_getkey_string(IMG_TYPE)", log_msg_code, opf);			
			write_mini_output(PROCSTAT_NO_LF_DB_READ_PB,record_in,record_out,opf,debug); 
			*status=ERR_DRMS_READ_MISSING_KW;   
			return(ERR_DRMS_READ_MISSING_KW);   
		}

		sprintf(log_msg," Selection: %s %d %s", hwltnset,missvals,imgtype);
		lf_logmsg("INFO", "APP", 0, 0, log_msg, log_msg_code, opf);			
		if (debug) printf("selection: %s %d %s\n",hwltnset,missvals,imgtype);

   		if ((strncmp(hwltnset,"CLOSE",5)==0) && (missvals==0) && (strncmp(imgtype,"LIGHT",5)==0)) 
		{
			if (do_one_limbfit(fsn,record_in,record_out,tmp_dir,opf,debug,&rstatus))
			{
				if (rstatus < 0 && rstatus > -300)
				{
					drms_close_records(drs_out, DRMS_INSERT_RECORD);
					drms_close_records(drs_in, DRMS_FREE_RECORD);
					lf_logmsg("ERROR", "APP", rstatus, 0, "to be aborted", log_msg_code, opf);			
					return(rstatus);
				}
			}
		}
		else
		{
			sprintf(errcode,"%s",PROCSTAT_NOK);
			if (strncmp(hwltnset,"CLOSE",5)!=0) sprintf(errcode,"%s",PROCSTAT_NO_LF_OPENLOOP); 
				else if (missvals!=0) sprintf(errcode,"%s",PROCSTAT_NO_LF_MISSVALS);
					else if (strncmp(imgtype,"LIGHT",5)!=0) sprintf(errcode,"%s",PROCSTAT_NO_LF_DARKIMG);
			write_mini_output(errcode,record_in,record_out,opf,debug);
		}
	}
	drms_close_records(drs_out, DRMS_INSERT_RECORD);
	drms_close_records(drs_in, DRMS_FREE_RECORD);
	lf_logmsg("INFO", "APP", 0, 0, "Records saved", log_msg_code, opf);			
 	if (debug) printf("records saved\n");
	
return(0);
}


//************************************************************************
// MAIN
// 			usage ./lfwrp bfsn efsn [dsin] [dsout] [debug]
//************************************************************************

int DoIt(void)
{

	CmdParams_t *params = &cmdparams;
	
	int  debug = params_get_int (params, "debug");
	char* dsin = params_get_str (params, "dsin");
	char* dsout = params_get_str (params, "dsout");
	char* log_dir = params_get_str (params, "logdir");
	char* tmp_dir = params_get_str (params, "tmpdir");
	long long bfsn = params_get_int (params, "bfsn");
	long long efsn = params_get_int (params, "efsn");

	static char *log_msg_code="DoIt";
	char log_msg[120];
	int result;
	
	static char open_dsname[256];
	char recrange[128];

    if(bfsn == 0 || efsn == 0) 
    {
		fprintf(stderr, "bfsn and efsn must be given for fsn mode. 0 not allowed\n");
		return(0);
    }
    if(bfsn > efsn) 
    {
		fprintf(stderr, "bfsn must be <= efsn\n");
		return(0);
    }

	//**********************************************************************
	//	LOG FILE: MAKE A GLOBAL LOG FILE FOR THE WHOLE SESSION: 
	//**********************************************************************
	static char flogname[128];
	static char sdate[32];
	get_sdate(sdate);

	FILE *opf;
	sprintf(flogname, "%slimbfit_%s_%lld_%lld.log",log_dir,sdate,bfsn,efsn);
    if((opf=fopen(flogname, "w")) == NULL)
    {
		fprintf(stderr, "**Can't open the log file %s\n", flogname);
		return(0);
	}

	lf_logmsg("INFO", "APP", 0, 0, "Begin... ", log_msg_code, opf);

	//------------------------------------------------------
	// process NUMRECLEV1 at a time during one transaction
	//------------------------------------------------------
	long long numofrecs, frec, lrec;
	int numrec, numofchunks, i;    
	numofrecs = (efsn - bfsn) + 1;
	numrec = NUMRECLEV1;
	numofchunks = numofrecs/numrec;
	if((numofrecs % numrec) != 0) numofchunks++; //extra loop for partial chunk
	lrec = bfsn-1;
	for(i = 0; i < numofchunks; i++) 
	{
		frec = lrec+1; 
		lrec = (frec + numrec)-1;
		if(lrec > efsn) lrec=efsn;
		sprintf(recrange, "%lld-%lld", frec, lrec);
		sprintf(open_dsname, "%s[%s]", dsin, recrange);	
		sprintf(log_msg,"open %s", open_dsname);
		lf_logmsg("INFO", "APP", 0, 0, log_msg, log_msg_code, opf);
		if(process_n_records(open_dsname, dsout, tmp_dir,opf, debug, &result)) 
		{  //do a chunk to get files from the lev0
			if (result < 0 && result > -300)
			{
				lf_logmsg("ERROR", "ABORT", result, 0, "", log_msg_code, opf);
				fprintf(opf,"lfwrp abort\nSee log: %s\n", flogname); 
				//send_mail("build_lev1 abort\nSee log: %s\n", logname); 
				fclose(opf);
				return(0);
			}
		}
		drms_server_end_transaction(drms_env,0,0);
		drms_server_begin_transaction(drms_env);
	}
	fclose(opf);
	
	lf_logmsg("INFO", "APP", 0, 0, "End... ", log_msg_code, opf);

return(0);
} //end doit()
    
