/* 
	I.Scholl 

		High level wrapper for the limbfit code
		Can be used as a stand-alone module (includes a 'DoIt' function)
			calls: do_one_limbfit() that calls:
										limbfit() and reformat output data to be rewritten in the DB
			do_one_limbfit() can be plugged into lev0lev1 module
			
	
	#define CODE_NAME 		"limbfit"
	#define CODE_VERSION 	"V2.0r0" 
	#define CODE_DATE 		"Fri Mar  2 13:28:38 PST 2012" 
*/

#include "limbfit.h"

char *module_name = "lfwrp";

ModuleArgs_t module_args[] = {
  {ARG_STRING, 	"dsin", "hmi.lev1", "input data set"},
  {ARG_STRING, 	"dsout", "su_scholl.limbfit", "output data set"},
  {ARG_STRING, 	"logdir", "./", "logs directory"},
  {ARG_STRING, 	"tmpdir", "./", "tmp directory"},
  {ARG_INTS, 	"bfsn", "0", "first lev1 fsn# to process"},
  {ARG_INTS, 	"efsn", "0", "last  lev1 fsn# to process"},
  {ARG_INT, 	"cam", "0", "camera selection: 0: both = default, otherwise: 1 or 2"},
  {ARG_INT, 	"fid", "0", "default=all, otherwise one FID number or mask"},
//  {ARG_INT, 	"smpl", " ", "sampling factor: default=no sampling, otherwise: number/unit eg: 10 first of the hour: 10/h. this is combind with cam and fid"},
  {ARG_INT, 	"spe", "0", "used to activate only 1 iteration of the FORTRAN code and AHI low for roll analaysis (default=0 normal processing - =1 activate)"},
  {ARG_INT, 	"cc",  "0", "used to activate center calculation instead of using X0/YO_LF (default=0 no center calculation, =1 calculation)"},
  {ARG_INT, 	"iter", "2", "to change the number of iterations of the limb.f code (default=2)"},
  {ARG_INT, 	"fldf", "0", "0: no full ldf processing (default), 1: with full ldf proc."},
  {ARG_INT, 	"src", "0", "0: hmi.lev1 (default), 1: hmi.lev1_nrt  (quality test is different)"},
//  {ARG_INT, 	"sav", "0", "0: do not save when the process ends with failure, otherwise: 1"},
  {ARG_STRING, 	"comment", " ", "to add a comment in a dataset"},
  {ARG_STRING, 	"preset", "none", "current possible values: roll_cont, roll_all, pipeline, (default=none)"},
  {ARG_INT, 	"debug", "0", "debug level (default: 0 no debug info)"},
  {ARG_END}

};  // ADD A SWITCH FOR A DIFFERENT QUALITY KW

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

int process_n_records(char * open_dsname, char *dsout, char *tmp_dir, FILE *opf, int cc, int spe, int iter, int fldf, char *dsin, char *comment, char *bld_vers, int debug, int *status)    
{
	static char *log_msg_code="process_n_records";
	char log_msg[200];
	sprintf(log_msg,"doing process for %s -> %s",open_dsname,dsout);
	lf_logmsg("INFO", "APP", 0,0, log_msg, log_msg_code, opf);			
	static char errcode[20]=PROCSTAT_NOK;

	// for error management purposes
	static LIMBFIT_OUTPUT limbfit_res_t ;
	static LIMBFIT_OUTPUT *lfr_t = &limbfit_res_t;
	lfr_t->code_name=CODE_NAME;
	lfr_t->code_version=CODE_VERSION;
	lfr_t->code_date=CODE_DATE;
	lfr_t->comment=comment;
	lfr_t->dsin=dsin;
	lfr_t->bld_vers=bld_vers;
	
	//************************************************************************
	//  Open DRMS connexion, get the range of data, decide what to do, 
	//			and either call the next step or skip
	//************************************************************************
	
	DRMS_RecordSet_t *drs_in,*drs_out;
	DRMS_Record_t *record_in,*record_out;
	int rstatus, ncnt,r;

    drs_in = drms_open_records(drms_env, open_dsname, &rstatus);
	if (!drs_in) {
		sprintf(log_msg,"unable to open record set %s\n",open_dsname);
		lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ, rstatus, log_msg, log_msg_code, opf);			
		fprintf (stderr, log_msg);
		*status=ERR_DRMS_READ;
		return(ERR_DRMS_READ);
	}

	ncnt = drs_in->n;
	if (ncnt < 1) {
		sprintf(log_msg,"no records in selected set %s\n",open_dsname);
		lf_logmsg("WARNING", "DRMS", WAR_DRMS_NORECORD, rstatus, log_msg, log_msg_code, opf);			
		fprintf (stderr, log_msg);
		*status=WAR_DRMS_NORECORD;
		drms_close_records(drs_in, DRMS_FREE_RECORD);
		return(WAR_DRMS_NORECORD);
	}

    drs_out = drms_create_records(drms_env, ncnt, dsout, DRMS_PERMANENT,&rstatus);
	if (!drs_out) {
		sprintf(log_msg,"unable to create record set %s\n",dsout);
		lf_logmsg("ERROR", "DRMS", ERR_DRMS_WRITE, rstatus, log_msg, log_msg_code, opf);			
		fprintf (stderr, log_msg);
		*status=ERR_DRMS_WRITE;
		drms_close_records(drs_in, DRMS_FREE_RECORD);
		return(ERR_DRMS_WRITE);
	}
    
	unsigned int fsn = 0;
	
    for(r=0; r < ncnt; r++) 
    {
		record_in = drs_in->records[r];
		record_out = drs_out->records[r]; 

		fsn = drms_getkey_int(record_in, "FSN", &rstatus);       
		if(rstatus) {
			lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_KW, rstatus, "drms_getkey_string(FSN)", log_msg_code, opf);			
			write_mini_output(PROCSTAT_NO_LF_DB_READ_PB,record_in,record_out,opf,0,lfr_t,debug);
			*status=ERR_DRMS_READ_MISSING_KW;   
			drms_close_records(drs_out, DRMS_INSERT_RECORD);
			drms_close_records(drs_in, DRMS_FREE_RECORD);
			return(0);   
		}
		sprintf(log_msg,"FSN:%u processing",fsn);
		lf_logmsg("INFO", "APP", 0,0, log_msg, log_msg_code, opf);			
		printf("selection: #%u\n",fsn);

		if (do_one_limbfit(fsn,record_in,record_out,tmp_dir,opf,cc,spe,iter,fldf,dsin,comment,bld_vers,debug,&rstatus))
		{
			if (rstatus < 0 && rstatus > -300)
			{
				drms_close_records(drs_out, DRMS_INSERT_RECORD);
				drms_close_records(drs_in, DRMS_FREE_RECORD);
				lf_logmsg("ERROR", "APP", rstatus, 0, "to be aborted", log_msg_code, opf);			
				return(rstatus);
			}
		}
		fflush(opf);
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
	
	int  debug 	= params_get_int (params, "debug");
	int  spe 	= params_get_int (params, "spe");
	int  cc 	= params_get_int (params, "cc");
	int  fldf 	= params_get_int (params, "fldf");
	int  iter	= params_get_int (params, "iter");
	int  cam 	= params_get_int (params, "cam");
	int  fid 	= params_get_int (params, "fid");
	int  src 	= params_get_int (params, "src");
//	int  sav 	= params_get_int (params, "sav");
	char* dsin 	= params_get_str (params, "dsin");
	char* dsout = params_get_str (params, "dsout");
	char* log_dir = params_get_str (params, "logdir");
	char* tmp_dir = params_get_str (params, "tmpdir");
	char* preset = params_get_str (params, "preset");
	char* comment = params_get_str (params, "comment");
	long long bfsn = params_get_int (params, "bfsn");
	long long efsn = params_get_int (params, "efsn");

	static char *log_msg_code="DoIt";
	char log_msg[200];
	int result;
		
	static char open_dsname[400];
	char recrange[128];
	static char qual[15];

	if (bfsn == 0 || efsn == 0) 
    {
		fprintf(stderr, "bfsn and efsn must be given for fsn mode. \n");
		return(0);
    }
	if(bfsn > efsn) 
    {
		fprintf(stderr, "bfsn must be <= efsn\n");
		return(0);
    }
    if(cam < 0 || cam > 2) 
    {
		fprintf(stderr, "cam must be equal to 0, 1, or 2\n");
		return(0);
    }
	if(src == 0) 
		sprintf(qual,"0");
	else 
		sprintf(qual,"1073741824"); //= 0x40000000 
		
	char ext_dsin[256];
//	if !(strncmp(dates,"-",1) 
//	{
//		sprintf(ext_dsin, "%s[%s]", dsin,dates); 
//	}
//	else
//	{
		sprintf(ext_dsin, "%s[]", dsin);
//	}	

	static char bld_vers[16];
  	sprintf(bld_vers, "%s", jsoc_version);
 	
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
	//prepa base name
	static char tcam[20];
	static char tfid[60];
	static char tfil[20];
	static char tbase[128];
	sprintf(tfil, "[? quality = %s ", qual);
		
	if(strcmp(preset,"pipeline"))
	{
		sprintf(tcam, " ");
		sprintf(tfid, "AND (FID<10010 or CAST(FID AS TEXT) LIKE '%8')");
	}
	else 
	{
		if(strcmp(preset,"roll_all"))
		{
			sprintf(tcam, " ");
			sprintf(tfid, " ");
		}
		else
		{
			if(strcmp(preset,"roll_cont"))
			{
				sprintf(tcam, " ");
				sprintf(tfid, "AND FID<10010 ");
			}
			else /*if(strcmp(preset,"none") */
			{
				if (cam == 1 || cam == 2) sprintf(tcam, " and camera = %d ", cam);	else sprintf(tcam, " ");
				if (fid != 0) sprintf(tfid, " and fid = %d ", fid);	else sprintf(tfid, " ");
			}
		}
	}
	sprintf(tbase, "%s%s%s ?]", tfil,tcam,tfid);	
printf("%s\n",tbase);
	for(i = 0; i < numofchunks; i++) 
	{
		frec = lrec+1; 
		lrec = (frec + numrec)-1;
		if(lrec > efsn) lrec=efsn;
		sprintf(recrange, "%lld-%lld", frec, lrec);
		sprintf(open_dsname, "%s[%s]%s", ext_dsin, recrange,tbase);	
		
		sprintf(log_msg,"open %s -> %s (logdir: %s, tmpdir: %s, cam: %d, spe: %d, cc: %d, iter: %d, comment: %s , debug: %d)",
								open_dsname,dsout,log_dir,tmp_dir,cam,spe,cc,iter,comment,debug);

		lf_logmsg("INFO", "APP", 0, 0, log_msg, log_msg_code, opf);
		if(process_n_records(open_dsname, dsout, tmp_dir,opf, cc, spe, iter, fldf, dsin, comment, bld_vers, debug, &result)) 
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
		sprintf(log_msg,"close %s", open_dsname);
		lf_logmsg("INFO", "APP", 0, 0, log_msg, log_msg_code, opf);
		drms_server_end_transaction(drms_env,0,0);
		drms_server_begin_transaction(drms_env);
	}
	lf_logmsg("INFO", "APP", 0, 0, "End... ", log_msg_code, opf);
	fclose(opf);

return(0);
} //end doit()
    
