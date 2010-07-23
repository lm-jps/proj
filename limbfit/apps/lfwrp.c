/* 
	I.Scholl 

		High level wrapper for the limbfit code
		Can be used as a stand-alone module (includes a 'DoIt' function)
			calls: do_one_limbfit() that calls:
										limbfit() and reformat output data to be rewritten in the DB
			do_one_limbfit() can be plugged into lev0lev1 module
			
	
	#define CODE_NAME 		"limbfit"
	#define CODE_VERSION 	"V1.2r0" 
	#define CODE_DATE 		"Wed Jul 21 09:32:15 PDT 2010" 
*/

#include "limbfit.h"


ModuleArgs_t module_args[] = {
  {ARG_STRING, 	"dsin", "hmi.lev1c_nrt[]", "input data set"},
  {ARG_STRING, 	"dsout", "su_scholl.limbfit", "output data set"},
  {ARG_INTS, 	"bfsn", "0", "first lev1 fsn# to process"},
  {ARG_INTS, 	"efsn", "0", "last  lev1 fsn# to process."},
  {ARG_INT, 	"debug", "1", "debug level (default: 1)"},
  {ARG_END}

};

char *module_name = "lfwrp";


int process_n_records(char * open_dsname, char *dsout, FILE *opf, int debug)    
{
	fprintf(opf,"doing process for %s\n", open_dsname); 

	/***********************************************************************
		Open DRMS connexion, get the range of data, load structure and run lf
	************************************************************************/
	DRMS_RecordSet_t *drs,*drs_out;
	DRMS_Record_t *record_in,*record_out;
	DRMS_Segment_t *segment;
	int rstatus, ncnt,r;

    drs = drms_open_records(drms_env, open_dsname, &rstatus);
	if (!drs) {
		fprintf (stderr, "Error: unable to open record_in set %s\n", open_dsname);
 	   return(0);
	}
	ncnt = drs->n;
	if (ncnt < 1) {
		fprintf (stderr, "No records in selected set %s\n", open_dsname);
		return(0);
	}
    drs_out = drms_create_records(drms_env, ncnt, dsout, DRMS_PERMANENT,&rstatus);

	
	int lf_retcode;
	
/***********************************************************************

************************************************************************/
    char *hwltnset;
    char s_hwltnset[10];
	unsigned int fsnx = 0;
	long long recnum0;
	int seg_ct;
	double x0_lf, y0_lf, rsun_lf;
	
	int cam,fid,hpltid,missvals;
	TIME tobs;
	float sat_rot;//,inst_rot;
//	char s_tobs[64];


    for(r=0; r < ncnt; r++) 
    {
		record_in = drs->records[r];
		record_out = drs_out->records[r]; 

		recnum0 = record_in->recnum;
		fsnx = drms_getkey_int(record_in, "FSN", &rstatus);       
		printf ("rec_num:%lld FSN:%u\n", recnum0, fsnx);
		fid = drms_getkey_int(record_in, "FID", &rstatus);
		if(rstatus) {
			fprintf (stderr, "ERROR: in drms_getkey_int(FID) fid=%d\n", fid);
		}
		//printf("FID: %d\n", fid);
		hwltnset = drms_getkey_string(record_in, "HWLTNSET", &rstatus);
		if(rstatus) {
			fprintf (stderr, "ERROR: in drms_getkey_int(HWLTNSET) hwltnset=%s\n", hwltnset);
		}
		missvals = drms_getkey_int(record_in, "MISSVALS", &rstatus);
		if(rstatus) {
			fprintf (stderr, "ERROR: in drms_getkey_int(MISSVALS) missvals=%d\n", missvals);
		}
		sprintf(s_hwltnset,hwltnset);
		//printf("HWLTNSET: -%s-\n", s_hwltnset);
		int a=(strncmp(hwltnset,"CLOS",4));
   		//printf("cmp hwltnset fid %d %s %d\n",a,hwltnset,fid);
   		if (a==0 && (missvals==0) )
   		{
			//only those needed by limbfit
			tobs = drms_getkey_time(record_in, "t_obs", &rstatus);
			if(rstatus) {
				fprintf(opf,"Error on drms_getkey_time() fsn=%u. Use DRMS_MISSING_TIME\n", fsnx);
			}
			cam = drms_getkey_int(record_in, "CAMERA", &rstatus);
			if(rstatus) {
				fprintf(opf,"ERROR: in drms_getkey_int(CAMERA) cam=%d\n", cam);
			}
			hpltid = drms_getkey_int(record_in, "hpltid", &rstatus);
			  if(rstatus) {
				fprintf(opf,"ERROR: in drms_getkey_int(HPLTID) hpltid=%d\n", hpltid);
			}
			x0_lf = drms_getkey_float(record_in, "X0_LF", &rstatus);
			if(rstatus) {
				fprintf(opf,"ERROR: in drms_getkey_int(X0_LF) x0_lf=%f\n", x0_lf);
			}
			y0_lf = drms_getkey_float(record_in, "Y0_LF", &rstatus);
			if(rstatus) {
				fprintf(opf,"ERROR: in drms_getkey_int(Y0_LF) y0_lf=%f\n", y0_lf);
			}
			rsun_lf = drms_getkey_float(record_in, "RSUN_LF", &rstatus);
			  if(rstatus) {
				fprintf(opf,"ERROR: in drms_getkey_int(RSUN_LF) rsun_lf=%f\n", rsun_lf);
			}
			sat_rot = drms_getkey_float(record_in, "SAT_ROT", &rstatus);
			  if(rstatus) {
				fprintf(opf,"ERROR: in drms_getkey_int(SAT_ROT) sat_rot=%f\n", sat_rot);
			}
			seg_ct = drms_record_numsegments (record_in);
			if (seg_ct < 1) 
			{
				fprintf (opf,"  no segments in selected record_in\n");
				//what to do then?
			}
			segment = drms_segment_lookupnum(record_in, 0);
			img = drms_segment_read(segment, DRMS_TYPE_FLOAT, &rstatus);
			if(!img) {
				fprintf(opf,"Can't do drms_segment_read() status=%d\n", rstatus);
				return(1);   //what to do exactly?
			} else
			{
				lfv->data=img->data;
				lfv->img_sz0=img->axis[0];
				lfv->img_sz1=img->axis[1];
				lfv->ix=x0_lf;
				lfv->iy=y0_lf;
				lfv->ir=rsun_lf;
				lfv->camera=cam;
				lfv->fid=fid;
				lfv->fsn=fsnx;
				lfv->hpltid=hpltid;
				lfv->dsout=dsout;
				lfv->opf=opf;
//				lfv->tobs=tobs;
				lfv->sat_rot=sat_rot;
		
				lf_retcode=do_one_limbfit(lfv,lfr,record_in,record_out,debug);
				drms_free_array (img);
			}
		} // if loop ok & no missvals
		else fprintf(opf,"REJECTED:(1) rec_num:%lld FSN:%u\n", recnum0, fsnx);
	} //end process one DRMS record_in : for(r=0; r < ncnt; r++) 

	drms_close_records(drs_out, DRMS_INSERT_RECORD);
	drms_close_records(drs, DRMS_FREE_RECORD);
	fprintf(opf,"records save\n");
	printf("records save\n");
/**********************************************************************/
	//free(data); 
	//free(filenameout);
	
	//+ free all variables coming from limbfit.c
/*	printf("ok1\n");	
	free(lfr->fits_ldfs_data);
	printf("ok1\n");	
	free(lfr->fits_as);
	printf("ok1\n");	
	free(lfr->fits_es);
	printf("ok1\n");	
	free(lfr->fits_r);
	printf("ok1\n");	
	free(lfr->fits_ab_dataa);
	printf("ok1\n");	
	free(lfr->fits_ab_datab);
	printf("ok1\n");	
	free(lfr->fits_ip);
	printf("ok1\n");	
*/	
	fprintf(opf,"end of main\n");
	
return(0); // to change that
}

// ajouter un retcode/quality if sthg wrong happened but not serious


int DoIt(void)
{
	//************************************************************************
	// usage ./lfwrp bfsn efsn [dsin] [dsout] [debug]
	//************************************************************************
	CmdParams_t *params = &cmdparams;
	
	int  debug = params_get_int (params, "debug");
	char* dsin = params_get_str (params, "dsin");
	char* dsout = params_get_str (params, "dsout");
	long long bfsn = params_get_int (params, "bfsn");
	long long efsn = params_get_int (params, "efsn");

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

	/***********************************************************************
		LOG FILE: MAKE A GLOBAL LOG FILE FOR THE WHOLE SESSION: 
			change the close statement everywhere....
				change free too
	************************************************************************/

	static char flogname[128];
	FILE *opf;
//	sprintf(flogname, "%slog_%d_%d.log", LOG_DIR,bfsn,efsn);
	sprintf(flogname, "./log_%d_%d.log",bfsn,efsn);
    if((opf=fopen(flogname, "w")) == NULL)
    {
		fprintf(stderr, "**Can't open the log file %s\n", flogname);
		return(0);
	}
	// make it quick for now but add a time tag to the file name %slog_%s_%s
	// if can't open then stop

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
		if(process_n_records(open_dsname, dsout, opf, debug)) 
		{  //do a chunk to get files from the lev0
			fprintf(opf,"lfwrp abort\nSee log: %s\n", flogname); 
			//send_mail("build_lev1 abort\nSee log: %s\n", logname); 
			return(0);
		}
		drms_server_end_transaction(drms_env,0,0);
		drms_server_begin_transaction(drms_env);
	}
	fclose(opf);

return(0);
} //end doit()
    
