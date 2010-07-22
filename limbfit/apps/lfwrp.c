/* 
	I.Scholl 

		High level wrapper for the limbfit code
		Can be used as a stand-alone module (includes a 'DoIt' function)
			calls: do_one_limbfit() that calls:
										limbfit() and reformat output data to be rewritten in the DB
			do_one_limbfit() can be plugged into lev0lev1 module
			
	
	#define CODE_NAME 		"limbfit"
	#define CODE_VERSION 	"V0.2r02" 
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
//	long long bnumx, enumx;

	sprintf(recrange, "%lld-%lld", bfsn, efsn);
	sprintf(open_dsname, "%s[%s]", dsin, recrange);
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
	Open DRMS connexion, get the range of data, load structure and run lf
************************************************************************/
	DRMS_RecordSet_t *drs;
	DRMS_Record_t *record;
	DRMS_Segment_t *segment;
	int rstatus, ncnt,r;

    drs = drms_open_records(drms_env, open_dsname, &rstatus);
	if (!drs) {
		fprintf (stderr, "Error: unable to open record set %s\n", open_dsname);
 	   return(0);
	}
	ncnt = drs->n;
	if (ncnt < 1) {
		fprintf (stderr, "No records in selected set %s\n", open_dsname);
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
	//make it quick for now but add a time tag to the file name %slog_%s_%s

	// if can't open then stop
	
/***********************************************************************
	Initialization of the structure to pass to the limbfit routine
************************************************************************/
	static LIMBFIT_INPUT limbfit_vars ;
	static LIMBFIT_INPUT *lfv = &limbfit_vars;
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
		record = drs->records[r];
		recnum0 = record->recnum;
		fsnx = drms_getkey_int(record, "FSN", &rstatus);       
		printf ("rec_num:%lld FSN:%u\n", recnum0, fsnx);
		fid = drms_getkey_int(record, "FID", &rstatus);
		if(rstatus) {
			fprintf (stderr, "ERROR: in drms_getkey_int(FID) fid=%d\n", fid);
		}
		printf("FID: %d\n", fid);
		hwltnset = drms_getkey_string(record, "HWLTNSET", &rstatus);
		if(rstatus) {
			fprintf (stderr, "ERROR: in drms_getkey_int(HWLTNSET) hwltnset=%s\n", hwltnset);
		}
		missvals = drms_getkey_int(record, "MISSVALS", &rstatus);
		if(rstatus) {
			fprintf (stderr, "ERROR: in drms_getkey_int(MISSVALS) missvals=%d\n", missvals);
		}
		sprintf(s_hwltnset,hwltnset);
		printf("HWLTNSET: -%s-\n", s_hwltnset);
		int a=(strncmp(hwltnset,"CLOS",4));
   		printf("cmp hwltnset fid %d %s %d\n",a,hwltnset,fid);
   		if (a==0 && (missvals==0) )
   		{
			//only those needed by limbfit
			tobs = drms_getkey_time(record, "t_obs", &rstatus);
			if(rstatus) {
				fprintf(opf,"Error on drms_getkey_time() fsn=%u. Use DRMS_MISSING_TIME\n", fsnx);
			}
			cam = drms_getkey_int(record, "CAMERA", &rstatus);
			if(rstatus) {
				fprintf(opf,"ERROR: in drms_getkey_int(CAMERA) cam=%d\n", cam);
			}
			hpltid = drms_getkey_int(record, "hpltid", &rstatus);
			  if(rstatus) {
				fprintf(opf,"ERROR: in drms_getkey_int(HPLTID) hpltid=%d\n", hpltid);
			}
			x0_lf = drms_getkey_float(record, "X0_LF", &rstatus);
			if(rstatus) {
				fprintf(opf,"ERROR: in drms_getkey_int(X0_LF) x0_lf=%f\n", x0_lf);
			}
			y0_lf = drms_getkey_float(record, "Y0_LF", &rstatus);
			if(rstatus) {
				fprintf(opf,"ERROR: in drms_getkey_int(Y0_LF) y0_lf=%f\n", y0_lf);
			}
			rsun_lf = drms_getkey_float(record, "RSUN_LF", &rstatus);
			  if(rstatus) {
				fprintf(opf,"ERROR: in drms_getkey_int(RSUN_LF) rsun_lf=%f\n", rsun_lf);
			}
			sat_rot = drms_getkey_float(record, "SAT_ROT", &rstatus);
			  if(rstatus) {
				fprintf(opf,"ERROR: in drms_getkey_int(SAT_ROT) sat_rot=%f\n", sat_rot);
			}
			seg_ct = drms_record_numsegments (record);
			if (seg_ct < 1) 
			{
				fprintf (opf,"  no segments in selected record\n");
				//what to do then?
			}
			segment = drms_segment_lookupnum(record, 0);
			DRMS_Array_t *img;
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
		
				lf_retcode=do_one_limbfit(lfv,record,debug);
			}
		} // if loop ok & no missvals
		else fprintf(opf,"REJECTED:(1) rec_num:%lld FSN:%u\n", recnum0, fsnx);
	} //end process one DRMS record : for(r=0; r < ncnt; r++) 
						
		
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
	fclose(opf);
	
return(0); //change that
}

