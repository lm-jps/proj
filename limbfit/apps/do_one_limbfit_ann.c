/* 
	I.Scholl 

		Second level wrapper for the limbfit code
		Can be used as a callable module (to be plugged into lev0lev1 module or by the main wrapper)
			calls: limbfit() and reformat output data to be rewritten in the DB
			
	
	#define results->code_name 		"limbfit_tas"
	#define CODE_VERSION 			"V5.04" 
	#define CODE_DATE 				"Tue May  7 21:27:32 PDT 2013" 
*/

#include "limbfit_tas.h"

int do_one_limbfit(DRMS_Record_t *record_in,DRMS_Record_t *record_out, 
					LIMBFIT_INPUT *input, LIMBFIT_OUTPUT *results, LIMBFIT_IO_PUT *ios, int *status)
{
	static char *log_msg_code="do_one_limbfit";
	static char *series_name="limbfit data";
	results->series_name=series_name;
	//static char *proc_stat=PROCSTAT_OK;

 	int lf_retcode;
	int rstatus;
	char log_msg[200];
	
	int seg_cnt;
	DRMS_Segment_t *segment_in;
	DRMS_Segment_t *segment_out;
	static DRMS_Array_t *img;

	//if (debug) lf_logmsg("DEBUG", "INFO", 0, 0, "IN do___", log_msg_code, results->opf);

	seg_cnt = drms_record_numsegments (record_in);
	if (seg_cnt < 1) 
	{
		sprintf(log_msg,"no segments in selected record_in");
		lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_DATA, 0, log_msg, log_msg_code,results->opf);
		write_mini_output(PROCSTAT_NO_LF_DB_READ_PB,record_in,record_out,0,results);
		return(ERR_DRMS_READ_MISSING_DATA);   
	}
	segment_in = drms_segment_lookupnum(record_in, 0);
	img = drms_segment_read(segment_in, DRMS_TYPE_FLOAT, &rstatus);
	if(!img) {
		sprintf(log_msg,"Can't do drms_segment_read() status=%d", rstatus);
		lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_DATA, 0, log_msg, log_msg_code, results->opf);
		write_mini_output(PROCSTAT_NO_LF_DB_READ_PB,record_in,record_out,0,results);
		return(ERR_DRMS_READ_MISSING_DATA);   
	} 
	else
	{
		//if (debug) lf_logmsg("DEBUG", "INFO", 0, 0, "IN do___getting data", log_msg_code, results->opf);
		input->ix = drms_getkey_float(record_in, "X0_LF", &rstatus);
		if(rstatus) {
			lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_DATA,rstatus, "drms_getkey_float(X0_LF)", log_msg_code, results->opf);
			write_mini_output(PROCSTAT_NO_LF_DB_READ_PB,record_in,record_out,0,results);
			return(ERR_DRMS_READ_MISSING_DATA);   
		}
		if(isnan(input->ix)) {
			lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_XYR_LF,rstatus, "X0_LF missing", log_msg_code, results->opf);
			write_mini_output(PROCSTAT_NO_LF_XYR_LF_MISSING,record_in,record_out,0,results);
			return(ERR_DRMS_READ_MISSING_XYR_LF);   
		}
		input->iy = drms_getkey_float(record_in, "Y0_LF", &rstatus);
		if(rstatus) {
			lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_DATA,rstatus, "drms_getkey_float(Y0_LF)", log_msg_code, results->opf);
			write_mini_output(PROCSTAT_NO_LF_DB_READ_PB,record_in,record_out,0,results);
			return(ERR_DRMS_READ_MISSING_DATA);   
		}
		if(isnan(input->iy)) {
			lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_XYR_LF,rstatus, "Y0_LF missing", log_msg_code, results->opf);
			write_mini_output(PROCSTAT_NO_LF_XYR_LF_MISSING,record_in,record_out,0,results);
			return(ERR_DRMS_READ_MISSING_XYR_LF);   
		}
		input->ir = drms_getkey_float(record_in, "RSUN_LF", &rstatus);
		  if(rstatus) {
			lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_DATA,rstatus, "drms_getkey_float(RSUN_LF)", log_msg_code, results->opf);
			write_mini_output(PROCSTAT_NO_LF_DB_READ_PB,record_in,record_out,0,results);
			return(ERR_DRMS_READ_MISSING_DATA);   
		}
		if(isnan(input->ir)) {
			lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_XYR_LF,rstatus, "RSUN_LF missing", log_msg_code, results->opf);
			write_mini_output(PROCSTAT_NO_LF_XYR_LF_MISSING,record_in,record_out,0,results);
			return(ERR_DRMS_READ_MISSING_XYR_LF);   
		}
		input->img_sz0=img->axis[0];
		input->img_sz1=img->axis[1];
		input->data=img->data;

		if (results->debug) 
		{
			sprintf(log_msg,"X0_LF %f", input->ix);
			lf_logmsg("DEBUG", "INFO", 0, 0, log_msg, log_msg_code, results->opf);
			sprintf(log_msg,"Y0_LF %f", input->iy);
			lf_logmsg("DEBUG", "INFO", 0, 0, log_msg, log_msg_code, results->opf);
			sprintf(log_msg,"RSUN_LF %f", input->ir);
			lf_logmsg("DEBUG", "INFO", 0, 0, log_msg, log_msg_code, results->opf);
		}
/* add annulus code here! */

		lf_retcode=limbfit(input,results,ios);
		
		drms_free_array(img);

		if (results->debug) 
		{
			sprintf(log_msg,"XC %f", results->cenx);
			lf_logmsg("DEBUG", "INFO", 0, 0, log_msg, log_msg_code, results->opf);
			sprintf(log_msg,"YC %f", results->ceny);
			lf_logmsg("DEBUG", "INFO", 0, 0, log_msg, log_msg_code, results->opf);
			sprintf(log_msg,"R %f", results->radius);
			lf_logmsg("DEBUG", "INFO", 0, 0, log_msg, log_msg_code, results->opf);
			sprintf(log_msg,"Returned code %d", lf_retcode);
			lf_logmsg("DEBUG", "INFO", 0, 0, log_msg, log_msg_code, results->opf);
			sprintf(log_msg,"err1 %d", results->error1);
			lf_logmsg("DEBUG", "INFO", 0, 0, log_msg, log_msg_code, results->opf);
			sprintf(log_msg,"err2 %d", results->error2);
			lf_logmsg("DEBUG", "INFO", 0, 0, log_msg, log_msg_code, results->opf);
		}
		
		// >= 0 => all
		// ERR_LIMBFIT_FIT_FAILED => all
		// ERR_LIMBFIT_FAILED => mini

		if (input->fldf==0) results->fldfr=4;

		if (lf_retcode >= 0 || lf_retcode == ERR_LIMBFIT_FIT_FAILED)
		{	
			if (results->debug) printf("write in the db\n");
			//**********************************************************************
			//	Write into the DB
			//**********************************************************************
			
			// get/set all keywords

			drms_copykeys(record_out, record_in, 1, kDRMS_KeyClass_Explicit); 		
			// change DATE to our proc_date, do not keep level 1 DATE:
			drms_keyword_setdate(record_out); 
			write_lf_keywords(VOID,record_out, results,0);

			// write the segments 
			*status=0;
			rstatus=0;

			int naxes[2];
			
			naxes[0]=results->fits_ldfs_naxis1;
			naxes[1]=results->fits_ldfs_naxis2;
			DRMS_Array_t *array1=drms_array_create(DRMS_TYPE_FLOAT,2,naxes,results->fits_ldfs_data,&rstatus);			
			if (rstatus) 
			{
				lf_logmsg("ERROR", "DRMS", ERR_DRMS_ARRAY_CREATE,rstatus, "ERR_DRMS_ARRAY_CREATE(ldfs)", log_msg_code, results->opf);
				return ERR_DRMS_ARRAY_CREATE;
			}
			segment_out = drms_segment_lookupnum(record_out, 0); 
			if (segment_out)
			{
				if ( drms_segment_write(segment_out,array1,0) ) 
				{
					lf_logmsg("ERROR", "DRMS", ERR_DRMS_SEGMENT_WRITE,rstatus, "ERR_DRMS_SEGMENT_WRITE(ldfs)", log_msg_code, results->opf);
					return ERR_DRMS_SEGMENT_WRITE;
				}	
			}
			else
			{
				lf_logmsg("ERROR", "DRMS", ERR_DRMS_SEGMENT_LOOKUPNUM,rstatus, "ERR_DRMS_SEGMENT_LOOKUPNUM(ldfs)", log_msg_code, results->opf);
				return ERR_DRMS_SEGMENT_LOOKUPNUM;
			}
			drms_free_array (array1);
			// ????to check-> NOTE: because of that, from now tbf MUST BE EQUAL to 1
		//add more KWs
		//ORIGIN kw ok?			
		//TIME too?
			/*
			TIME tobs;
			static char s_tobs[50];
			tobs = drms_getkey_time(record_in, "DATE__OBS", &rstatus);
			sprint_time(s_tobs, tobs, "UTC", 0);
			drms_setkey_string(record_out, "DATE__OBS",	&s_tobs);

			tobs = drms_getkey_time(record_in, "T_OBS", &rstatus);
			sprint_time(s_tobs, tobs, "UTC", 0);
			drms_setkey_string(record_out, "T_OBS",	s_tobs);
			*/
			
			rstatus=0;
		
			// AB ---EXT#1---------------------------------------------------------------------
			//if (results->nb_iter==1)
			//{
			// attention only one beta (the last one) => to be changed in limbfit.c too
				naxes[0]=results->fits_ab_naxis1;
				naxes[1]=results->fits_ab_naxis2;
				DRMS_Array_t *array2=drms_array_create(DRMS_TYPE_FLOAT,2,naxes,results->fits_alpha_beta,&rstatus);			
				if (rstatus) 
				{
					lf_logmsg("ERROR", "DRMS", ERR_DRMS_ARRAY_CREATE,rstatus, "ERR_DRMS_ARRAY_CREATE(ab)", log_msg_code, results->opf);
					return ERR_DRMS_ARRAY_CREATE;
				}
				segment_out = drms_segment_lookupnum(record_out, 1); 
				if (segment_out)
				{
					if ( drms_segment_write(segment_out,array2,0) ) 
					{
					lf_logmsg("ERROR", "DRMS", ERR_DRMS_SEGMENT_WRITE,rstatus, "ERR_DRMS_SEGMENT_WRITE(ab)", log_msg_code, results->opf);
					return ERR_DRMS_SEGMENT_WRITE;
					}	
				}
				else
				{
					lf_logmsg("ERROR", "DRMS", ERR_DRMS_SEGMENT_LOOKUPNUM,rstatus, "ERR_DRMS_SEGMENT_LOOKUPNUM(ab)", log_msg_code, results->opf);
					return ERR_DRMS_SEGMENT_LOOKUPNUM;
				}
				drms_free_array (array2);

			
			if (input->fldf==1)
			{
				// Full LDF ---EXT#3-------------------------------------------------------------------------
				naxes[0]=results->fits_fldfs_nrows;
				naxes[1]=results->fits_fldfs_tfields;
				DRMS_Array_t *array4=drms_array_create(DRMS_TYPE_FLOAT,2,naxes,results->fits_fulldfs,&rstatus);			
				if (rstatus) 
				{
					lf_logmsg("ERROR", "DRMS", ERR_DRMS_ARRAY_CREATE,rstatus, "ERR_DRMS_ARRAY_CREATE(fullldfs)", log_msg_code, results->opf);
					return ERR_DRMS_ARRAY_CREATE;
				}
				segment_out = drms_segment_lookupnum(record_out, 3); 
				if (segment_out)
				{
					if ( drms_segment_write(segment_out,array4,0) ) 
					{
						lf_logmsg("ERROR", "DRMS", ERR_DRMS_SEGMENT_WRITE,rstatus, "ERR_DRMS_SEGMENT_WRITE(fullldfs)", log_msg_code, results->opf);
						return ERR_DRMS_SEGMENT_WRITE;
					}	
				}
				else
				{
					lf_logmsg("ERROR", "DRMS", ERR_DRMS_SEGMENT_LOOKUPNUM,rstatus, "ERR_DRMS_SEGMENT_LOOKUPNUM(fullldfs)", log_msg_code, results->opf);
					return ERR_DRMS_SEGMENT_LOOKUPNUM;
				}
				drms_free_array (array4);
			}

			//---------------------------------------
			//		5) close & free
			//---------------------------------------
			rstatus=drms_setkey_string(record_out, "PROCSTAT", PROCSTAT_OK);

			if (results->debug) printf("end write in the db %d\n",rstatus);
			
			/* no need this because of the drms_free_array(array)	
			free(results->fits_ldfs_data);
			free(results->fits_alpha_beta);
			free(results->fits_params);
			if (results->fldfr<3) free(results->fits_fulldfs);
			*/
			if (results->debug) printf("done...%u\n",input->fsn);							
			lf_logmsg("INFO", "APP", 0, 0, "End writing in the DB", log_msg_code, results->opf);
		}
		else 
		{
			sprintf(log_msg,"NO LDFS file created (ifail >0) or retcode <0 for fsn# %u", input->fsn);
			lf_logmsg("WARNING", "APP", ERR_LIMBFIT_FAILED, lf_retcode, log_msg, log_msg_code, results->opf);			
			write_mini_output(PROCSTAT_NO_LF_FAILED,record_in,record_out,0,results);
			// is this the right retcode? or should I make a test???
			
			// free 
			if (results->error1 == 0)
			{
				free(results->fits_ldfs_data);
				free(results->fits_alpha_beta);
				if (results->fldfr<3) free(results->fits_fulldfs);
			}
			return(0);
		}
	}

return 0;
}

int	write_mini_output(char * errcode, DRMS_Record_t *record_in,DRMS_Record_t *record_out, 
			int tbf, LIMBFIT_OUTPUT *results)
{
		//int status;
		char *log_msg_code="write_mini_output";
		lf_logmsg("INFO", "APP", 0, 0, "Writing only header", log_msg_code, results->opf);			
		// write only the header of segment in the DB
		drms_copykeys(record_out, record_in, 1, kDRMS_KeyClass_Explicit); 	
		// change DATE to our proc_date, do not keep level 1 DATE:
		drms_keyword_setdate(record_out); 

		write_lf_keywords(errcode,record_out, results,1);
		if (tbf >=1)
		{
			free(results->fits_alpha_beta);
			if (results->fldfr<3) free(results->fits_fulldfs);
			if (tbf ==2) free(results->fits_ldfs_data);	
		}
return(0);
}

int	write_lf_keywords(char * errcode, DRMS_Record_t *record_out, LIMBFIT_OUTPUT *results, 
				int pass)
{
// need to test status...
		drms_setkey_string(record_out, "SERIESCN",  results->series_name);
		drms_setkey_string(record_out, "INSERIES", 	results->dsin);
		drms_setkey_string(record_out, "CODENAME",	results->code_name);
		drms_setkey_string(record_out, "CODEVERS",	results->code_version);
		drms_setkey_string(record_out, "CODEDATE",	results->code_date);
		drms_setkey_float (record_out, "X_LFS",		results->cenx);
		drms_setkey_float (record_out, "Y_LFS",		results->ceny);
		drms_setkey_double(record_out, "R_LFS",		results->radius);
		drms_setkey_int   (record_out, "QUAL_LFS",	results->quality);
		drms_setkey_int   (record_out, "ERRO_LFS",	results->error1);
		drms_setkey_int   (record_out, "ERRO_FIT",	results->error2);
		drms_setkey_int   (record_out, "NB_FBINS",	results->nb_fbins);
		drms_setkey_int   (record_out, "FLDFR",		results->fldfr);
		drms_setkey_float (record_out, "ES0",		results->fits_es[0]);
		drms_setkey_float (record_out, "ES1",		results->fits_es[1]);
		drms_setkey_float (record_out, "ES2",		results->fits_es[2]);
		drms_setkey_float (record_out, "ES3",		results->fits_es[3]);
		drms_setkey_float (record_out, "ES4",		results->fits_es[4]);
		drms_setkey_float (record_out, "ES5",		results->fits_es[5]);
		drms_setkey_float (record_out, "AS0",		results->fits_as[0]);
		drms_setkey_float (record_out, "AS1",		results->fits_as[1]);
		drms_setkey_float (record_out, "AS2",		results->fits_as[2]);
		drms_setkey_float (record_out, "AS3",		results->fits_as[3]);
		drms_setkey_float (record_out, "AS4",		results->fits_as[4]);
		drms_setkey_float (record_out, "AS5",		results->fits_as[5]);
		drms_setkey_double(record_out, "CMEAN",		results->cmean);
		drms_setkey_int   (record_out, "ANN_WD",	results->ann_wd);
		drms_setkey_int   (record_out, "MXSZANNV",	results->mxszannv);
		drms_setkey_int   (record_out, "NB_LDF",	results->nb_ldf);
		drms_setkey_int   (record_out, "NB_RDB",	results->nb_rdb);
		drms_setkey_int   (record_out, "NB_ABB",	results->nb_abb);
		drms_setkey_double(record_out, "UP_LIMIT",	results->up_limit);
		drms_setkey_double(record_out, "LO_LIMIT",	results->lo_limit);
		drms_setkey_double(record_out, "INC_X",		results->inc_x);
		drms_setkey_double(record_out, "INC_Y",		results->inc_y);
		drms_setkey_int   (record_out, "NFITPNTS",	results->nfitpnts);
		drms_setkey_int   (record_out, "NB_ITER",	results->nb_iter);
		drms_setkey_int   (record_out, "CEN_CALC",	results->cc);
		drms_setkey_double(record_out, "AHI",		results->ahi);
		drms_setkey_string(record_out, "COMMENTS",	results->comment);			
		drms_setkey_string(record_out, "BLD_VERS",  results->bld_vers);

//where this one comes from? if(get_set_kw(1,"ORIGIN","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);

		if (pass==1)
			drms_setkey_string(record_out, "PROCSTAT",	errcode);

return(0);
}
