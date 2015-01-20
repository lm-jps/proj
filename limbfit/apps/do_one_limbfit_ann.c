/* 
	I.Scholl 

		Second level wrapper for the limbfit code
		Can be used as a callable module (to be plugged into lev0lev1 module or by the main wrapper)
			calls: limbfit() and reformat output data to be rewritten in the DB
			
	
	#define results->code_name 		"limbfit_ann"
	#define CODE_VERSION 			"V1.00" 
	#define CODE_DATE 				"Mon Sep 15 15:14:19 PDT 2014" 
*/

#include "limbfit_ann.h"

int do_one_limbfit(DRMS_Record_t *record_in,DRMS_Record_t *record_out, 
					LIMBFIT_INPUT *input, LIMBFIT_OUTPUT *results, int *status)
{
	static char *log_msg_code="do_one_limbfit_ann";
	static char *series_name="limbfit annulus data";
	results->series_name=series_name;
	//static char *proc_stat=PROCSTAT_OK;

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
		write_mini_output(PROCSTAT_NO_LA_DB_READ_PB,record_in,record_out,0,results);
		return(ERR_DRMS_READ_MISSING_DATA);   
	}
	segment_in = drms_segment_lookupnum(record_in, 0);
	img = drms_segment_read(segment_in, DRMS_TYPE_INT, &rstatus);
	if(!img) {
		sprintf(log_msg,"Can't do drms_segment_read() status=%d", rstatus);
		lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_DATA, 0, log_msg, log_msg_code, results->opf);
		write_mini_output(PROCSTAT_NO_LA_DB_READ_PB,record_in,record_out,0,results);
		return(ERR_DRMS_READ_MISSING_DATA);   
	} 
	else
	{
		input->ix = drms_getkey_float(record_in, "X0_LF", &rstatus);
		if(rstatus) {
			lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_DATA,rstatus, "drms_getkey_float(X0_LF)", log_msg_code, results->opf);
			write_mini_output(PROCSTAT_NO_LA_DB_READ_PB,record_in,record_out,0,results);
			return(ERR_DRMS_READ_MISSING_DATA);   
		}
		if(isnan(input->ix)) {
			lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_XYR_LA,rstatus, "X0_LF missing", log_msg_code, results->opf);
			write_mini_output(PROCSTAT_NO_LA_XYR_LF_MISSING,record_in,record_out,0,results);
			return(ERR_DRMS_READ_MISSING_XYR_LA);   
		}
		input->iy = drms_getkey_float(record_in, "Y0_LF", &rstatus);
		if(rstatus) {
			lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_DATA,rstatus, "drms_getkey_float(Y0_LF)", log_msg_code, results->opf);
			write_mini_output(PROCSTAT_NO_LA_DB_READ_PB,record_in,record_out,0,results);
			return(ERR_DRMS_READ_MISSING_DATA);   
		}
		if(isnan(input->iy)) {
			lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_XYR_LA,rstatus, "Y0_LF missing", log_msg_code, results->opf);
			write_mini_output(PROCSTAT_NO_LA_XYR_LF_MISSING,record_in,record_out,0,results);
			return(ERR_DRMS_READ_MISSING_XYR_LA);   
		}
		/* add annulus code here! */

		/************************************************************************/
		/*                        set parameters                               */
		/************************************************************************/
		int		naxis_row=img->axis[0];
		int		naxis_col=img->axis[1];
		long i,j;
		int *data=img->data;

		/************************************************************************/
		/*          Compute the mean of the center of the image                 */
		/************************************************************************/
		// doesn't work anymore: how to convert the INT to FLOAT????

		float cmean,ctot=0.0;
		long nbp=0;
		float limx_m=(float)input->ix-250;
		float limx_p=(float)input->ix+250;
		float limy_m=(float)input->iy-250;
		float limy_p=(float)input->iy+250;
		for (i=limx_m;i<limx_p;i++)
		{
			for (j=limy_m;j<limy_p;j++)
			{
				ctot=ctot+data[i*naxis_col+j]; 
				nbp++;
			}
		}
		cmean=ctot/nbp;
		results->cmean=cmean;
		if (results->debug)
		{
			sprintf(log_msg," cmean = %6.4f (ctot= %6.4f , nbp=%ld)", cmean,ctot,nbp);
			lf_logmsg("DEBUG", "APP", 0, 0, log_msg, log_msg_code, results->opf);
		}
		/************************************************************************/
		/*                        Make annulus data                             */
		/************************************************************************/

		/* Select Points */

			int *p_data=&data[0];		
			int *p_mask=input->pf_mask;		
			while(p_mask<=input->pl_mask) 
			{
				if (*(p_mask) == 0)
				{
					*(p_data)=EQNANVAL;
				}	
				*(p_data++);
				*(p_mask++);
			}

/* end annulus */		

		
			if (results->debug) printf("write in the db\n");
			//**********************************************************************
			//	Write into the DB
			//**********************************************************************
			
			// get/set all keywords

			drms_copykeys(record_out, record_in, 1, kDRMS_KeyClass_Explicit); 		
			// change DATE to our proc_date, do not keep level 1 DATE:
			drms_keyword_setdate(record_out); 
			write_lf_keywords(VOID,record_out, results,0);

			// write the segment
			*status=0;
			rstatus=0;

			int naxes[2];
			naxes[0]=naxis_row;
			naxes[1]=naxis_col;
			
			DRMS_Array_t *array1=drms_array_create(DRMS_TYPE_INT,2,naxes,data,&rstatus);			
			if (rstatus) 
			{
				lf_logmsg("ERROR", "DRMS", ERR_DRMS_ARRAY_CREATE,rstatus, "ERR_DRMS_ARRAY_CREATE(annulus)", log_msg_code, results->opf);
				return ERR_DRMS_ARRAY_CREATE;
			}
			segment_out = drms_segment_lookupnum(record_out, 0); 
			if (segment_out)
			{
				if ( drms_segment_write(segment_out,array1,0) ) 
				{
					lf_logmsg("ERROR", "DRMS", ERR_DRMS_SEGMENT_WRITE,rstatus, "ERR_DRMS_SEGMENT_WRITE(annulus)", log_msg_code, results->opf);
					return ERR_DRMS_SEGMENT_WRITE;
				}	
			}
			else
			{
				lf_logmsg("ERROR", "DRMS", ERR_DRMS_SEGMENT_LOOKUPNUM,rstatus, "ERR_DRMS_SEGMENT_LOOKUPNUM(annulus)", log_msg_code, results->opf);
				return ERR_DRMS_SEGMENT_LOOKUPNUM;
			}
			drms_free_array (array1);
//			drms_free_array (img);
			
			rstatus=0;
		

			//---------------------------------------
			//		5) close & free
			//---------------------------------------
			rstatus=drms_setkey_string(record_out, "PROCSTAT", PROCSTAT_OK);

			if (results->debug) printf("end write in the db %d\n",rstatus);
			
			if (results->debug) printf("done...%u\n",input->fsn);							
			lf_logmsg("INFO", "APP", 0, 0, "End writing in the DB", log_msg_code, results->opf);
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
		drms_setkey_float (record_out, "CENT_X",	results->cenx);
		drms_setkey_float (record_out, "CENT_Y",	results->ceny);
		drms_setkey_float (record_out, "R_MIN",		results->r_min);
		drms_setkey_float (record_out, "R_MAX",		results->r_max);
		drms_setkey_double(record_out, "CMEAN",		results->cmean);
		drms_setkey_string(record_out, "COMMENTS",	results->comment);			
		drms_setkey_string(record_out, "BLD_VERS",  results->bld_vers);

//where this one comes from? if(get_set_kw(1,"ORIGIN","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);

		if (pass==1)
			drms_setkey_string(record_out, "PROCSTAT",	errcode);

return(0);
}
