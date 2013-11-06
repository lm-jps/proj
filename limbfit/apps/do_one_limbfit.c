/* 
	I.Scholl 

		Second level wrapper for the limbfit code
		Can be used as a callable module (to be plugged into lev0lev1 module or by the main wrapper)
			calls: limbfit() and reformat output data to be rewritten in the DB
			
	
	#define results->code_name 		"limbfit"
	#define results->code_version 	"V4.0r8" 
	#define results->code_date 		"Wed Nov  6 09:42:51 HST 2013" 
*/

#include "limbfit.h"

void lf_logmsg4fitsio(char *log_msg,char *log_msg_code,char *kw,unsigned int fsn,int status, FILE *opf)
{	
	sprintf(log_msg,"can't fits_update_key[%s] for fsn# %u", kw,fsn);
	lf_logmsg("ERROR", "FITSIO", ERR_FITSIO, status, log_msg, log_msg_code, opf);			
	fits_report_error(opf, status); 
}

int get_set_kw(int typ, char *kw, char *kw_txt, unsigned int fsn, 
		DRMS_Record_t *record_in,DRMS_Record_t *record_out, fitsfile *outfptr, int tbf, LIMBFIT_OUTPUT *results, int *status)
{
	float kw_f;
	double kw_d;
	int kw_i;
	char *kw_c;
	int rstatus;
	char log_msg[200];
	char *log_msg_code="get_set_kw";

	switch( typ ) 
	{
	    case 1: 
			kw_c = drms_getkey_string(record_in, kw, &rstatus);
			break;
	    case 2: 
			kw_i = drms_getkey_int(record_in, kw, &rstatus);
			break;
	    case 3: 
			kw_f = drms_getkey_float(record_in, kw, &rstatus);
			break;
	    case 4: 
			kw_d = drms_getkey_double(record_in, kw, &rstatus);
	}
	if(rstatus) {
		sprintf(log_msg,"drms_getkey_xxxx(%s)",kw);
		lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_DATA,rstatus, log_msg, log_msg_code, results->opf);
		write_mini_output(PROCSTAT_NO_LF_DB_READ_PB,record_in,record_out,tbf,results);
		*status=ERR_DRMS_READ_MISSING_DATA;   
		return(ERR_DRMS_READ_MISSING_DATA);   
	}
	switch( typ ) 
	{
	    case 1: 
			fits_update_key(outfptr, TSTRING, kw, kw_c, kw_txt, &rstatus);
			break;
	    case 2: 
			fits_update_key(outfptr, TINT, kw, &kw_i, kw_txt, &rstatus);
			break;
	    case 3: 
			fits_update_key(outfptr, TFLOAT, kw, &kw_f, kw_txt, &rstatus);
			break;
	    case 4: 
			fits_update_key(outfptr, TDOUBLE, kw, &kw_d, kw_txt, &rstatus);
			break;
	}

	if ( rstatus != 0)
	{
		lf_logmsg4fitsio(log_msg, log_msg_code,kw,fsn,rstatus,results->opf);			
		write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,tbf,results);
		*status=ERR_DRMS_READ_MISSING_DATA;   
		return ERR_FITSIO;
	}
return(0);
}

int do_one_limbfit(unsigned int fsn, DRMS_Record_t *record_in,DRMS_Record_t *record_out, 
					LIMBFIT_INPUT *input, LIMBFIT_OUTPUT *results, LIMBFIT_IO_PUT *ios, int *status)
{
	static char *log_msg_code="do_one_limbfit";
	static char *series_name="limbfit data";
	static char *file_type="full FITS";
	static char *proc_stat=PROCSTAT_OK;

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
		// add a test if these values are ok, otherwise terminate this observation
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
			
			//---------------------------------------
			//		1) get all keywords
			//---------------------------------------
			drms_copykeys(record_out, record_in, 1, kDRMS_KeyClass_All); 	
	
			// change DATE to our proc_date, do not keep level 1 DATE:
			drms_keyword_setdate(record_out); 
			//---------------------------------------
			//		2) set all keywords
			//---------------------------------------
// need to test status...
			rstatus=drms_setkey_string(record_out, "SERIESCN", 	series_name);
			rstatus=drms_setkey_string(record_out, "INSERIES", 	results->dsin);
			rstatus=drms_setkey_string(record_out, "CODENAME",	results->code_name);
			rstatus=drms_setkey_string(record_out, "CODEVERS",	results->code_version);
			rstatus=drms_setkey_string(record_out, "CODEDATE",	results->code_date);
			rstatus=drms_setkey_float (record_out, "X_LFS",		results->cenx);
			rstatus=drms_setkey_float (record_out, "Y_LFS",		results->ceny);
			rstatus=drms_setkey_double(record_out, "R_LFS",		results->radius);
			rstatus=drms_setkey_int   (record_out, "QUAL_LFS",	results->quality);
			rstatus=drms_setkey_int   (record_out, "ERRO_LFS",	results->error1);
			rstatus=drms_setkey_int   (record_out, "ERRO_FIT",	results->error2);
			rstatus=drms_setkey_int   (record_out, "NB_FBINS",	results->nb_fbins);
			rstatus=drms_setkey_int   (record_out, "FLDFR",		results->fldfr);
			rstatus=drms_setkey_double(record_out, "MAX_LIMB",	results->max_limb);
			rstatus=drms_setkey_double(record_out, "CMEAN",		results->cmean);
			rstatus=drms_setkey_int   (record_out, "ANN_WD",	results->ann_wd);
			rstatus=drms_setkey_int   (record_out, "MXSZANNV",	results->mxszannv);
			rstatus=drms_setkey_int   (record_out, "NB_LDF",	results->nb_ldf);
			rstatus=drms_setkey_int   (record_out, "NB_RDB",	results->nb_rdb);
			rstatus=drms_setkey_int   (record_out, "NB_ABB",	results->nb_abb);
			rstatus=drms_setkey_double(record_out, "UP_LIMIT",	results->up_limit);
			rstatus=drms_setkey_double(record_out, "LO_LIMIT",	results->lo_limit);
			rstatus=drms_setkey_double(record_out, "INC_X",		results->inc_x);
			rstatus=drms_setkey_double(record_out, "INC_Y",		results->inc_y);
			rstatus=drms_setkey_int   (record_out, "NFITPNTS",	results->nfitpnts);
			rstatus=drms_setkey_int   (record_out, "NB_ITER",	results->nb_iter);
			rstatus=drms_setkey_int   (record_out, "CEN_CALC",	input->cc);
			rstatus=drms_setkey_double(record_out, "AHI",		results->ahi);
			rstatus=drms_setkey_string(record_out, "COMMENTS",	results->comment);			
			rstatus=drms_setkey_string(record_out, "BLD_VERS",  results->bld_vers);
			//---------------------------------------
			//		3) write the generic segment 
			//			(as the full FITS file)
			//---------------------------------------
// something left...: in case of problem, filenameout not closed and therefore not removed...

			long nelements, firstrow, firstelem;
			static char filenameout[128];
			sprintf(filenameout,"%s/%u_full.fits",results->tmp_dir,fsn);
		
			firstrow=1;
			firstelem=1;
			*status=0;
			rstatus=0;
			remove(filenameout); 
			fitsfile *outfptr;
			if (fits_create_file(&outfptr, filenameout, &rstatus)) 
			{
				fits_report_error(results->opf, rstatus); 
				lf_logmsg4fitsio(log_msg, log_msg_code, "fits_create_file()",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,2,results);
				return ERR_FITSIO;
			}

			long l_naxes[2];
			l_naxes[0]=results->fits_ldfs_naxis1;
			l_naxes[1]=results->fits_ldfs_naxis2;
			nelements=l_naxes[0]*l_naxes[1];
			int nax=2;
			if ( fits_create_img(outfptr, FLOAT_IMG, nax, l_naxes, &rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "fits_create_img()",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,2,results);
				return ERR_FITSIO;
			}
			if ( fits_write_img(outfptr, TFLOAT, 1, nelements, results->fits_ldfs_data, &rstatus) ) 
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "fits_write_img()",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,2,results);
				return ERR_FITSIO;
			}	
//			drms_free_array (data_array1); // NOTE: because of that, from now tbf MUST BE EQUAL to 1

		//add more KWs
			if(get_set_kw(1,"ORIGIN","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			
			TIME tobs;
			static char s_tobs[50];
			tobs = drms_getkey_time(record_in, "DATE__OBS", &rstatus);
			sprint_time(s_tobs, tobs, "UTC", 0);
			if ( fits_update_key(outfptr, TSTRING, "DATE-OBS", &s_tobs, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(DATE__OBS)",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}
			tobs = drms_getkey_time(record_in, "T_OBS", &rstatus);
			sprint_time(s_tobs, tobs, "UTC", 0);
			if ( fits_update_key(outfptr, TSTRING, "T_OBS", s_tobs, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(T_OBS)",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}
				
			if(get_set_kw(2,"CAMERA","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if(get_set_kw(4,"EXPTIME","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if(get_set_kw(3,"EXPSDEV","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if(get_set_kw(2,"WAVELNTH","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if(get_set_kw(1,"WAVEUNIT","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if ( fits_update_key(outfptr, TINT, "FSN", &fsn, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(FSN)",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}
			if(get_set_kw(2,"FID","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if(get_set_kw(2,"DATAMIN","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if(get_set_kw(2,"DATAMAX","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if(get_set_kw(2,"DATAMEDN","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if(get_set_kw(3,"DATAMEAN","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if(get_set_kw(3,"DATARMS","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if(get_set_kw(1,"FLAT_REC","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if(get_set_kw(1,"CTYPE1","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if(get_set_kw(1,"CUNIT1","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if(get_set_kw(3,"CRVAL1","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if(get_set_kw(3,"CDELT1","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if(get_set_kw(3,"CRPIX1","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if(get_set_kw(1,"CTYPE2","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if(get_set_kw(1,"CUNIT2","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if(get_set_kw(3,"CRVAL2","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if(get_set_kw(3,"CDELT2","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if(get_set_kw(3,"CRPIX2","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if(get_set_kw(3,"CROTA2","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if(get_set_kw(3,"R_SUN","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if(get_set_kw(1,"MPO_REC","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if(get_set_kw(3,"INST_ROT","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);			
			if (isnan(input->ix)) { input->ix=0.0;}
			if ( fits_update_key(outfptr, TDOUBLE, "X0_LF", &input->ix, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(X0_LF)",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}
			if (isnan(input->iy)) { input->iy=0.0;}
			if ( fits_update_key(outfptr, TDOUBLE, "Y0_LF", &input->iy, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(Y0_LF)",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}
			if (isnan(input->ir)) { input->ir=0.0;}
			if ( fits_update_key(outfptr, TDOUBLE, "RSUN_LF", &input->ir, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(RSUN_LF)",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}
			if(get_set_kw(1,"ASD_REC","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if(get_set_kw(3,"SAT_ROT","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if(get_set_kw(1,"ORB_REC","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if(get_set_kw(4,"DSUN_OBS","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if(get_set_kw(4,"RSUN_OBS","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if(get_set_kw(4,"HAEX_OBS","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if(get_set_kw(4,"HAEY_OBS","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if(get_set_kw(4,"HAEZ_OBS","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if(get_set_kw(2,"HPLTID","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if(get_set_kw(2,"HFTSACID","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if(get_set_kw(2,"HWLTID","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if(get_set_kw(2,"HCFTID","",fsn,record_in,record_out,outfptr,1,results,status)) return(*status);
			if ( fits_update_key(outfptr, TSTRING, "SERIESCN", series_name, "Series Name",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(SERIESCN)",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TSTRING, "FILETYPE", file_type, "File content",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(FILETYPE)",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TSTRING, "CODENAME", results->code_name, "Code Name",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(CODENAME)",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TSTRING, "CODEVERS", results->code_version, "Code Version",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(CODEVERS)",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TSTRING, "CODEDATE", results->code_date, "Code Date",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(CODEDATE)",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TSTRING, "PROCSTAT", proc_stat, "Processing rstatus code",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(PROCSTAT)",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}
		
			tobs=drms_keyword_getdate(record_out);		
			sprint_time(s_tobs, tobs, "UTC", 0);
			if ( fits_update_key(outfptr, TSTRING, "DATE", s_tobs, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(DATE)",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TFLOAT, "X_LFS", &results->cenx, "Center Position X",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(X_LFS)",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TFLOAT, "Y_LFS", &results->ceny, "Center Position Y",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(Y_LFS)",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TDOUBLE, "R_LFS", &results->radius, "Solar Radius",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(R_LFS)",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TDOUBLE, "CMEAN", &results->cmean, "Mean value of the 500x500 pixels center",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(CMEAN)",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TDOUBLE, "MAX_LIMB", &results->cmean, "Maximum value on the limb",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(MAX_LIMB)",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "QUAL_LFS", &results->quality, "Processing quality indicator",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(QUAL_LFS)",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "ERRO_LFS", &results->error1, "Limb.f Processing error code",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(ERRO_LFS)",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "ERRO_FIT", &results->error2, "Fitting Processing error code",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(ERRO_FIT)",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "FLDFR", &results->fldfr, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(FLDFR)",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "NB_FBINS", &results->nb_fbins, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(NB_FBINS)",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "ANN_WD", &results->ann_wd, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(ANN_WD)",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "MXSZANNV", &results->mxszannv, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(MXSZANNV)",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "NB_LDF", &results->nb_ldf, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(NB_LDF)",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "NB_RDB", &results->nb_rdb, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(NB_RDB)",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "NB_ABB", &results->nb_abb, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(NB_ABB)",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TDOUBLE, "UP_LIMIT", &results->up_limit, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(UP_LIMIT)",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TDOUBLE, "LO_LIMIT", &results->lo_limit, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(LO_LIMIT)",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TDOUBLE, "INC_X", &results->inc_x, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(INC_X)",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TDOUBLE, "INC_Y", &results->inc_y, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(INC_Y)",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "NFITPNTS", &results->nfitpnts, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(NFITPNTS)",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "NB_ITER", &results->nb_iter, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(NB_ITER)",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}

			rstatus=0;
			int i;			
		
			// AB ---EXT#1---------------------------------------------------------------------
			if (results->nb_iter==1)
			{
				char *tunitAB[] = { "","" };
				char *ttypeAB[] = { "Alpha", "Beta"};
				char *tformAB[] = { "1E","1E" };
				char extnameAB[]="LDF info / Alpha&Beta";
				if ( fits_create_tbl( outfptr, BINARY_TBL, results->fits_ab_nrows, results->fits_ab_tfields, 
													ttypeAB, tformAB, tunitAB, extnameAB, &rstatus) )
				{
					fits_report_error(results->opf, rstatus); 
					return ERR_FITSIO;
				}
			}
			else
			{
				char *tunitAB[] = { "","",""};
				char *ttypeAB[] = { "Alpha", "Beta","Beta Prev"};
				char *tformAB[] = { "1E","1E","1E" };
				char extnameAB[]="LDF info / Alpha&Beta[&Beta Prev]";
				if ( fits_create_tbl( outfptr, BINARY_TBL, results->fits_ab_nrows, results->fits_ab_tfields, 
													ttypeAB, tformAB, tunitAB, extnameAB, &rstatus) )
				{
					fits_report_error(results->opf, rstatus); 
					return ERR_FITSIO;
				}			
			}
			for (i=1;i<=results->fits_ab_tfields;i++)
			{
				if (fits_write_col(outfptr, TFLOAT, i, firstrow, firstelem, results->fits_ab_nrows, 
												results->fits_alpha_beta+(i-1)*results->fits_ab_nrows, &rstatus))
				{
					lf_logmsg4fitsio(log_msg, log_msg_code, "fits_write_col(AB)",fsn,rstatus,results->opf);			
					write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
					return ERR_FITSIO;
				}
			}

			// All Params ---EXT#2-------------------------------------------------------------------------
			char extnameA[] = "LDF info / As / Es / Radius / IP";
			char *tunitAE[] = { "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "" };
			char *tformAE[] = { "1D", "1D", "1D", "1D", "1D", "1D", "1D", "1D", "1D", "1D", "1D", "1D", "1D", "1D"  };
			char *ttypeA[]  = { "A0", "A1", "A2", "A3", "A4", "A5", "E0", "E1", "E2", "E3", "E4", "E5", "Radius", "IP1"};
		
			if ( fits_create_tbl(outfptr, BINARY_TBL, results->fits_params_nrows, results->fits_params_tfields, 
												ttypeA, tformAE,tunitAE, extnameA, &rstatus) )
			{
					fits_report_error(results->opf, rstatus); 
					return ERR_FITSIO;
			}

			for (i=1;i<=results->fits_params_tfields;i++)
			{
				if (fits_write_col(outfptr, TDOUBLE, i, firstrow, firstelem, results->fits_params_nrows, 
												results->fits_params+(i-1)*results->fits_params_nrows, &rstatus))
				{
					lf_logmsg4fitsio(log_msg, log_msg_code, "fits_write_col(PARAMS)",fsn,rstatus,results->opf);			
					write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
					return ERR_FITSIO;
				}
			}	
			
			if (input->fldf==1)
			{
				// Full LDF ---EXT#3-------------------------------------------------------------------------
				char *tunitF[] = { "","" };
				char *ttypeF[] = { "Intensity", "Radius"};
				char *tformF[] = { "1E","1E" };
				char extnameF[]="Full LDF ";
				if ( fits_create_tbl( outfptr, BINARY_TBL, results->fits_fldfs_nrows, results->fits_fldfs_tfields, 
														ttypeF, tformF, tunitF, extnameF, &rstatus) )
				{
					fits_report_error(results->opf, rstatus); 
					return ERR_FITSIO;
				}
				for (i=1;i<=results->fits_fldfs_tfields;i++)
				{
					if (fits_write_col(outfptr, TFLOAT, i, firstrow, firstelem, results->fits_fldfs_nrows, 
													results->fits_fulldfs+(i-1)*results->fits_fldfs_nrows, &rstatus))
					{
						lf_logmsg4fitsio(log_msg, log_msg_code, "fits_write_col(FLDF)",fsn,rstatus,results->opf);			
						write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
						return ERR_FITSIO;
					}
				}
			}
			//---------------------------------------
			//		5) close & free
			//---------------------------------------
			if ( fits_update_key(outfptr, TSTRING, "PROCSTAT", proc_stat, "Processing rstatus code",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(PROCSTAT)",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}

			// Close file & attach to DB & remove
			if ( fits_close_file(outfptr, &rstatus) ) 
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "fits_close_file()",fsn,rstatus,results->opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,1,results);
				return ERR_FITSIO;
			}						

			//warning change segment number !!!
			segment_out = drms_segment_lookupnum (record_out, 0); //WARNING change this number if the number of segment changes...

			rstatus=drms_segment_write_from_file(segment_out,filenameout);
			if(rstatus) {
					sprintf(log_msg,"can't drms_segment_write_from_file() for fsn# %u", fsn);
					lf_logmsg("ERROR", "DRMS", ERR_DRMS_WRITE, rstatus, log_msg, log_msg_code, results->opf);
					return(ERR_DRMS_WRITE);
			}

			rstatus=drms_setkey_string(record_out, "PROCSTAT",	PROCSTAT_OK);
			if (!rstatus) remove(filenameout); // test if this is ok too!!!
			else 
				{
					sprintf(log_msg,"can't drms_setkey_string(PROCSTAT) for fsn# %u", fsn);
					lf_logmsg("ERROR", "DRMS", ERR_DRMS_WRITE, rstatus, log_msg, log_msg_code, results->opf);
					return(ERR_DRMS_WRITE);
				}
			if (results->debug) printf("end write in the db %d\n",rstatus);
				
			free(results->fits_ldfs_data);
			free(results->fits_alpha_beta);
			free(results->fits_params);
			if (results->fldfr<3) free(results->fits_fulldfs);
			
			if (results->debug) printf("done...%u\n",fsn);							
			lf_logmsg("INFO", "APP", 0, 0, "End writing in the DB", log_msg_code, results->opf);
		}
		else 
		{
			sprintf(log_msg,"NO LDFS file created (ifail >0) or retcode <0 for fsn# %u", fsn);
			lf_logmsg("WARNING", "APP", ERR_LIMBFIT_FAILED, lf_retcode, log_msg, log_msg_code, results->opf);			
			write_mini_output(PROCSTAT_NO_LF_FAILED,record_in,record_out,0,results);
			// is this the right retcode? or should I make a test???
			
			// free 
			if (results->error1 == 0)
			{
				free(results->fits_ldfs_data);
				free(results->fits_alpha_beta);
				free(results->fits_params);
				if (results->fldfr<3) free(results->fits_fulldfs);
			}
			return(0);
		}
	}

return 0;
}
//pb here: quand this is called then no free are executed!!! 
int	write_mini_output(char * errcode, DRMS_Record_t *record_in,DRMS_Record_t *record_out, 
			int tbf, LIMBFIT_OUTPUT *results)
{
		//int status;
		char *log_msg_code="write_mini_output";
		lf_logmsg("INFO", "APP", 0, 0, "Writing only header", log_msg_code, results->opf);			
		// write only the header of segment in the DB
		drms_copykeys(record_out, record_in, 1, kDRMS_KeyClass_All); 	
		// change DATE to our proc_date, do not keep level 1 DATE:
		drms_keyword_setdate(record_out); 

// need to test status...
		drms_setkey_string(record_out, "SERIESCN", "Limbfit data");
		drms_setkey_string(record_out, "CODENAME",	results->code_name);
		drms_setkey_string(record_out, "CODEVERS",	results->code_version);
		drms_setkey_string(record_out, "CODEDATE",	results->code_date);
		drms_setkey_string(record_out, "INSERIES",	results->dsin);
		drms_setkey_string(record_out, "COMMENTS",	results->comment);
		drms_setkey_string(record_out, "PROCSTAT",	errcode);
		drms_setkey_int   (record_out, "QUAL_LFS",	results->quality);
		drms_setkey_int   (record_out, "ERRO_LFS",	results->error1);
		drms_setkey_int   (record_out, "ERRO_FIT",	results->error2);
		drms_setkey_int   (record_out, "ANN_WD",	results->ann_wd);
		drms_setkey_int   (record_out, "MXSZANNV",	results->mxszannv);
		drms_setkey_int   (record_out, "NB_LDF",	results->nb_ldf);
		drms_setkey_int   (record_out, "NB_RDB",	results->nb_rdb);
		drms_setkey_int   (record_out, "NB_ABB",	results->nb_abb);
		drms_setkey_double(record_out, "UP_LIMIT",	results->up_limit);
		drms_setkey_double(record_out, "LO_LIMIT",	results->lo_limit);
		drms_setkey_double(record_out, "INC_X",		results->inc_x);
		drms_setkey_double(record_out, "INC_Y",		results->inc_y);
		drms_setkey_int   (record_out, "NB_ITER",	results->nb_iter);
		drms_setkey_int   (record_out, "CEN_CALC",	results->cc);
		drms_setkey_double(record_out, "AHI",		results->ahi);		
		drms_setkey_int   (record_out, "FLDFR",		results->fldfr);
		drms_setkey_int   (record_out, "NFITPNTS",	results->nfitpnts);
		drms_setkey_string(record_out, "BLD_VERS", 	results->bld_vers);
		if (tbf >=1)
		{
			free(results->fits_alpha_beta);
			free(results->fits_params);
			if (results->fldfr<3) free(results->fits_fulldfs);
			if (tbf ==2) free(results->fits_ldfs_data);	
		}
return(0);
}
