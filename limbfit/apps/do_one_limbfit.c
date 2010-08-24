/* 
	I.Scholl 

		Second level wrapper for the limbfit code
		Can be used as a callable module (to be plugged into lev0lev1 module or by the main wrapper)
			calls: limbfit() and reformat output data to be rewritten in the DB
			
	
	#define CODE_NAME 		"limbfit"
	#define CODE_VERSION 	"V1.4r0" 
	#define CODE_DATE 		"Tue Aug 24 10:19:25 HST 2010" 
*/

#include "limbfit.h"

void lf_logmsg4fitsio(char *log_msg,char *log_msg_code,char *kw,unsigned int fsn,int status, FILE *opf)
{	
	sprintf(log_msg,"can't fits_update_key[%s] for fsn# %u", kw,fsn);
	lf_logmsg("ERROR", "FITSIO", ERR_FITSIO, status, log_msg, log_msg_code, opf);			
	fits_report_error(opf, status); 
}

int get_set_kw(int typ, char *kw, char *kw_txt, unsigned int fsn, DRMS_Record_t *record_in,DRMS_Record_t *record_out, fitsfile *outfptr, FILE *opf, int debug, int *status)
{
	float kw_f;
	double kw_d;
	int kw_i;
	char *kw_c;
	int rstatus;
	char log_msg[120];
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
		lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_DATA,rstatus, log_msg, log_msg_code, opf);
		write_mini_output(fsn,PROCSTAT_NO_LF_DB_READ_PB,record_in,record_out,opf,debug); 
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
		lf_logmsg4fitsio(log_msg, log_msg_code,kw,fsn,rstatus,opf);			
		write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
		*status=ERR_DRMS_READ_MISSING_DATA;   
		return ERR_FITSIO;
	}
return(0);
}

int do_one_limbfit(unsigned int fsn, DRMS_Record_t *record_in,DRMS_Record_t *record_out,FILE *opf, int debug, int *status)
{
	static char *log_msg_code="do_one_limbfit";
	static char *series_name="limbfit data";
	static char *file_type="full FITS";
	static char *proc_stat=PROCSTAT_OK;

	static LIMBFIT_INPUT limbfit_vars ;
	static LIMBFIT_INPUT *lfv = &limbfit_vars;
	static LIMBFIT_OUTPUT limbfit_res ;
	static LIMBFIT_OUTPUT *lfr = &limbfit_res;
 	int lf_retcode;
	int rstatus;
	char log_msg[120];
	
	int seg_cnt;
	DRMS_Segment_t *segment_in;
	DRMS_Segment_t *segment_out;
	static DRMS_Array_t *img;

	seg_cnt = drms_record_numsegments (record_in);
	if (seg_cnt < 1) 
	{
		sprintf(log_msg,"no segments in selected record_in");
		lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_DATA, 0, log_msg, log_msg_code,opf);
		write_mini_output(fsn,PROCSTAT_NO_LF_DB_READ_PB,record_in,record_out,opf,debug); 
		return(ERR_DRMS_READ_MISSING_DATA);   
	}
	segment_in = drms_segment_lookupnum(record_in, 0);
	img = drms_segment_read(segment_in, DRMS_TYPE_FLOAT, &rstatus);
	if(!img) {
		sprintf(log_msg,"Can't do drms_segment_read() status=%d", rstatus);
		lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_DATA, 0, log_msg, log_msg_code, opf);
		write_mini_output(fsn,PROCSTAT_NO_LF_DB_READ_PB,record_in,record_out,opf,debug); 
		return(ERR_DRMS_READ_MISSING_DATA);   
	} 
	else
	{
		// add a test if these values are ok, otherwise terminate this observation
		lfv->ix = drms_getkey_float(record_in, "X0_LF", &rstatus);
		if(rstatus) {
			lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_DATA,rstatus, "drms_getkey_float(X0_LF)", log_msg_code, opf);
			write_mini_output(fsn,PROCSTAT_NO_LF_DB_READ_PB,record_in,record_out,opf,debug); 
			return(ERR_DRMS_READ_MISSING_DATA);   
		}
		if(isnan(lfv->ix)) {
			lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_XYR,rstatus, "X0_LF missing", log_msg_code, opf);
			write_mini_output(fsn,PROCSTAT_NO_LF_DB_READ_PB,record_in,record_out,opf,debug); 
			return(ERR_DRMS_READ_MISSING_XYR);   
		}
		lfv->iy = drms_getkey_float(record_in, "Y0_LF", &rstatus);
		if(rstatus) {
			lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_DATA,rstatus, "drms_getkey_float(Y0_LF)", log_msg_code, opf);
			write_mini_output(fsn,PROCSTAT_NO_LF_DB_READ_PB,record_in,record_out,opf,debug); 
			return(ERR_DRMS_READ_MISSING_DATA);   
		}
		if(isnan(lfv->iy)) {
			lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_XYR,rstatus, "Y0_LF missing", log_msg_code, opf);
			write_mini_output(fsn,PROCSTAT_NO_LF_DB_READ_PB,record_in,record_out,opf,debug); 
			return(ERR_DRMS_READ_MISSING_XYR);   
		}
		lfv->ir = drms_getkey_float(record_in, "RSUN_LF", &rstatus);
		  if(rstatus || (isnan(lfv->ir)) ) {
			lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_DATA,rstatus, "drms_getkey_float(RSUN_LF)", log_msg_code, opf);
			write_mini_output(fsn,PROCSTAT_NO_LF_DB_READ_PB,record_in,record_out,opf,debug); 
			return(ERR_DRMS_READ_MISSING_DATA);   
		}
		if(isnan(lfv->ir)) {
			lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_XYR,rstatus, "RSUN_LF missing", log_msg_code, opf);
			write_mini_output(fsn,PROCSTAT_NO_LF_DB_READ_PB,record_in,record_out,opf,debug); 
			return(ERR_DRMS_READ_MISSING_XYR);   
		}
		lfv->img_sz0=img->axis[0];
		lfv->img_sz1=img->axis[1];
		lfv->data=img->data;

		if (debug) 
		{
			sprintf(log_msg,"X0_LF %f", lfv->ix);
			lf_logmsg("DEBUG", "INFO", 0, 0, log_msg, log_msg_code, opf);
			sprintf(log_msg,"Y0_LF %f", lfv->iy);
			lf_logmsg("DEBUG", "INFO", 0, 0, log_msg, log_msg_code, opf);
			sprintf(log_msg,"RSUN_LF %f", lfv->ir);
			lf_logmsg("DEBUG", "INFO", 0, 0, log_msg, log_msg_code, opf);
		}
		
		lf_retcode=limbfit(lfv,lfr,opf,debug);
		free(img->data);
		

		if (debug) 
		{
			sprintf(log_msg,"XC %f", lfr->cenx);
			lf_logmsg("DEBUG", "INFO", 0, 0, log_msg, log_msg_code, opf);
			sprintf(log_msg,"YC %f", lfr->ceny);
			lf_logmsg("DEBUG", "INFO", 0, 0, log_msg, log_msg_code, opf);
			sprintf(log_msg,"R %f", lfr->radius);
			lf_logmsg("DEBUG", "INFO", 0, 0, log_msg, log_msg_code, opf);
			sprintf(log_msg,"Returned code %d", lf_retcode);
			lf_logmsg("DEBUG", "INFO", 0, 0, log_msg, log_msg_code, opf);
			sprintf(log_msg,"err1 %d", lfr->error1);
			lf_logmsg("DEBUG", "INFO", 0, 0, log_msg, log_msg_code, opf);
			sprintf(log_msg,"err2 %d", lfr->error2);
			lf_logmsg("DEBUG", "INFO", 0, 0, log_msg, log_msg_code, opf);
		}
		
		// >= 0 => all
		// ERR_LIMBFIT_FIT_FAILED => all
		// ERR_LIMBFIT_FAILED => mini
		
		if (lf_retcode >= 0 || lf_retcode == ERR_LIMBFIT_FIT_FAILED)
		{	
			if (debug) printf("write in the db\n");
			//**********************************************************************
			//	Write into the DB
			//**********************************************************************
			int status;
			
			//---------------------------------------
			//		1) get all keywords
			//---------------------------------------
			drms_copykeys(record_out, record_in, 1, kDRMS_KeyClass_All); 	
	
			// change DATE to our proc_date, do not keep level 1 DATE:
			drms_keyword_setdate(record_out); 
			//---------------------------------------
			//		2) set all keywords
			//---------------------------------------
			char code_name[10]=CODE_NAME;
			char code_version[10]=CODE_VERSION;
			char code_date[30]=CODE_DATE;
			int num_ext=2;//to change with full_ldfs
			
// need to test status...
			status=drms_setkey_string(record_out, "SERIESCN", 	series_name);
			status=drms_setkey_string(record_out, "CODENAME",	code_name);
			status=drms_setkey_string(record_out, "CODEVERS",	code_version);
			status=drms_setkey_string(record_out, "CODEDATE",	code_date);
			status=drms_setkey_int	 (record_out, "NUM_EXT",	num_ext);  
			status=drms_setkey_float (record_out, "X_LFS",		lfr->cenx);
			status=drms_setkey_float (record_out, "Y_LFS",		lfr->ceny);
			status=drms_setkey_double(record_out, "R_LFS",		lfr->radius);
			status=drms_setkey_int   (record_out, "QUAL_LFS",	lfr->quality);
			status=drms_setkey_int   (record_out, "ERRO_LFS",	lfr->error1);
			status=drms_setkey_int   (record_out, "ERRO_FIT",	lfr->error2);
			status=drms_setkey_double(record_out, "MAX_LIMB",	lfr->max_limb);
			status=drms_setkey_double(record_out, "CMEAN",		lfr->cmean);
			status=drms_setkey_int   (record_out, "ANN_WD",		lfr->ann_wd);
			status=drms_setkey_int   (record_out, "MXSZANNV",	lfr->mxszannv);
			status=drms_setkey_int   (record_out, "NB_LDF",		lfr->nb_ldf);
			status=drms_setkey_int   (record_out, "NB_RDB",		lfr->nb_rdb);
			status=drms_setkey_int   (record_out, "NB_ABB",		lfr->nb_abb);
			status=drms_setkey_double(record_out, "UP_LIMIT",	lfr->up_limit);
			status=drms_setkey_double(record_out, "LO_LIMIT",	lfr->lo_limit);
			status=drms_setkey_double(record_out, "INC_X",		lfr->inc_x);
			status=drms_setkey_double(record_out, "INC_Y",		lfr->inc_y);
			status=drms_setkey_int   (record_out, "NFITPNTS",	lfr->nfitpnts);
			status=drms_setkey_int   (record_out, "NB_ITER",	lfr->nb_iter);
			status=drms_setkey_int   (record_out, "SKIPGC",		lfr->skipgc);

			//---------------------------------------
			//		3) write segments
			//---------------------------------------
			int i_naxes[2];
			i_naxes[0]=lfr->fits_ldfs_naxis1;
			i_naxes[1]=lfr->fits_ldfs_naxis2;
			//ldfs
			DRMS_Array_t *data_array1;
			segment_out = drms_segment_lookupnum (record_out, 0);
			drms_setkey_string(record_out, "FILETYPE[0]", "LDFS");
			data_array1 = drms_array_create (DRMS_TYPE_FLOAT, 2, i_naxes, (void *)lfr->fits_ldfs_data, &status);
			data_array1->bscale = 1.0;
			data_array1->bzero = 0.0;
			status = drms_segment_write (segment_out, data_array1, 0);
			if (status) 
			{
				sprintf(log_msg,"can't create segment 0 (LDF) for fsn# %u in series", fsn);
				lf_logmsg("ERROR", "DRMS", ERR_DRMS_WRITE, status, log_msg, log_msg_code, opf);			
				close_on_error(	record_in,record_out,data_array1,opf);
				return ERR_DRMS_WRITE;
			}

			//full LDF
			DRMS_Array_t *data_array2;
			segment_out = drms_segment_lookupnum (record_out, 1);
			i_naxes[0]=lfr->fits_fldfs_tfields;
			i_naxes[1]=lfr->fits_fldfs_nrows; 					// to check row/col
			drms_setkey_string(record_out, "FILETYPE[1]", "FULL_LDFS");
			data_array2 = drms_array_create (DRMS_TYPE_FLOAT, 2, i_naxes, (void *)lfr->fits_fulldfs, &status) ;
			data_array2->bscale = 1.0;
			data_array2->bzero = 0.0;
			status = drms_segment_write (segment_out, data_array2, 1);
			if (status) 
			{
				sprintf(log_msg,"can't create segment 0 (fldf) for fsn# %u in series", fsn);
				lf_logmsg("ERROR", "DRMS", ERR_DRMS_WRITE, status, log_msg, log_msg_code, opf);			
				close_on_error(	record_in,record_out,data_array2,opf);
				return ERR_DRMS_WRITE;
			}
		
			//alpha_betas
			DRMS_Array_t *data_array3;
			segment_out = drms_segment_lookupnum (record_out, 2);
			i_naxes[0]=lfr->fits_ab_tfields;
			i_naxes[1]=lfr->fits_ab_nrows;
			drms_setkey_string(record_out, "FILETYPE[2]", "AB");
			data_array3 = drms_array_create (DRMS_TYPE_FLOAT, 2 ,i_naxes, (void *)lfr->fits_alpha_beta1, &status);
			data_array3->bscale = 1.0;
			data_array3->bzero = 0.0;
			status = drms_segment_write (segment_out, data_array3, 2);
			if (status) 
			{
				sprintf(log_msg,"can't create segment 0 (AB) for fsn# %u in series", fsn);
				lf_logmsg("ERROR", "DRMS", ERR_DRMS_WRITE, status, log_msg, log_msg_code, opf);			
				close_on_error(	record_in,record_out,data_array3,opf);
				return ERR_DRMS_WRITE;
			}
	
			//PARAMS
			DRMS_Array_t *data_array4;
			segment_out = drms_segment_lookupnum (record_out, 3);
			i_naxes[0]=lfr->fits_params_tfields;
			i_naxes[1]=lfr->fits_params_nrows;
			drms_setkey_string(record_out, "FILETYPE[3]", "PARAMS");
			data_array4 = drms_array_create (DRMS_TYPE_DOUBLE, 2, i_naxes, (void *)lfr->fits_params1, &status) ;
			data_array4->bscale = 1.0;
			data_array4->bzero = 0.0;
			status = drms_segment_write (segment_out, data_array4, 3);
			if (status) 
			{
				sprintf(log_msg,"can't create segment 0 (PARAMS) for fsn# %u in series", fsn);
				lf_logmsg("ERROR", "DRMS", ERR_DRMS_WRITE, status, log_msg, log_msg_code, opf);			
				close_on_error(	record_in,record_out,data_array4,opf);
				return ERR_DRMS_WRITE;
			}
	
			//---------------------------------------
			//		4) write the generic segment 
			//			(as the full FITS file)
			//---------------------------------------
// something left...: in case of problem, filenameout not closed and therefore not removed...

			long nelements, firstrow, firstelem;
			static char filenameout[128];
			//sprintf(filenameout,"%s%u_full.fits",TMP_DIR,fsn);
			sprintf(filenameout,"./%u_full.fits",fsn);
		
			firstrow=1;
			firstelem=1;
			status=0;
		
			remove(filenameout); 
			fitsfile *outfptr;
			if (fits_create_file(&outfptr, filenameout, &status)) 
			{
				fits_report_error(opf, status); 
				lf_logmsg4fitsio(log_msg, log_msg_code, "fits_create_file()",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}

			long l_naxes[2];
			l_naxes[0]=lfr->fits_ldfs_naxis1;
			l_naxes[1]=lfr->fits_ldfs_naxis2;
			nelements=l_naxes[0]*l_naxes[1];
			int nax=2;
			if ( fits_create_img(outfptr, FLOAT_IMG, nax, l_naxes, &status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "fits_create_img()",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}
			//add more KWs
			if(get_set_kw(1,"ORIGIN","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(1,"TELESCOP","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(1,"INSTRUME","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status); 
			
			TIME tobs;
			static char s_tobs[50];
			tobs = drms_getkey_time(record_in, "DATE__OBS", &rstatus);
			sprint_time(s_tobs, tobs, "UTC", 0);
			if ( fits_update_key(outfptr, TSTRING, "DATE-OBS", &s_tobs, "",&status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(DATE__OBS)",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}
			tobs = drms_getkey_time(record_in, "T_OBS", &rstatus);
			sprint_time(s_tobs, tobs, "UTC", 0);
			if ( fits_update_key(outfptr, TSTRING, "T_OBS", s_tobs, "",&status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(T_OBS)",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}
				
			if(get_set_kw(2,"CAMERA","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(1,"IMG_TYPE","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(4,"EXPTIME","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(3,"EXPSDEV","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(2,"WAVELNTH","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(1,"WAVEUNIT","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if ( fits_update_key(outfptr, TINT, "FSN", &fsn, "",&status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(FSN)",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}
			if(get_set_kw(2,"FID","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(2,"QUALLEV0","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(2,"QUALITY","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(2,"DATAMIN","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(2,"DATAMAX","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(2,"DATAMEDN","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(3,"DATAMEAN","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(3,"DATARMS","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(2,"MISSVALS","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(1,"FLAT_REC","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(1,"CTYPE1","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(1,"CUNIT1","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(3,"CRVAL1","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(3,"CDELT1","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(3,"CRPIX1","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(1,"CTYPE2","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(1,"CUNIT2","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(3,"CRVAL2","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(3,"CDELT2","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(3,"CRPIX2","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(3,"CROTA2","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(3,"R_SUN","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(1,"MPO_REC","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(3,"INST_ROT","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);			
			if (isnan(lfv->ix)) { lfv->ix=0.0;}
			if ( fits_update_key(outfptr, TDOUBLE, "X0_LF", &lfv->ix, "",&status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(X0_LF)",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}
			if (isnan(lfv->iy)) { lfv->iy=0.0;}
			if ( fits_update_key(outfptr, TDOUBLE, "Y0_LF", &lfv->iy, "",&status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(Y0_LF)",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}
			if (isnan(lfv->ir)) { lfv->ir=0.0;}
			if ( fits_update_key(outfptr, TDOUBLE, "RSUN_LF", &lfv->ir, "",&status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(RSUN_LF)",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}
			if(get_set_kw(1,"ASD_REC","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(3,"SAT_ROT","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(1,"ORB_REC","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(4,"DSUN_OBS","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(4,"RSUN_OBS","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(4,"HCIEC_X","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(4,"HCIEC_Y","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(4,"HCIEC_Z","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(4,"HPLTID","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(1,"HWLTNSET","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(4,"HCFTID","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if(get_set_kw(4,"HFTSACID","",fsn,record_in,record_out,outfptr,opf,debug,&status)) return(status);
			if ( fits_update_key(outfptr, TSTRING, "SERIESCN", series_name, "Series Name",&status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(SERIESCN)",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TSTRING, "FILETYPE", file_type, "File content",&status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(FILETYPE)",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TSTRING, "CODENAME", code_name, "Code Name",&status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(CODENAME)",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TSTRING, "CODEVERS", code_version, "Code Version",&status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(CODEVERS)",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TSTRING, "CODEDATE", code_date, "Code Date",&status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(CODEDATE)",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TSTRING, "PROCSTAT", proc_stat, "Processing status code",&status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(PROCSTAT)",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}
		
			tobs=drms_keyword_getdate(record_out);		
			sprint_time(s_tobs, tobs, "UTC", 0);
			if ( fits_update_key(outfptr, TSTRING, "DATE", s_tobs, "",&status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(DATE)",fsn,status,opf);			
				write_mini_output(fsn,ERR_FITSIO,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "NUM_EXT", &num_ext, "Number of extensions",&status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(NUM_EXT)",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TFLOAT, "X_LFS", &lfr->cenx, "Center Position X",&status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(X_LFS)",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TFLOAT, "Y_LFS", &lfr->ceny, "Center Position Y",&status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(Y_LFS)",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TDOUBLE, "R_LFS", &lfr->radius, "Solar Radius",&status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(R_LFS)",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TDOUBLE, "CMEAN", &lfr->cmean, "Mean value of the 500x500 pixels center",&status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(CMEAN)",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TDOUBLE, "MAX_LIMB", &lfr->cmean, "Maximum value on the limb",&status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(MAX_LIMB)",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "QUAL_LFS", &lfr->quality, "Processing quality indicator",&status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(QUAL_LFS)",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "ERRO_LFS", &lfr->error1, "Limb.f Processing error code",&status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(ERRO_LFS)",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "ERRO_FIT", &lfr->error2, "Fitting Processing error code",&status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(ERRO_FIT)",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "ANN_WD", &lfr->ann_wd, "",&status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(ANN_WD)",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "MXSZANNV", &lfr->mxszannv, "",&status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(MXSZANNV)",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "NB_LDF", &lfr->nb_ldf, "",&status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(NB_LDF)",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "NB_RDB", &lfr->nb_rdb, "",&status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(NB_RDB)",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "NB_ABB", &lfr->nb_abb, "",&status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(NB_ABB)",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TDOUBLE, "UP_LIMIT", &lfr->up_limit, "",&status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(UP_LIMIT)",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TDOUBLE, "LO_LIMIT", &lfr->lo_limit, "",&status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(LO_LIMIT)",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TDOUBLE, "INC_X", &lfr->inc_x, "",&status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(INC_X)",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TDOUBLE, "INC_Y", &lfr->inc_y, "",&status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(INC_Y)",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "NFITPNTS", &lfr->nfitpnts, "",&status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(NFITPNTS)",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "NB_ITER", &lfr->nb_iter, "",&status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(NB_ITER)",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "SKIPGC", &lfr->skipgc, "",&status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(SKIPGC)",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}

			if ( fits_write_img(outfptr, TFLOAT, 1, nelements, lfr->fits_ldfs_data, &status) ) 
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "fits_write_img()",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}	

			status=0;
			int i;			
		
			// AB ---EXT#1---------------------------------------------------------------------
			char *tunitAB[] = { "","" };
			char *ttypeAB[] = { "Alpha", "Beta"};
			char *tformAB[] = { "1E","1E" };
			char extnameAB[]="LDF info / Alpha&Beta";
			if ( fits_create_tbl( outfptr, BINARY_TBL, lfr->fits_ab_nrows, lfr->fits_ab_tfields, 
													ttypeAB, tformAB, tunitAB, extnameAB, &status) )
			{
				fits_report_error(opf, status); 
				return ERR_FITSIO;
			}
			for (i=1;i<=lfr->fits_ab_tfields;i++)
			{
				if (fits_write_col(outfptr, TFLOAT, i, firstrow, firstelem, lfr->fits_ab_nrows, 
												lfr->fits_alpha_beta2+(i-1)*lfr->fits_ab_nrows, &status))
				{
					lf_logmsg4fitsio(log_msg, log_msg_code, "fits_write_col(AB)",fsn,status,opf);			
					write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
					return ERR_FITSIO;
				}
			}

			// All Params ---EXT#2-------------------------------------------------------------------------
			char extnameA[] = "LDF info / As / Es / Radius / IP";
			char *tunitAE[] = { "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "" };
			char *tformAE[] = { "1D", "1D", "1D", "1D", "1D", "1D", "1D", "1D", "1D", "1D", "1D", "1D", "1D", "1D"  };
			char *ttypeA[]  = { "A0", "A1", "A2", "A3", "A4", "A5", "E0", "E1", "E2", "E3", "E4", "E5", "Radius", "IP1"};
		
			if ( fits_create_tbl(outfptr, BINARY_TBL, lfr->fits_params_nrows, lfr->fits_params_tfields, 
												ttypeA, tformAE,tunitAE, extnameA, &status) )
			{
					fits_report_error(opf, status); 
					return ERR_FITSIO;
			}

			for (i=1;i<=lfr->fits_params_tfields;i++)
			{
				if (fits_write_col(outfptr, TDOUBLE, i, firstrow, firstelem, lfr->fits_params_nrows, 
												lfr->fits_params2+(i-1)*lfr->fits_params_nrows, &status))
				{
					lf_logmsg4fitsio(log_msg, log_msg_code, "fits_write_col(PARAMS)",fsn,status,opf);			
					write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
					return ERR_FITSIO;
				}
			}
	

			// Full LDF ---EXT#3-------------------------------------------------------------------------
			// one array with 2 columns
			

			//---------------------------------------
			//		5) close & free
			//---------------------------------------
			if ( fits_update_key(outfptr, TSTRING, "PROCSTAT", proc_stat, "Processing status code",&status) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(PROCSTAT)",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}

			// Close file & attach to DB & remove
			if ( fits_close_file(outfptr, &status) ) 
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "fits_close_file()",fsn,status,opf);			
				write_mini_output(fsn,PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,debug); 
				return ERR_FITSIO;
			}						
			//warning change segment number !!!
			segment_out = drms_segment_lookupnum (record_out, 4);
			status=drms_segment_write_from_file(segment_out,filenameout);
			// a tester ! abort if pb
			status=drms_setkey_string(record_out, "PROCSTAT",	PROCSTAT_OK);
			// a tester ! abort if pb

			if (debug) printf("end write in the db %d\n",status);

			if (!status) remove(filenameout); // test if this is ok too!!!
			else 
				{
					sprintf(log_msg,"can't drms_segment_write_from_file() for fsn# %u", fsn);
					lf_logmsg("ERROR", "DRMS", ERR_DRMS_WRITE, status, log_msg, log_msg_code, opf);
					return(ERR_DRMS_WRITE);
				}
				
			drms_free_array (data_array1);
			drms_free_array (data_array2);
			drms_free_array (data_array3);
			drms_free_array (data_array4);
			
			if (debug) printf("done...%u\n",fsn);							
			lf_logmsg("INFO", "APP", 0, 0, "End writing in the DB", log_msg_code, opf);
		}
		else 
		{
			sprintf(log_msg,"NO LDFS file created (ifail >0) or retcode <0 for fsn# %u", fsn);
			lf_logmsg("WARNING", "APP", ERR_LIMBFIT_FAILED, lf_retcode, log_msg, log_msg_code, opf);			
			write_mini_output(fsn,PROCSTAT_NO_LF_FAILED,record_in,record_out,opf,debug); 
			// is this the right retcode? or should I make a test???
			return(0);
		}
	}
	
return 0;
}

int	write_mini_output(unsigned int fsn, char * errcode, DRMS_Record_t *record_in,DRMS_Record_t *record_out,FILE *opf, int debug)
{
		int status;
		char *log_msg_code="write_mini_output";
		lf_logmsg("INFO", "APP", 0, 0, "Writing only header", log_msg_code, opf);			
		// write only the header of segment in the DB
		drms_copykeys(record_out, record_in, 1, kDRMS_KeyClass_All); 	
		// change DATE to our proc_date, do not keep level 1 DATE:
		drms_keyword_setdate(record_out); 
		char code_name[10]=CODE_NAME;
		char code_version[10]=CODE_VERSION;
		char code_date[30]=CODE_DATE;

		status=drms_setkey_string(record_out, "SERIESCN", "Limbfit data");
		status=drms_setkey_string(record_out, "CODENAME",	code_name);
		status=drms_setkey_string(record_out, "CODEVERS",	code_version);
		status=drms_setkey_string(record_out, "CODEDATE",	code_date);
		status=drms_setkey_string(record_out, "PROCSTAT",	errcode);
// need to test status...


return(0); // change to 1 to see what happened and   to -1
}
