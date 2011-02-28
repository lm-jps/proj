/* 
	I.Scholl 

		Second level wrapper for the limbfit code
		Can be used as a callable module (to be plugged into lev0lev1 module or by the main wrapper)
			calls: limbfit() and reformat output data to be rewritten in the DB
			
	
	#define lfr->code_name 		"limbfit"
	#define lfr->code_version 	"V1.9r0" 
	#define lfr->code_date 		"Mon Feb 28 11:45:21 PST 2011" 
*/

#include "limbfit.h"

void lf_logmsg4fitsio(char *log_msg,char *log_msg_code,char *kw,unsigned int fsn,int status, FILE *opf)
{	
	sprintf(log_msg,"can't fits_update_key[%s] for fsn# %u", kw,fsn);
	lf_logmsg("ERROR", "FITSIO", ERR_FITSIO, status, log_msg, log_msg_code, opf);			
	fits_report_error(opf, status); 
}

int get_set_kw(int typ, char *kw, char *kw_txt, unsigned int fsn, 
		DRMS_Record_t *record_in,DRMS_Record_t *record_out, fitsfile *outfptr, FILE *opf, int debug, int tbf, LIMBFIT_OUTPUT *lfr, int *status)
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
		write_mini_output(PROCSTAT_NO_LF_DB_READ_PB,record_in,record_out,opf,tbf,lfr,debug);
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
		write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,tbf,lfr,debug);
		*status=ERR_DRMS_READ_MISSING_DATA;   
		return ERR_FITSIO;
	}
return(0);
}

int do_one_limbfit(unsigned int fsn, DRMS_Record_t *record_in,DRMS_Record_t *record_out,char *tmp_dir,FILE *opf, int spe, char *dsin, char *comment, int debug, int *status)
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

	lfr->code_name=CODE_NAME;
	lfr->code_version=CODE_VERSION;
	lfr->code_date=CODE_DATE;
	lfr->comment=comment;
	lfr->dsin=dsin;

	seg_cnt = drms_record_numsegments (record_in);
	if (seg_cnt < 1) 
	{
		sprintf(log_msg,"no segments in selected record_in");
		lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_DATA, 0, log_msg, log_msg_code,opf);
		write_mini_output(PROCSTAT_NO_LF_DB_READ_PB,record_in,record_out,opf,0,NULL,debug);
		return(ERR_DRMS_READ_MISSING_DATA);   
	}
	segment_in = drms_segment_lookupnum(record_in, 0);
	img = drms_segment_read(segment_in, DRMS_TYPE_FLOAT, &rstatus);
	if(!img) {
		sprintf(log_msg,"Can't do drms_segment_read() status=%d", rstatus);
		lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_DATA, 0, log_msg, log_msg_code, opf);
		write_mini_output(PROCSTAT_NO_LF_DB_READ_PB,record_in,record_out,opf,0,NULL,debug);
		return(ERR_DRMS_READ_MISSING_DATA);   
	} 
	else
	{
		// add a test if these values are ok, otherwise terminate this observation
		lfv->ix = drms_getkey_float(record_in, "X0_LF", &rstatus);
		if(rstatus) {
			lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_DATA,rstatus, "drms_getkey_float(X0_LF)", log_msg_code, opf);
			write_mini_output(PROCSTAT_NO_LF_DB_READ_PB,record_in,record_out,opf,0,lfr,debug);
			return(ERR_DRMS_READ_MISSING_DATA);   
		}
		if(isnan(lfv->ix)) {
			lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_XYR_LF,rstatus, "X0_LF missing", log_msg_code, opf);
			write_mini_output(PROCSTAT_NO_LF_XYR_LF_MISSING,record_in,record_out,opf,0,lfr,debug);
			return(ERR_DRMS_READ_MISSING_XYR_LF);   
		}
		lfv->iy = drms_getkey_float(record_in, "Y0_LF", &rstatus);
		if(rstatus) {
			lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_DATA,rstatus, "drms_getkey_float(Y0_LF)", log_msg_code, opf);
			write_mini_output(PROCSTAT_NO_LF_DB_READ_PB,record_in,record_out,opf,0,lfr,debug);
			return(ERR_DRMS_READ_MISSING_DATA);   
		}
		if(isnan(lfv->iy)) {
			lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_XYR_LF,rstatus, "Y0_LF missing", log_msg_code, opf);
			write_mini_output(PROCSTAT_NO_LF_XYR_LF_MISSING,record_in,record_out,opf,0,lfr,debug);
			return(ERR_DRMS_READ_MISSING_XYR_LF);   
		}
		lfv->ir = drms_getkey_float(record_in, "RSUN_LF", &rstatus);
		  if(rstatus) {
			lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_DATA,rstatus, "drms_getkey_float(RSUN_LF)", log_msg_code, opf);
			write_mini_output(PROCSTAT_NO_LF_DB_READ_PB,record_in,record_out,opf,0,lfr,debug);
			return(ERR_DRMS_READ_MISSING_DATA);   
		}
		if(isnan(lfv->ir)) {
			lf_logmsg("ERROR", "DRMS", ERR_DRMS_READ_MISSING_XYR_LF,rstatus, "RSUN_LF missing", log_msg_code, opf);
			write_mini_output(PROCSTAT_NO_LF_XYR_LF_MISSING,record_in,record_out,opf,0,lfr,debug);
			return(ERR_DRMS_READ_MISSING_XYR_LF);   
		}
		lfv->img_sz0=img->axis[0];
		lfv->img_sz1=img->axis[1];
		lfv->spe=spe;
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
			
			//---------------------------------------
			//		1) get all keywords
			//---------------------------------------
			drms_copykeys(record_out, record_in, 1, kDRMS_KeyClass_All); 	
	
			// change DATE to our proc_date, do not keep level 1 DATE:
			drms_keyword_setdate(record_out); 
			//---------------------------------------
			//		2) set all keywords
			//---------------------------------------
			int num_ext=3;
// need to test status...
			rstatus=drms_setkey_string(record_out, "SERIESCN", 	series_name);
			rstatus=drms_setkey_string(record_out, "INSERIES", 	lfr->dsin);
			rstatus=drms_setkey_string(record_out, "CODENAME",	lfr->code_name);
			rstatus=drms_setkey_string(record_out, "CODEVERS",	lfr->code_version);
			rstatus=drms_setkey_string(record_out, "CODEDATE",	lfr->code_date);
			rstatus=drms_setkey_int	  (record_out, "FITS_NXT",	num_ext);  
			rstatus=drms_setkey_float (record_out, "X_LFS",		lfr->cenx);
			rstatus=drms_setkey_float (record_out, "Y_LFS",		lfr->ceny);
			rstatus=drms_setkey_double(record_out, "R_LFS",		lfr->radius);
			rstatus=drms_setkey_int   (record_out, "QUAL_LFS",	lfr->quality);
			rstatus=drms_setkey_int   (record_out, "ERRO_LFS",	lfr->error1);
			rstatus=drms_setkey_int   (record_out, "ERRO_FIT",	lfr->error2);
			rstatus=drms_setkey_int   (record_out, "NB_FBINS",	lfr->nb_fbins);
			rstatus=drms_setkey_double(record_out, "MAX_LIMB",	lfr->max_limb);
			rstatus=drms_setkey_double(record_out, "CMEAN",		lfr->cmean);
			rstatus=drms_setkey_int   (record_out, "ANN_WD",	lfr->ann_wd);
			rstatus=drms_setkey_int   (record_out, "MXSZANNV",	lfr->mxszannv);
			rstatus=drms_setkey_int   (record_out, "NB_LDF",	lfr->nb_ldf);
			rstatus=drms_setkey_int   (record_out, "NB_RDB",	lfr->nb_rdb);
			rstatus=drms_setkey_int   (record_out, "NB_ABB",	lfr->nb_abb);
			rstatus=drms_setkey_double(record_out, "UP_LIMIT",	lfr->up_limit);
			rstatus=drms_setkey_double(record_out, "LO_LIMIT",	lfr->lo_limit);
			rstatus=drms_setkey_double(record_out, "INC_X",		lfr->inc_x);
			rstatus=drms_setkey_double(record_out, "INC_Y",		lfr->inc_y);
			rstatus=drms_setkey_int   (record_out, "NFITPNTS",	lfr->nfitpnts);
			rstatus=drms_setkey_int   (record_out, "NB_ITER",	lfr->nb_iter);
			rstatus=drms_setkey_double(record_out, "AHI",		lfr->ahi);
			rstatus=drms_setkey_int   (record_out, "SKIPGC",	lfr->skipgc);
			rstatus=drms_setkey_string(record_out, "COMMENTS",	lfr->comment);

			//---------------------------------------
			//		3) write segments
			//---------------------------------------
/*
			int i_naxes[2];
			i_naxes[0]=lfr->fits_ldfs_naxis1;
			i_naxes[1]=lfr->fits_ldfs_naxis2;
			//ldfs
			DRMS_Array_t *data_array1;
			segment_out = drms_segment_lookupnum (record_out, 0);
			drms_setkey_string(record_out, "FILETYPE[0]", "LDFS");
			data_array1 = drms_array_create (DRMS_TYPE_FLOAT, 2, i_naxes, (void *)lfr->fits_ldfs_data, &rstatus);
			data_array1->bscale = 1.0;
			data_array1->bzero = 0.0;
			rstatus = drms_segment_write (segment_out, data_array1, 0);
			if (rstatus) 
			{
				sprintf(log_msg,"can't create segment 0 (LDF) for fsn# %u in series", fsn);
				lf_logmsg("ERROR", "DRMS", ERR_DRMS_WRITE, rstatus, log_msg, log_msg_code, opf);			
				close_on_error(	record_in,record_out,data_array1);//,opf);
				return ERR_DRMS_WRITE;
			}
*/
			//---------------------------------------
			//		4) write the generic segment 
			//			(as the full FITS file)
			//---------------------------------------
// something left...: in case of problem, filenameout not closed and therefore not removed...

			long nelements, firstrow, firstelem;
			static char filenameout[128];
			sprintf(filenameout,"%s%u_full.fits",tmp_dir,fsn);
		
			firstrow=1;
			firstelem=1;
			*status=0;
			rstatus=0;
			remove(filenameout); 
			fitsfile *outfptr;
			if (fits_create_file(&outfptr, filenameout, &rstatus)) 
			{
				fits_report_error(opf, rstatus); 
				lf_logmsg4fitsio(log_msg, log_msg_code, "fits_create_file()",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,2,lfr,debug);
				return ERR_FITSIO;
			}

			long l_naxes[2];
			l_naxes[0]=lfr->fits_ldfs_naxis1;
			l_naxes[1]=lfr->fits_ldfs_naxis2;
			nelements=l_naxes[0]*l_naxes[1];
			int nax=2;
			if ( fits_create_img(outfptr, FLOAT_IMG, nax, l_naxes, &rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "fits_create_img()",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,2,lfr,debug);
				return ERR_FITSIO;
			}
			if ( fits_write_img(outfptr, TFLOAT, 1, nelements, lfr->fits_ldfs_data, &rstatus) ) 
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "fits_write_img()",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,2,lfr,debug);
				return ERR_FITSIO;
			}	
//			drms_free_array (data_array1); // NOTE: because of that, from now tbf MUST BE EQUAL to 1

		//add more KWs
			if(get_set_kw(1,"ORIGIN","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(1,"TELESCOP","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(1,"INSTRUME","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status); 
			
			TIME tobs;
			static char s_tobs[50];
			tobs = drms_getkey_time(record_in, "DATE__OBS", &rstatus);
			sprint_time(s_tobs, tobs, "UTC", 0);
			if ( fits_update_key(outfptr, TSTRING, "DATE-OBS", &s_tobs, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(DATE__OBS)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}
			tobs = drms_getkey_time(record_in, "T_OBS", &rstatus);
			sprint_time(s_tobs, tobs, "UTC", 0);
			if ( fits_update_key(outfptr, TSTRING, "T_OBS", s_tobs, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(T_OBS)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}
				
			if(get_set_kw(2,"CAMERA","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(1,"IMG_TYPE","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(4,"EXPTIME","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(3,"EXPSDEV","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(2,"WAVELNTH","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(1,"WAVEUNIT","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if ( fits_update_key(outfptr, TINT, "FSN", &fsn, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(FSN)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}
			if(get_set_kw(2,"FID","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(2,"QUALLEV0","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(2,"QUALITY","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(2,"DATAMIN","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(2,"DATAMAX","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(2,"DATAMEDN","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(3,"DATAMEAN","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(3,"DATARMS","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(2,"MISSVALS","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(1,"FLAT_REC","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(1,"CTYPE1","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(1,"CUNIT1","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(3,"CRVAL1","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(3,"CDELT1","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(3,"CRPIX1","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(1,"CTYPE2","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(1,"CUNIT2","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(3,"CRVAL2","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(3,"CDELT2","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(3,"CRPIX2","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(3,"CROTA2","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(3,"R_SUN","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(1,"MPO_REC","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(3,"INST_ROT","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);			
			if(get_set_kw(3,"HPL1POS","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);			
			if(get_set_kw(3,"HPL2POS","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);			
			if(get_set_kw(3,"HPL3POS","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);			
			if (isnan(lfv->ix)) { lfv->ix=0.0;}
			if ( fits_update_key(outfptr, TDOUBLE, "X0_LF", &lfv->ix, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(X0_LF)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}
			if (isnan(lfv->iy)) { lfv->iy=0.0;}
			if ( fits_update_key(outfptr, TDOUBLE, "Y0_LF", &lfv->iy, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(Y0_LF)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}
			if (isnan(lfv->ir)) { lfv->ir=0.0;}
			if ( fits_update_key(outfptr, TDOUBLE, "RSUN_LF", &lfv->ir, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(RSUN_LF)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}
			if(get_set_kw(1,"ASD_REC","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(3,"SAT_ROT","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(1,"ORB_REC","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(4,"DSUN_OBS","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(4,"RSUN_OBS","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(4,"HAEX_OBS","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(4,"HAEY_OBS","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(4,"HAEZ_OBS","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(4,"HPLTID","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(1,"HWLTNSET","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(4,"HCFTID","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if(get_set_kw(4,"HFTSACID","",fsn,record_in,record_out,outfptr,opf,debug,1,lfr,status)) return(*status);
			if ( fits_update_key(outfptr, TSTRING, "SERIESCN", series_name, "Series Name",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(SERIESCN)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TSTRING, "FILETYPE", file_type, "File content",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(FILETYPE)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TSTRING, "CODENAME", lfr->code_name, "Code Name",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(CODENAME)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TSTRING, "CODEVERS", lfr->code_version, "Code Version",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(CODEVERS)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TSTRING, "CODEDATE", lfr->code_date, "Code Date",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(CODEDATE)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TSTRING, "PROCSTAT", proc_stat, "Processing rstatus code",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(PROCSTAT)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}
		
			tobs=drms_keyword_getdate(record_out);		
			sprint_time(s_tobs, tobs, "UTC", 0);
			if ( fits_update_key(outfptr, TSTRING, "DATE", s_tobs, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(DATE)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "FITS_NXT", &num_ext, "Number of extensions",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(FITS_NXT)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TFLOAT, "X_LFS", &lfr->cenx, "Center Position X",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(X_LFS)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TFLOAT, "Y_LFS", &lfr->ceny, "Center Position Y",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(Y_LFS)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TDOUBLE, "R_LFS", &lfr->radius, "Solar Radius",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(R_LFS)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TDOUBLE, "CMEAN", &lfr->cmean, "Mean value of the 500x500 pixels center",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(CMEAN)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TDOUBLE, "MAX_LIMB", &lfr->cmean, "Maximum value on the limb",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(MAX_LIMB)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "QUAL_LFS", &lfr->quality, "Processing quality indicator",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(QUAL_LFS)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "ERRO_LFS", &lfr->error1, "Limb.f Processing error code",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(ERRO_LFS)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "ERRO_FIT", &lfr->error2, "Fitting Processing error code",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(ERRO_FIT)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "NB_FBINS", &lfr->nb_fbins, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(NB_FBINS)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "ANN_WD", &lfr->ann_wd, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(ANN_WD)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "MXSZANNV", &lfr->mxszannv, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(MXSZANNV)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "NB_LDF", &lfr->nb_ldf, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(NB_LDF)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "NB_RDB", &lfr->nb_rdb, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(NB_RDB)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "NB_ABB", &lfr->nb_abb, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(NB_ABB)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TDOUBLE, "UP_LIMIT", &lfr->up_limit, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(UP_LIMIT)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TDOUBLE, "LO_LIMIT", &lfr->lo_limit, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(LO_LIMIT)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TDOUBLE, "INC_X", &lfr->inc_x, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(INC_X)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TDOUBLE, "INC_Y", &lfr->inc_y, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(INC_Y)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "NFITPNTS", &lfr->nfitpnts, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(NFITPNTS)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "NB_ITER", &lfr->nb_iter, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(NB_ITER)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "SKIPGC", &lfr->skipgc, "",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(SKIPGC)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}

			rstatus=0;
			int i;			
		
			// AB ---EXT#1---------------------------------------------------------------------
			char *tunitAB[] = { "","" };
			char *ttypeAB[] = { "Alpha", "Beta"};
			char *tformAB[] = { "1E","1E" };
			char extnameAB[]="LDF info / Alpha&Beta";
			if ( fits_create_tbl( outfptr, BINARY_TBL, lfr->fits_ab_nrows, lfr->fits_ab_tfields, 
													ttypeAB, tformAB, tunitAB, extnameAB, &rstatus) )
			{
				fits_report_error(opf, rstatus); 
				return ERR_FITSIO;
			}
			for (i=1;i<=lfr->fits_ab_tfields;i++)
			{
				if (fits_write_col(outfptr, TFLOAT, i, firstrow, firstelem, lfr->fits_ab_nrows, 
												lfr->fits_alpha_beta+(i-1)*lfr->fits_ab_nrows, &rstatus))
				{
					lf_logmsg4fitsio(log_msg, log_msg_code, "fits_write_col(AB)",fsn,rstatus,opf);			
					write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
					return ERR_FITSIO;
				}
			}

			// All Params ---EXT#2-------------------------------------------------------------------------
			char extnameA[] = "LDF info / As / Es / Radius / IP";
			char *tunitAE[] = { "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "" };
			char *tformAE[] = { "1D", "1D", "1D", "1D", "1D", "1D", "1D", "1D", "1D", "1D", "1D", "1D", "1D", "1D"  };
			char *ttypeA[]  = { "A0", "A1", "A2", "A3", "A4", "A5", "E0", "E1", "E2", "E3", "E4", "E5", "Radius", "IP1"};
		
			if ( fits_create_tbl(outfptr, BINARY_TBL, lfr->fits_params_nrows, lfr->fits_params_tfields, 
												ttypeA, tformAE,tunitAE, extnameA, &rstatus) )
			{
					fits_report_error(opf, rstatus); 
					return ERR_FITSIO;
			}

			for (i=1;i<=lfr->fits_params_tfields;i++)
			{
				if (fits_write_col(outfptr, TDOUBLE, i, firstrow, firstelem, lfr->fits_params_nrows, 
												lfr->fits_params+(i-1)*lfr->fits_params_nrows, &rstatus))
				{
					lf_logmsg4fitsio(log_msg, log_msg_code, "fits_write_col(PARAMS)",fsn,rstatus,opf);			
					write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
					return ERR_FITSIO;
				}
			}	
/*
			// Full LDF ---EXT#3-------------------------------------------------------------------------
			char *tunitF[] = { "","" };
			char *ttypeF[] = { "Intensity", "Radius"};
			char *tformF[] = { "1E","1E" };
			char extnameF[]="Full LDF ";
			if ( fits_create_tbl( outfptr, BINARY_TBL, lfr->fits_fldfs_nrows, lfr->fits_fldfs_tfields, 
													ttypeF, tformF, tunitF, extnameF, &rstatus) )
			{
				fits_report_error(opf, rstatus); 
				return ERR_FITSIO;
			}
			for (i=1;i<=lfr->fits_fldfs_tfields;i++)
			{
				if (fits_write_col(outfptr, TFLOAT, i, firstrow, firstelem, lfr->fits_fldfs_nrows, 
												lfr->fits_fulldfs+(i-1)*lfr->fits_fldfs_nrows, &rstatus))
				{
					lf_logmsg4fitsio(log_msg, log_msg_code, "fits_write_col(FLDF)",fsn,rstatus,opf);			
					write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
					return ERR_FITSIO;
				}
			}
*/			
			//---------------------------------------
			//		5) close & free
			//---------------------------------------
			if ( fits_update_key(outfptr, TSTRING, "PROCSTAT", proc_stat, "Processing rstatus code",&rstatus) )
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "(PROCSTAT)",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}

			// Close file & attach to DB & remove
			if ( fits_close_file(outfptr, &rstatus) ) 
			{
				lf_logmsg4fitsio(log_msg, log_msg_code, "fits_close_file()",fsn,rstatus,opf);			
				write_mini_output(PROCSTAT_NO_LF_FITS_WRITE_PB,record_in,record_out,opf,1,lfr,debug);
				return ERR_FITSIO;
			}						
			//warning change segment number !!!
			segment_out = drms_segment_lookupnum (record_out, 1); //WARNING change this number if the number of segment changes...
			rstatus=drms_segment_write_from_file(segment_out,filenameout);
			// a tester ! abort if pb
			rstatus=drms_setkey_string(record_out, "PROCSTAT",	PROCSTAT_OK);
			// a tester ! abort if pb

			if (debug) printf("end write in the db %d\n",rstatus);

			if (!rstatus) remove(filenameout); // test if this is ok too!!!
			else 
				{
					sprintf(log_msg,"can't drms_segment_write_from_file() for fsn# %u", fsn);
					lf_logmsg("ERROR", "DRMS", ERR_DRMS_WRITE, rstatus, log_msg, log_msg_code, opf);
					return(ERR_DRMS_WRITE);
				}
				
			free(lfr->fits_alpha_beta);
			free(lfr->fits_params);
			free(lfr->fits_fulldfs);
			
			if (debug) printf("done...%u\n",fsn);							
			lf_logmsg("INFO", "APP", 0, 0, "End writing in the DB", log_msg_code, opf);
		}
		else 
		{
			sprintf(log_msg,"NO LDFS file created (ifail >0) or retcode <0 for fsn# %u", fsn);
			lf_logmsg("WARNING", "APP", ERR_LIMBFIT_FAILED, lf_retcode, log_msg, log_msg_code, opf);			
			write_mini_output(PROCSTAT_NO_LF_FAILED,record_in,record_out,opf,0,lfr,debug);
			// is this the right retcode? or should I make a test???
			
			// free 
			if (lfr->error1 == 0)
			{
				free(lfr->fits_ldfs_data);
				free(lfr->fits_alpha_beta);
				free(lfr->fits_params);
				free(lfr->fits_fulldfs);
			}
			return(0);
		}
	}

return 0;
}
//pb here: quand this is called then no free are executed!!! 
int	write_mini_output(char * errcode, DRMS_Record_t *record_in,DRMS_Record_t *record_out,FILE *opf, int tbf, LIMBFIT_OUTPUT *lfr, int debug)
{
		//int status;
		char *log_msg_code="write_mini_output";
		lf_logmsg("INFO", "APP", 0, 0, "Writing only header", log_msg_code, opf);			
		// write only the header of segment in the DB
		drms_copykeys(record_out, record_in, 1, kDRMS_KeyClass_All); 	
		// change DATE to our proc_date, do not keep level 1 DATE:
		drms_keyword_setdate(record_out); 

// need to test status...
		drms_setkey_string(record_out, "SERIESCN", "Limbfit data");
		drms_setkey_string(record_out, "CODENAME",	lfr->code_name);
		drms_setkey_string(record_out, "CODEVERS",	lfr->code_version);
		drms_setkey_string(record_out, "CODEDATE",	lfr->code_date);
		drms_setkey_string(record_out, "INSERIES",	lfr->dsin);
		drms_setkey_string(record_out, "COMMENTS",	lfr->comment);
		drms_setkey_string(record_out, "PROCSTAT",	errcode);
		
		if (tbf >=1)
		{
			free(lfr->fits_alpha_beta);
			free(lfr->fits_params);
			free(lfr->fits_fulldfs);
			if (tbf ==2) free(lfr->fits_ldfs_data);	
		}
return(0);
}
