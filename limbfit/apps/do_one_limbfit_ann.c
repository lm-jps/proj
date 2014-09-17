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
	img = drms_segment_read(segment_in, DRMS_TYPE_INT, &rstatus);
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
//		input->img_sz0=img->axis[0];
//		input->img_sz1=img->axis[1];
//		input->data=img->data;

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
//int 	w		 = ANNULUS_WIDTH;
long 	S		 = MAX_SIZE_ANN_VARS;
int		naxis_row=img->axis[0];
int		naxis_col=img->axis[1];

/************************************************************************/
/*                        set parameters                               */
/************************************************************************/
long npixels=naxis_row*naxis_col;
long ii, jj, jk;//, i,j; 
 float cmx, cmy,r;//, guess_cx, guess_cy, guess_r;
// int nitr=0, ncut=0;
// 
//r  = (float)naxis_row/2;
long nanval=EQNANVAL;	
/* Initial guess estimate of center position from Richard's code*/
// 

int *data=img->data;

/*
int* out_img = (int *) malloc(sizeof(int)*(naxis_row)*(naxis_col));

//init
	for(ii = 0; ii < naxis_row; ii++)
	{
 			for(jj = 0; jj < naxis_col; jj++) 
			{
					out_img[ii,jj]=data[ii,jj];
			}
			
	}		

*/

/************************************************************************/
/*                        Make annulus data                             */
/************************************************************************/
int nbc=3;
double iimcmy2,jjmcmx;
float d,w2p,w2m;
jk=-1;

//float *p_anls=ios->pf_anls;		
// 
 	if (ios->is_firstobs == 0) // = if fixed size
 	{
		cmx=2048.;
		cmy=2048.;
 		w2p=1985.;
 		w2m=1825.;
		r=(w2p+w2m)/2;
 	}
 	else // = variable size
 	{
		cmx=(float)input->ix;
		cmy=(float)input->iy;
		r=(float)input->ir;
 		w2p=r+40.;
 		w2m=r-40.;
 	}
 	/* Select Points */
	for(ii = 0; ii < naxis_row; ii++)
	{
 			iimcmy2=(ii-cmy)*(ii-cmy);
 			for(jj = 0; jj < naxis_col; jj++) 
 			{ 
 				jjmcmx=jj-cmx;
 				d=(float)sqrt(iimcmy2+(jjmcmx)*(jjmcmx));
 				if (d<w2m || d>w2p)
 				{
 					jk++;
// 					*(p_anls++)=(float)jj;
// 					*(p_anls++)=(float)ii;
// 					*(p_anls++)=data[ii*naxis_col+jj];					 
					//img->data[ii*naxis_col+jj]=nanval;
					data[ii*naxis_col+jj]=nanval;
					//out_img[ii*naxis_col+jj]=nanval;

 				}
 			}
	}

// 		ios->anls_nbpix=jk;
// 		ios->pl_anls=p_anls-1;
// 	
 		//if (results->debug)
 		//{
 			long t_is=(4096*4096)-jk;
 			sprintf(log_msg," total ann points = %ld (total points changed: jk = %ld), r = %f, w2m = %f, w2p = %f, ", t_is,jk,r,w2m,w2p);
 			lf_logmsg("DEBUG", "APP", 0, 0, log_msg, log_msg_code, results->opf);
 		//}
 		if ((jk*3) >= S) 
 		{
 			lf_logmsg("ERROR", "APP", ERR_SIZE_ANN_TOO_BIG, 0,"nbc>S", log_msg_code, results->opf);
 			return ERR_SIZE_ANN_TOO_BIG;
 		}
// 	}
// 	else 
// 	{
// 		jk=ios->anls_nbpix;
// 		ii=0;
// 		p_anls=p_anls+2;
// 		while(p_anls<=ios->pl_anls) 
// 		{
// 			*(p_anls)=data[(int)ios->anls[nbc*ii+1]*naxis_col+(int)ios->anls[nbc*ii]];
// 			ii++;
// 			p_anls=p_anls+3;
// 		}		
// 	}


/* end annulus */

//		lf_retcode=limbfit(input,results,ios);
		
//		drms_free_array(img);

		if (results->debug) 
		{
			sprintf(log_msg,"XC %f", results->cenx);
			lf_logmsg("DEBUG", "INFO", 0, 0, log_msg, log_msg_code, results->opf);
			sprintf(log_msg,"YC %f", results->ceny);
			lf_logmsg("DEBUG", "INFO", 0, 0, log_msg, log_msg_code, results->opf);
			sprintf(log_msg,"R %f", results->radius);
			lf_logmsg("DEBUG", "INFO", 0, 0, log_msg, log_msg_code, results->opf);
			//sprintf(log_msg,"Returned code %d", lf_retcode);
			//lf_logmsg("DEBUG", "INFO", 0, 0, log_msg, log_msg_code, results->opf);
			sprintf(log_msg,"err1 %d", results->error1);
			lf_logmsg("DEBUG", "INFO", 0, 0, log_msg, log_msg_code, results->opf);
			sprintf(log_msg,"err2 %d", results->error2);
			lf_logmsg("DEBUG", "INFO", 0, 0, log_msg, log_msg_code, results->opf);
		}
		
		// >= 0 => all
		// ERR_LIMBFIT_FIT_FAILED => all
		// ERR_LIMBFIT_FAILED => mini

//		if (input->fldf==0) results->fldfr=4;

//		if (lf_retcode >= 0 || lf_retcode == ERR_LIMBFIT_FIT_FAILED)
//		{	
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
			//DRMS_Array_t *array1=drms_array_create(DRMS_TYPE_INT,2,naxes,out_img,&rstatus);			
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
			drms_free_array(img);

			
			rstatus=0;
		

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
