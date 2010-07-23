/* 
	I.Scholl 

		Second level wrapper for the limbfit code
		Can be used as a callable module (to be plugged into lev0lev1 module or by the main wrapper)
			calls: limbfit() and reformat output data to be rewritten in the DB
			
	
	#define CODE_NAME 		"limbfit"
	#define CODE_VERSION 	"V1.2r0" 
	#define CODE_DATE 		"Wed Jul 21 09:32:15 PDT 2010" 
*/

#include "fitsio.h"
#include "limbfit.h"

int do_one_limbfit(LIMBFIT_INPUT *info,LIMBFIT_OUTPUT *results,DRMS_Record_t *record_in,DRMS_Record_t *record_out,int debug)
{

	int retcode;
//	char s_tobs[64];

	DRMS_Segment_t *segment_out;

//?	sprint_time (s_tobs, info->tobs, "UT", 3);
//?	fprintf(info->opf,"T_OBS %s \n",s_tobs);
	fprintf(info->opf,"CAMERA %d \n",info->camera);
	fprintf(info->opf,"FID %d \n",info->fid);
	fprintf(info->opf,"info->fsn %u \n",info->fsn);
	fprintf(info->opf,"HPLTID %d \n",info->hpltid);
	fprintf(info->opf,"X0_LF %f \n",info->ix);
	fprintf(info->opf,"Y0_LF %f \n",info->iy);
	fprintf(info->opf,"RSUN_LF %f \n",info->ir);
			
	/***********************************************************************
		Call the limbfit routine like build_lev1.c would do it
		- update the structure with parameters got from the command line
	***********************************************************************/
	
	retcode=limbfit(info,results,debug);
	
	fprintf(info->opf,"returned code %d\n",retcode);
	fprintf(info->opf,"after limbfit x %6.3f\n",results->cenx);
	fprintf(info->opf,"after limbfit y %6.3f\n",results->ceny);
	fprintf(info->opf,"after limbfit r %6.3f\n",results->radius);
	fprintf(info->opf,"after limbfit err %d\n",results->error);
	fprintf(info->opf,"after limbfit cmean %6.3f\n",results->cmean);
	
	if (retcode >= 0)
	{	
		/**********************************************************************
			Write into the DB
		**********************************************************************/
		int status;
		
		//---------------------------------------
		//		1) get all keywords
		//---------------------------------------


		//---------------------------------------
		//		2) set all keywords
		//---------------------------------------
		status=drms_setkey_string(record_out, "PROCDATE",	results->proc_date);
		status=drms_setkey_string(record_out, "SERIESCN", "Limbfit data");
		status=drms_setkey_string(record_out, "CODENAME",	results->code_name);
		status=drms_setkey_string(record_out, "CODEVERS",	results->code_version);
		status=drms_setkey_float (record_out, "X_LFS",		results->cenx);
		status=drms_setkey_float (record_out, "Y_LFS",		results->ceny);
		status=drms_setkey_double(record_out, "R_LFS",		results->radius);
		status=drms_setkey_int   (record_out, "QUAL_LFS",	results->quality);
		status=drms_setkey_int   (record_out, "ERRO_LFS",	results->error);
		status=drms_setkey_double(record_out, "MAX_LIMB",	results->max_limb);
		status=drms_setkey_int   (record_out, "ANN_WD",		results->ann_wd);
		status=drms_setkey_int   (record_out, "MXSZANNV",	results->mxszannv);
		status=drms_setkey_int   (record_out, "NB_LDF",		results->nb_ldf);
		status=drms_setkey_int   (record_out, "NB_RDB",		results->nb_rdb);
		status=drms_setkey_int   (record_out, "NB_ABB",		results->nb_abb);
		status=drms_setkey_double(record_out, "UP_LIMIT",	results->up_limit);
		status=drms_setkey_double(record_out, "LO_LIMIT",	results->lo_limit);
		status=drms_setkey_double(record_out, "INC_X",		results->inc_x);
		status=drms_setkey_double(record_out, "INC_Y",		results->inc_y);
		status=drms_setkey_int   (record_out, "NFITPNTS",	results->nfitpnts);
		status=drms_setkey_int   (record_out, "NB_ITER",	results->nb_iter);
		status=drms_setkey_int   (record_out, "SKIPGC",		results->skipgc);
		// + set those just got
		drms_copykeys(record_out, record_in, 1, kDRMS_KeyClass_All); 	
		//---------------------------------------
		//		3) write segments
		//---------------------------------------
		int i_naxes[2];
		i_naxes[0]=results->fits_ldfs_naxis1;
		i_naxes[1]=results->fits_ldfs_naxis2;
		//ldfs
		segment_out = drms_segment_lookupnum (record_out, 0);
		drms_setkey_string(record_out, "FILETYPE[0]", "LDFS");
		data_array = drms_array_create (DRMS_TYPE_FLOAT, 2, i_naxes, (void *)results->fits_ldfs_data, &status);
		data_array->bscale = 1.0;
		data_array->bzero = 0.0;
		status = drms_segment_write (segment_out, data_array, 0);
		if (status) 
		{
			fprintf(info->opf,"can't create segment 0 (LDF) for info->fsn# %u in series %s\n",info->fsn,info->dsout);
			drms_free_array (data_array);
			drms_close_record (record_out, DRMS_FREE_RECORD);
			drms_close_record (record_in, DRMS_FREE_RECORD);
			fclose(info->opf);
			return ERR_DRMS_WR;
		}
		drms_free_array (data_array);

		//full LDF
		segment_out = drms_segment_lookupnum (record_out, 1);
		i_naxes[0]=results->fits_fldfs_tfields;
		i_naxes[1]=results->fits_fldfs_nrows; 					// to check row/col
		drms_setkey_string(record_out, "FILETYPE[1]", "FULL_LDFS");
		data_array = drms_array_create (DRMS_TYPE_FLOAT, 2, i_naxes, (void *)results->fits_fulldfs, &status) ;
		data_array->bscale = 1.0;
		data_array->bzero = 0.0;
		status = drms_segment_write (segment_out, data_array, 1);
		if (status) 
		{
			fprintf(info->opf,"can't create segment 0 (fldf) for info->fsn# %u in series %s\n",info->fsn,info->dsout);
			drms_free_array (data_array);
			drms_close_record (record_out, DRMS_FREE_RECORD);
			drms_close_record (record_in, DRMS_FREE_RECORD);
			return ERR_DRMS_WR;
		}
		drms_free_array (data_array);
		
		//alpha_betas
		segment_out = drms_segment_lookupnum (record_out, 2);
		i_naxes[0]=results->fits_ab_tfields;
		i_naxes[1]=results->fits_ab_nrows;
		drms_setkey_string(record_out, "FILETYPE[2]", "AB");
		data_array = drms_array_create (DRMS_TYPE_FLOAT, 2 ,i_naxes, (void *)results->fits_alpha_beta1, &status);
		data_array->bscale = 1.0;
		data_array->bzero = 0.0;
		status = drms_segment_write (segment_out, data_array, 2);
		if (status) 
		{
			fprintf(info->opf,"can't create segment 0 (AB) for info->fsn# %u in series %s\n",info->fsn,info->dsout);
			drms_free_array (data_array);
			drms_close_record (record_out, DRMS_FREE_RECORD);
			drms_close_record (record_in, DRMS_FREE_RECORD);
			fclose(info->opf);
			return ERR_DRMS_WR;
		}
		drms_free_array (data_array);

		//PARAMS
		segment_out = drms_segment_lookupnum (record_out, 3);
		i_naxes[0]=results->fits_params_tfields;
		i_naxes[1]=results->fits_params_nrows;
		drms_setkey_string(record_out, "FILETYPE[3]", "PARAMS");
		data_array = drms_array_create (DRMS_TYPE_DOUBLE, 2, i_naxes, (void *)results->fits_params1, &status) ;
		data_array->bscale = 1.0;
		data_array->bzero = 0.0;
		status = drms_segment_write (segment_out, data_array, 3);
		if (status) 
		{
			fprintf(info->opf,"can't create segment 0 (PARAMS) for info->fsn# %u in series %s\n",info->fsn,info->dsout);
			drms_free_array (data_array);
			drms_close_record (record_out, DRMS_FREE_RECORD);
			drms_close_record (record_in, DRMS_FREE_RECORD);
			return ERR_DRMS_WR;
		}
		drms_free_array (data_array);

		//---------------------------------------
		//		4) write the generic segment 
		//			(as the full FITS file)
		//---------------------------------------

			long nelements, firstrow, firstelem;
			static char filenameout[128];
			//sprintf(filenameout,"%s%u_full.fits",TMP_DIR,info->fsn);
			sprintf(filenameout,"./%u_full.fits",info->fsn);
		
			firstrow=1;
			firstelem=1;
			status=0;
		
			remove(filenameout); 
			fitsfile *outfptr;
			if (fits_create_file(&outfptr, filenameout, &status)) 
			{
				fits_report_error(stdout, status); 
				fclose(info->opf);
				//return ERR_FITSIO;
			}

			long l_naxes[2];
   			l_naxes[0]=results->fits_ldfs_naxis1;
			l_naxes[1]=results->fits_ldfs_naxis2;
			nelements=l_naxes[0]*l_naxes[1];
			int nax=2;
			if ( fits_create_img(outfptr, FLOAT_IMG, nax, l_naxes, &status) )
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}
		
			if ( fits_update_key(outfptr, TSTRING, "PROCDATE", &results->proc_date, "Processing Date",&status) )
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TSTRING, "CODENAME", &results->code_name, "Code Name",&status) )
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TSTRING, "CODEVERS", &results->code_version, "Code Version",&status) )
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TFLOAT, "CENTER_X", &results->cenx, "Center Position X",&status) )
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TFLOAT, "CENTER_Y", &results->ceny, "Center Position Y",&status) )
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TDOUBLE, "RADIUS", &results->radius, "Solar Radius",&status) )
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TDOUBLE, "CMEAN", &results->cmean, "Mean value of the center 500x500 pixels",&status) )
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "QUALITY", &results->quality, "Processing quality indicator",&status) )
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "ERROR", &results->error, "Processing error code",&status) )
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "NUM_EXT", &results->numext, "Number of extensions",&status) )
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}
/*			if ( fits_update_key(outfptr, TSTRING, "T_OBS", &s_tobs, "Observing date (from original data file)",&status) )
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}
*/				
			if (isnan(info->ix)) { info->ix=0.0;}
			if ( fits_update_key(outfptr, TDOUBLE, "X0_LF", &info->ix, "Center Position X from lev1",&status) )
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}
			if (isnan(info->iy)) { info->iy=0.0;}
			if ( fits_update_key(outfptr, TDOUBLE, "Y0_LF", &info->iy, "Center Position Y from lev1",&status) )
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}
			if (isnan(info->ir)) { info->ir=0.0;}
			if ( fits_update_key(outfptr, TDOUBLE, "RSUN_LF", &info->ir, "Solar Radius from lev1",&status) )
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "CAMERA", &info->camera, "Camera ID",&status) )
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "FID", &info->fid, "FID from lev1",&status) )
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "info->fsn", &info->fsn, "info->fsn from lev1",&status) )
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "HPLTID", &info->hpltid, "HPLTID from lev1",&status) )
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}
/*			if ( fits_update_key(outfptr, TSTRING, "HWLTNSET", &s_hwltnset, "HWLTNSET from lev1",&status) )
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}
*/			if ( fits_update_key(outfptr, TFLOAT, "SAT_ROT", &info->sat_rot, "SAT_ROT from lev1",&status) )
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}
			/*
			if ( fits_update_key(outfptr, TFLOAT, "INST_ROT", &inst_rot, "INST_ROT from lev1",&status) )
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}
			*/
/*			if ( fits_update_key(outfptr, TINT, "DATAMAX", &maxval, "Maximum value of all pixels ",&status) )
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}
*/
			if ( fits_write_img(outfptr, TFLOAT, 1, nelements, results->fits_ldfs_data, &status) ) 
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}	

			status=0;
			int i;			
		
			// AB ---EXT#1---------------------------------------------------------------------
			char *tunitAB[] = { "","" };
			char *ttypeAB[] = { "Alpha", "Beta"};
			char *tformAB[] = { "1E","1E" };
			char extnameAB[]="LDF info / Alpha&Beta";
			if ( fits_create_tbl( outfptr, BINARY_TBL, results->fits_ab_nrows, results->fits_ab_tfields, 
													ttypeAB, tformAB, tunitAB, extnameAB, &status) )
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}
			for (i=1;i<=results->fits_ab_tfields;i++)
			{
				if (fits_write_col(outfptr, TFLOAT, i, firstrow, firstelem, results->fits_ab_nrows, 
												results->fits_alpha_beta2+(i-1)*results->fits_ab_nrows, &status))
				{
					fits_report_error(stdout, status); 
					//return ERR_FITSIO;
				}
			}

			// All Params ---EXT#2-------------------------------------------------------------------------
			char extnameA[] = "LDF info / As / Es / Radius / IP";
			char *tunitAE[] = { "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "" };
			char *tformAE[] = { "1D", "1D", "1D", "1D", "1D", "1D", "1D", "1D", "1D", "1D", "1D", "1D", "1D", "1D"  };
			char *ttypeA[]  = { "A0", "A1", "A2", "A3", "A4", "A5", "E0", "E1", "E2", "E3", "E4", "E5", "Radius", "IP1"};
		
			if ( fits_create_tbl(outfptr, BINARY_TBL, results->fits_params_nrows, results->fits_params_tfields, 
												ttypeA, tformAE,tunitAE, extnameA, &status) )
			{
					fits_report_error(stdout, status); 
					//return ERR_FITSIO;
			}

			for (i=1;i<=results->fits_params_tfields;i++)
			{
				if (fits_write_col(outfptr, TDOUBLE, i, firstrow, firstelem, results->fits_params_nrows, 
												results->fits_params2+(i-1)*results->fits_params_nrows, &status))
				{
					fits_report_error(stdout, status); 
					//return ERR_FITSIO;
				}
			}
	

			// Full LDF ---EXT#3-------------------------------------------------------------------------
			// don't know how to do yet because of the variable lengths of columns...
			
			
			// Close file & attach to DB & remove
			if ( fits_close_file(outfptr, &status) ) 
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}						
			//warning change segment number !!!
			segment_out = drms_segment_lookupnum (record_out, 4);
			status=drms_segment_write_from_file(segment_out,filenameout);
			if (!status) remove(filenameout); 
				else fprintf(info->opf,"drms_segment_write_from_file failed for fsn %u\n",info->fsn);
		//---------------------------------------
		//		5) close & free
		//---------------------------------------
		//free(data); 
		//free(filenameout);
		
		//+ free all variables coming from limbfit.c		
//printf("ici %p\n",results->fits_ldfs_data);
//		if (results->fits_ldfs_data) free(results->fits_ldfs_data);
/*printf("ici %p\n",results->fits_params1);
		if (results->fits_params1) free(results->fits_params1);
printf("ici %p\n",results->fits_params2);
		if (results->fits_params2) free(results->fits_params2);
printf("ici %p\n",results->fits_alpha_beta1);
		if (results->fits_alpha_beta1) free(results->fits_alpha_beta1);
printf("ici %p\n",results->fits_alpha_beta2);
		if (results->fits_alpha_beta2) free(results->fits_alpha_beta2);
printf("ici %p\n",results->fits_fulldfs);
		if (results->fits_fulldfs) free(results->fits_fulldfs);
*/				
		printf("done...%u\n",info->fsn);							
		fprintf(info->opf,"done...%u\n",info->fsn);							
	}
	else 
	{
		fprintf(info->opf,"NO LDFS file created (ifail >0) or retcode <0\n");
	}



return 0;
}
