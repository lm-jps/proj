#include "fitsio.h"
#include "limbfit.h"

int do_one_limbfit(LIMBFIT_INPUT *info,DRMS_Record_t *record_in,int debug)
{

	static LIMBFIT_OUTPUT limbfit_res ;
	static LIMBFIT_OUTPUT *lfr = &limbfit_res;
	int retcode;
//	char s_tobs[64];

	DRMS_Record_t *record_out;
	DRMS_Segment_t *segment_out;

//?	sprint_time (s_tobs, info->tobs, "UT", 3);
//?	fprintf(info->opf,"T_OBS %s \n",s_tobs);
	fprintf(info->opf,"CAMERA %d \n",info->camera);
	fprintf(info->opf,"FID %d \n",info->fid);
	fprintf(info->opf,"FSN %u \n",info->fsn);
	fprintf(info->opf,"HPLTID %d \n",info->hpltid);
	fprintf(info->opf,"X0_LF %f \n",info->ix);
	fprintf(info->opf,"Y0_LF %f \n",info->iy);
	fprintf(info->opf,"RSUN_LF %f \n",info->ir);
			
	/***********************************************************************
		Call the limbfit routine like build_lev1.c would do it
		- update the structure with parameters got from the command line
	***********************************************************************/
	
	retcode=limbfit(info,lfr,debug);
	
	fprintf(info->opf,"returned code %d\n",retcode);
	fprintf(info->opf,"after limbfit x %6.3f\n",lfr->cenx);
	fprintf(info->opf,"after limbfit y %6.3f\n",lfr->ceny);
	fprintf(info->opf,"after limbfit r %6.3f\n",lfr->radius);
	fprintf(info->opf,"after limbfit err %d\n",lfr->error);
	fprintf(info->opf,"after limbfit cmean %6.3f\n",lfr->cmean);
	
	if (retcode >= 0)
	{	
		/**********************************************************************
			Write the drms record_in
		**********************************************************************/
		DRMS_Array_t *data_array;
		int status;

		record_out=drms_create_record(drms_env, info->dsout, DRMS_PERMANENT, &status);
		if (!record_out) 
		{
			fprintf(info->opf,"can't create record_in for fsn# %u in series %s\n",info->fsn,info->dsout);
			fclose(info->opf);
			return ERR_DRMS_WR;
		}	
		
		//---------------------------------------
		//		1) get all keywords
		//---------------------------------------


		//---------------------------------------
		//		2) set all keywords
		//---------------------------------------
		status=drms_setkey_string(record_out, "PROCDATE",	lfr->proc_date);
		status=drms_setkey_string(record_out, "SERIESCN", "Limbfit data");
		status=drms_setkey_string(record_out, "CODENAME",	lfr->code_name);
		status=drms_setkey_string(record_out, "CODEVERS",	lfr->code_version);
		status=drms_setkey_float (record_out, "X_LFS",		lfr->cenx);
		status=drms_setkey_float (record_out, "Y_LFS",		lfr->ceny);
		status=drms_setkey_double(record_out, "R_LFS",		lfr->radius);
		status=drms_setkey_int   (record_out, "QUAL_LFS",	lfr->quality);
		status=drms_setkey_int   (record_out, "ERRO_LFS",	lfr->error);
		status=drms_setkey_double(record_out, "MAX_LIMB",	lfr->max_limb);
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
		// + set those just got
		drms_copykeys(record_out, record_in, 1, kDRMS_KeyClass_All); 	
		//---------------------------------------
		//		3) write segments
		//---------------------------------------
		int i_naxes[2];
		i_naxes[0]=lfr->fits_ldfs_naxis1;
		i_naxes[1]=lfr->fits_ldfs_naxis2;
		//ldfs
		segment_out = drms_segment_lookupnum (record_out, 0);
//		drms_setkey_string(record_out, "FILETYPE[0]", "LDFS");
		data_array = drms_array_create (DRMS_TYPE_FLOAT, 2, i_naxes, (void *)lfr->fits_ldfs_data, &status);
		data_array->bscale = 1.0;
		data_array->bzero = 0.0;
		status = drms_segment_write (segment_out, data_array, 0);
		if (status) 
		{
			fprintf(info->opf,"can't create segment 0 (LDF) for fsn# %u in series %s\n",info->fsn,info->dsout);
			drms_free_array (data_array);
			drms_close_record (record_out, DRMS_FREE_RECORD);
			drms_close_record (record_in, DRMS_FREE_RECORD);
			fclose(info->opf);
			return ERR_DRMS_WR;
		}
		drms_free_array (data_array);
/*
		//full LDF
		segment_out = drms_segment_lookupnum (record_out, 1);
		naxes[0]=lfr->fits_fldfs_naxis1;
		naxes[1]=lfr->fits_fldfs_naxis2;
		drms_setkey_string(record_out, "FILETYPE[1]", "FULL_LDFS");
		data_array = drms_array_create (DRMS_TYPE_FLOAT, 2, naxes, (void *)lfr->fits_fulldfs, &status) ;
		data_array->bscale = 1.0;
		data_array->bzero = 0.0;
		status = drms_segment_write (segment_out, data_array, 1);
		if (status) 
		{
			fprintf(info->opf,"can't create segment 0 (fldf) for fsn# %u in series %s\n",fsn,dsout);
			drms_free_array (data_array);
			drms_close_record (record_out, DRMS_FREE_RECORD);
			drms_close_record (record_in, DRMS_FREE_RECORD);
			return ERR_DRMS_WR;
		}
		drms_free_array (data_array);
		
		//alpha_betas
		segment_out = drms_segment_lookupnum (record_out, 2);
		naxes[0]=lfr->fits_alpha_beta_naxis1;
		naxes[1]=lfr->fits_alpha_beta_naxis2;
		drms_setkey_string(record_out, "FILETYPE[2]", "AB");
		data_array = drms_array_create (DRMS_TYPE_FLOAT, 2 ,naxes, (void *)lfr->fits_alpha_beta, &status);
		data_array->bscale = 1.0;
		data_array->bzero = 0.0;
		status = drms_segment_write (segment_out, data_array, 2);
		if (status) 
		{
			fprintf(info->opf,"can't create segment 0 (AB) for fsn# %u in series %s\n",fsn,dsout);
			drms_free_array (data_array);
			drms_close_record (record_out, DRMS_FREE_RECORD);
			drms_close_record (record_in, DRMS_FREE_RECORD);
			fclose(info->opf);
			return ERR_DRMS_WR;
		}
		drms_free_array (data_array);

		//PARAMS
		segment_out = drms_segment_lookupnum (record_out, 3);
		naxes[0]=lfr->fits_params_naxis1;
		naxes[1]=lfr->fits_params_naxis2;
		drms_setkey_string(record_out, "FILETYPE[3]", "PARAMS");
		data_array = drms_array_create (DRMS_TYPE_DOUBLE, 2, naxes, (void *)lfr->fits_params, &status) ;
		data_array->bscale = 1.0;
		data_array->bzero = 0.0;
		status = drms_segment_write (segment_out, data_array, 3);
		if (status) 
		{
			fprintf(info->opf,"can't create segment 0 (PARAMS) for fsn# %u in series %s\n",fsn,dsout);
			drms_free_array (data_array);
			drms_close_record (record_out, DRMS_FREE_RECORD);
			drms_close_record (record_in, DRMS_FREE_RECORD);
			return ERR_DRMS_WR;
		}
		drms_free_array (data_array);
*/
		//---------------------------------------
		//		4) make the full fits file
		//---------------------------------------
printf("create fits\n");

			long nelements, firstrow, firstelem;
			static char filenameout[128];
			//sprintf(filenameout,"%sfcur.fits",TMP_DIR);
			sprintf(filenameout,"./fcur.fits");
		
			firstrow=1;
			firstelem=1;
			status=0;
		
			fprintf(info->opf,"F1: -%s-\n",filenameout);	
printf("F1: -%s-\n",filenameout);	
		
			remove(filenameout); 
			fitsfile *outfptr;
			if (fits_create_file(&outfptr, filenameout, &status)) 
			{
				fits_report_error(stdout, status); 
				fclose(info->opf);
				//return ERR_FITSIO;
			}

			long l_naxes[2];
   			l_naxes[0]=lfr->fits_ldfs_naxis1;
			l_naxes[1]=lfr->fits_ldfs_naxis2;
			nelements=l_naxes[0]*l_naxes[1];
			int nax=2;
			if ( fits_create_img(outfptr, FLOAT_IMG, nax, l_naxes, &status) )
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}
		
			if ( fits_update_key(outfptr, TSTRING, "PROCDATE", &lfr->proc_date, "Processing Date",&status) )
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TSTRING, "CODENAME", &lfr->code_name, "Code Name",&status) )
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TSTRING, "CODEVERS", &lfr->code_version, "Code Version",&status) )
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TFLOAT, "CENTER_X", &lfr->cenx, "Center Position X",&status) )
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TFLOAT, "CENTER_Y", &lfr->ceny, "Center Position Y",&status) )
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TDOUBLE, "RADIUS", &lfr->radius, "Solar Radius",&status) )
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TDOUBLE, "CMEAN", &lfr->cmean, "Mean value of the center 500x500 pixels",&status) )
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "QUALITY", &lfr->quality, "Processing quality indicator",&status) )
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "ERROR", &lfr->error, "Processing error code",&status) )
			{
				fits_report_error(stdout, status); 
				//return ERR_FITSIO;
			}
			if ( fits_update_key(outfptr, TINT, "NUM_EXT", &lfr->numext, "Number of extensions",&status) )
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
			if ( fits_update_key(outfptr, TINT, "FSN", &info->fsn, "FSN from lev1",&status) )
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
*/				if ( fits_write_img(outfptr, TFLOAT, 1, nelements, lfr->fits_ldfs_data, &status) ) 
				{
					fits_report_error(stdout, status); 
					//return ERR_FITSIO;
				}	
printf(" after write img %d\n",status);			
				status=0;
				int i;			
				// As ---EXT#1-------------------------------------------------------------------------
				char extnameA[] = "LDF info / As";
				char *tunitAE[] = { "", "", "", "", "", "", "" };
				char *tformAE[] = { "1D", "1D", "1D", "1D", "1D", "1D" };
				char *ttypeA[]  = { "A0", "A1", "A2", "A3", "A4", "A5"};
			
				if ( fits_create_tbl(outfptr, BINARY_TBL, lfr->fits_as_nrows, lfr->fits_as_tfields, 
													ttypeA, tformAE,tunitAE, extnameA, &status) )
				{
						fits_report_error(stdout, status); 
						//return ERR_FITSIO;
				}
printf("status 1a %d\n",status);
				for (i=1;i<=lfr->fits_as_tfields;i++)
				{
					if (fits_write_col(outfptr, TDOUBLE, i, firstrow, firstelem, lfr->fits_as_nrows, 
													lfr->fits_as+(i-1)*lfr->fits_as_nrows, &status))
					{
						fits_report_error(stdout, status); 
						//return ERR_FITSIO;
					}
				}
printf("status 1b %d\n",status);
		
				// Es ----EXT#2------------------------------------------------------------------------
				char *ttypeE[] = { "E0", "E1", "E2", "E3", "E4", "E5"};
				char extnameE[] = "LDF info / Es";
			
				if ( fits_create_tbl(outfptr, BINARY_TBL, lfr->fits_es_nrows, lfr->fits_es_tfields, 
														ttypeE, tformAE, tunitAE, extnameE, &status) )
				   {
						fits_report_error(stdout, status); 
						//return ERR_FITSIO;
					}
				for (i=1;i<=lfr->fits_es_tfields;i++)
				{
					if (fits_write_col(outfptr, TDOUBLE, i, firstrow, firstelem, lfr->fits_es_nrows, 
													(lfr->fits_es+(i-1)*lfr->fits_es_nrows), &status))
					{
						fits_report_error(stdout, status); 
						//return ERR_FITSIO;
					}
				}
printf("status 2 %d\n",status);
			
				// Radius ---EXT#3---------------------------------------------------------------------
				char *tunitR[] = { "" };
				char *ttypeR[] = { "Radius"};
				char *tformR[] = { "1D" };
				char extnameR[] = "LDF info / Radius";
			
				if ( fits_create_tbl(outfptr, BINARY_TBL, lfr->fits_r_nrows, lfr->fits_r_tfields, 
															ttypeR, tformR, tunitR, extnameR, &status) )
				{
					fits_report_error(stdout, status); 
					//return ERR_FITSIO;
				}
				if (fits_write_col(outfptr, TDOUBLE, 1, firstrow, firstelem, lfr->fits_r_nrows, 
																				lfr->fits_r, &status))
				{
					fits_report_error(stdout, status); 
					//return ERR_FITSIO;
				}
printf("status 3 %d\n",status);
			
				// AB ---EXT#4---------------------------------------------------------------------
				char *tunitAB[] = { "","" };
				char *ttypeAB[] = { "Alpha", "Beta"};
				char *tformAB[] = { "1E","1E" };
				char extnameAB[]="LDF info / Alpha&Beta";
				if ( fits_create_tbl( outfptr, BINARY_TBL, lfr->fits_ab_nrows, lfr->fits_ab_tfields, 
														ttypeAB, tformAB, tunitAB, extnameAB, &status) )
				{
					fits_report_error(stdout, status); 
					//return ERR_FITSIO;
				}
				if (fits_write_col(outfptr, TFLOAT, 1, firstrow, firstelem, lfr->fits_ab_nrows, 
																		lfr->fits_ab_dataa, &status))
				{
					fits_report_error(stdout, status); 
					//return ERR_FITSIO;
				}
				if (fits_write_col(outfptr, TFLOAT, 2, firstrow, firstelem, lfr->fits_ab_nrows,
																		lfr->fits_ab_datab, &status))
				{
					fits_report_error(stdout, status); 
					//return ERR_FITSIO;
				}
printf("status 4 %d\n",status);

				// IP ---EXT#5---------------------------------------------------------------------
				char *tunitIP[] = { "" };
				char *ttypeIP[] = { "IP1"};
				char *tformIP[] = { "1D" };
				char extnameIP[]="LDF info / Inflection Point 1";
			
				if ( fits_create_tbl( outfptr, BINARY_TBL, lfr->fits_ip_nrows,  lfr->fits_ip_tfields, 
															ttypeIP, tformIP, tunitIP, extnameIP, &status) )
				{
					fits_report_error(stdout, status); 
					//return ERR_FITSIO;
				}
				if (fits_write_col(outfptr, TDOUBLE, 1, firstrow, firstelem, lfr->fits_ip_nrows, 
																				lfr->fits_ip, &status))
				{
					fits_report_error(stdout, status); 
					//return ERR_FITSIO;
				}
printf("status 5 %d\n",status);
				if ( fits_close_file(outfptr, &status) ) 
				{
					fits_report_error(stdout, status); 
					//return ERR_FITSIO;
				}		
				
				// attach it to the segment ---------------------------------------------------------------------
				//warning change segment number !!!
				segment_out = drms_segment_lookupnum (record_out, 1);
				status=drms_segment_write_from_file(segment_out,filenameout);

		//---------------------------------------
		//		5) close & free
		//---------------------------------------
		drms_close_record(record_out, DRMS_INSERT_RECORD);
		printf("end in the db...\n");							
	}
	else 
	{
		fprintf(info->opf,"NO LDFS file created (ifail >0) or retcode <0\n");
	}



return 0;
}
