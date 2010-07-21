// prog to simulate the call from the level1 or higher pipeline
// Thu Apr 22 10:25:20 HST 2010 //read data from drms, can use Richard' info as guess

/*

TODO:

*/

#include "fitsio.h"

#include "limbfit.h"
#define ERR_USAGE -1
#define ERR_FITSIO -20
#define ERR_SPECIAL -100

#define dsin "hmi.lev1c_nrt[]"
#define dsout "hmi_lf.lev1.n"

#include <jsoc_main.h>


ModuleArgs_t module_args[] = {
 {ARG_STRING, 	"ds", "hmi.lev1c_nrt[]", "data set"},
  {ARG_STRING, 	"dirout", "./", "output directory (default: ./)"},
  {ARG_INT, 	"fid", "0", "fid (no default)"},
  {ARG_INT, 	"debug", "1", "debug level (default: 1)"},
  {ARG_INT, 	"nb_ab_bin", "0", "NB AB bins"},
  {ARG_INT, 	"nb_ldf", "0", "NB LDFs"},
  {ARG_INT, 	"nb_rd_bin", "0", "NB radial bins"},
  {ARG_INT, 	"nb_iter", "1", "NB iterations"},
  {ARG_INT, 	"skipc", "0", "skip center determination"},
  {ARG_END}

};

char *module_name = "lfwrp";

int DoIt(void)
{
//************************************************************************
// usage ./lfwrp [0|1|2|3]  dstdir range_min range_max [-1 n] [-2 n] [...] [-11 n]
//************************************************************************
  CmdParams_t *params = &cmdparams;

	int status;
	int debug = params_get_int (params, "debug");
printf("DEBUG -%d-\n",debug);

	char *dirout = params_get_str (params, "dirout");
printf("Dir out -%s-\n",dirout);

	int ifid = params_get_int (params, "fid");
printf("selected FID -%d-\n",ifid);
	// changer ca en valmn-valmx ?

	int nb_ab_bin = params_get_int (params, "nb_ab_bin");
printf("selected NB_AB_BIN -%d-\n",nb_ab_bin);

	int nb_rd_bin = params_get_int (params, "nb_rd_bin");
printf("selected NB_RD_BIN -%d-\n",nb_rd_bin);

	int nb_ldf = params_get_int (params, "nb_ldf");
printf("selected NB_LDF -%d-\n",nb_ldf);

	int nb_iter = params_get_int (params, "nb_iter");
printf("selected NB_ITER -%d-\n",nb_iter);

	int skipc = params_get_int (params, "skipc");
printf("selected SKIPC -%d-\n",skipc);

	char *dsspec = params_get_str (params, "ds");
printf("ds -%s-\n",dsspec);
/***********************************************************************
	Initialization of the structure to pass to the limbfit routine
************************************************************************/
static LIMBFIT_INPUT limbfit_vars ;
static LIMBFIT_INPUT *lfv = &limbfit_vars;
static LIMBFIT_OUTPUT limbfit_res ;
static LIMBFIT_OUTPUT *lfr = &limbfit_res;
int retcode;

/***********************************************************************
	Open DRMS connexion, get the range of data, load structure and run lf
************************************************************************/
	DRMS_RecordSet_t *drs;
	DRMS_Record_t *record;
	DRMS_Segment_t *segment;

	static char filenameout[100];
	fitsfile *outfptr;
	static char filenameouttxt[100];
	FILE *opf;
	
    long naxes[2];
	int naxis_col;
	int naxis_row;
	
	int cam,fid,hpltid,missvals,maxval;
	float percentd;
	TIME tobs;
    char *hwltnset;
    char s_hwltnset[10];
	double x0_lf, y0_lf, rsun_lf;
	float sat_rot;//,inst_rot;
	long nelements, firstrow, firstelem;
	int i;
	char tmps[50];
	int rstatus, ncnt,r;
	unsigned int fsnx = 0;
	long long recnum0;
	int seg_ct;
	char s_tobs[64];

    drs = drms_open_records(drms_env, dsspec, &rstatus);
	if (!drs) {
		fprintf (stderr, "Error: unable to open record set %s\n", dsspec);
 	   return 0;
	}
	ncnt = drs->n;
	if (ncnt < 1) {
		fprintf (stderr, "No records in selected set %s\n", dsspec);
		return 0;
	}
	printf("fetch %d records\n",ncnt);
    for(r=0; r < ncnt; r++) 
    {
		record = drs->records[r];
		recnum0 = record->recnum;
		fsnx = drms_getkey_int(record, "FSN", &rstatus);       
		printf ("rec_num:%lld FSN:%u\n", recnum0, fsnx);
		//fprintf(opf,"rec_num:%lld FSN:%u\n", recnum0, fsnx);
		fid = drms_getkey_int(record, "FID", &rstatus);
		if(rstatus) {
			fprintf (stderr, "ERROR: in drms_getkey_int(FID) fid=%d\n", fid);
			//fprintf(opf,"ERROR: in drms_getkey_int(FID) fid=%d\n", fid);
		}
		//fprintf(opf,"FID: %d\n", fid);
		printf("FID: %d\n", fid);
		hwltnset = drms_getkey_string(record, "HWLTNSET", &rstatus);
		if(rstatus) {
			fprintf (stderr, "ERROR: in drms_getkey_int(HWLTNSET) hwltnset=%s\n", hwltnset);
			//fprintf(opf,"ERROR: in drms_getkey_int(HWLTNSET) hwltnset=%s\n", hwltnset);
		}
		missvals = drms_getkey_int(record, "MISSVALS", &rstatus);
		if(rstatus) {
			fprintf (stderr, "ERROR: in drms_getkey_int(MISSVALS) missvals=%d\n", missvals);
			//fprintf(opf,"ERROR: in drms_getkey_int(MISSVALS) missvals=%d\n", missvals);
		}
		sprintf(s_hwltnset,hwltnset);
		printf("HWLTNSET: -%s-\n", s_hwltnset);
		int a=(strncmp(hwltnset,"CLOS",4));
   		printf("cmp hwltnset fid %d %s %d\n",a,hwltnset,fid);
   		if (a==0 && (missvals==0) )
   		{
   			if ( (ifid==0) || (fid == ifid) )
   			{
				filenameouttxt[0] = '\0';
				strcat(filenameouttxt,dirout);
				sprintf(tmps, "%u", fsnx);
				strcat(filenameouttxt,tmps);
				strcat(filenameouttxt,"_res.txt");
				printf("-----------filenameouttxt %s\n",filenameouttxt);
		      	opf = fopen(filenameouttxt, "w");      
//only those needed by limbfit
				tobs = drms_getkey_time(record, "t_obs", &rstatus);
				if(rstatus) {
					fprintf(opf,"Error on drms_getkey_time() fsn=%u. Use DRMS_MISSING_TIME\n", fsnx);
					//tobs[i] = DRMS_MISSING_TIME;
				}
				cam = drms_getkey_int(record, "CAMERA", &rstatus);
				if(rstatus) {
					fprintf(opf,"ERROR: in drms_getkey_int(CAMERA) cam=%d\n", cam);
				}
				hpltid = drms_getkey_int(record, "hpltid", &rstatus);
				  if(rstatus) {
					fprintf(opf,"ERROR: in drms_getkey_int(HPLTID) hpltid=%d\n", hpltid);
				}
				percentd = drms_getkey_float(record, "PERCENTD", &rstatus);
				if(rstatus) {
					fprintf(opf,"ERROR: in drms_getkey_int(PERCENTD) percentd=%f\n", percentd);
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
				maxval = drms_getkey_int(record, "DATAMAX", &rstatus);
				if(rstatus) {
					fprintf (stderr, "ERROR: in drms_getkey_int(DATAMAX) maxval=%d\n", maxval);
					//fprintf(opf,"ERROR: in drms_getkey_int(MISSVALS) missvals=%d\n", missvals);
				}
				/*
				inst_rot = drms_getkey_float(record, "INST_ROT", &rstatus);
				  if(rstatus) {
					fprintf(opf,"ERROR: in drms_getkey_int(INST_ROT) inst_rot=%f\n", inst_rot);
				}
				*/
				seg_ct = drms_record_numsegments (record);
				if (seg_ct < 1) printf ("  no segments in selected record\n");
				printf ("  nb segments in selected record %d\n",seg_ct);
	//			record_segment = drms_segment_lookup (record, seg_name);
	//				if (!record_segment) fprintf (stderr, "Error, unable to locate segment %s of record %d\n",seg_name, recn);
					
				segment = drms_segment_lookupnum(record, 0);
				DRMS_Array_t *img;
				img = drms_segment_read(segment, DRMS_TYPE_FLOAT, &rstatus);
				if(!img) {
					fprintf(opf,"Can't do drms_segment_read() status=%d\n", rstatus);
					//return(1);
				} else
				{
					naxis_col=img->axis[0];
					naxis_row=img->axis[1];
					fprintf(opf,"  naxis %d %d %d\n",img->naxis,img->axis[0],img->axis[1]);
					lfv->data=img->data;
					sprint_time (s_tobs, tobs, "UT", 3);
					fprintf(opf,"T_OBS %s \n",s_tobs);
					fprintf(opf,"CAMERA %d \n",cam);
					fprintf(opf,"FID %d \n",fid);
					fprintf(opf,"FSN %u \n",fsnx);
					fprintf(opf,"HPLTID %d \n",hpltid);
					fprintf(opf,"PERCENTD %3.2f \n",percentd);
					fprintf(opf,"X0_LF %f \n",x0_lf);
					fprintf(opf,"Y0_LF %f \n",y0_lf);
					fprintf(opf,"RSUN_LF %f \n",rsun_lf);
					fprintf(opf,"HWLTNSET %s\n",hwltnset);
				
					/***********************************************************************
						Call the limbfit routine like build_lev1.c would do it
						- update the structure with parameters got from the command line
					***********************************************************************/
					fprintf(opf,"init vars...................................\n");
					lfv->fd=opf;			
					// init variables if not changed by argv
					if(nb_ab_bin == 0)
						lfv->nb_abb=NUM_AB_BINS;
					else
						lfv->nb_abb=nb_ab_bin;
					fprintf(opf," nb_abb %d\n",lfv->nb_abb);
					if(nb_rd_bin == 0)
						lfv->nb_rdb=NUM_RADIAL_BINS;
					else
						lfv->nb_rdb=nb_rd_bin;
					fprintf(opf," nb_rdb %d\n",lfv->nb_rdb);			
					if(nb_ldf == 0)
						lfv->nb_ldf=NUM_LDF;
					else
						lfv->nb_ldf=nb_ldf;
					fprintf(opf," nb_ldf %d\n",lfv->nb_ldf);
					if(nb_iter == 1)
						lfv->iter=1;
					else
						lfv->iter=nb_iter;
					fprintf(opf," nb_iter %d\n",lfv->iter);
		
					lfv->annw=ANNULUS_WIDTH;
						fprintf(opf," annw %d\n",lfv->annw);
					lfv->mxsz_annv=MAX_SIZE_ANN_VARS;
						fprintf(opf," mxsz_annv %ld\n",lfv->mxsz_annv);		
					lfv->lol=LO_LIMIT;
						fprintf(opf," lol %f\n",lfv->lol);
					lfv->upl=UP_LIMIT;
						fprintf(opf," upl %f\n",lfv->upl);
					lfv->incx=INC_X;
						fprintf(opf," incx %f\n",lfv->incx);
					lfv->incy=INC_Y;
						fprintf(opf," incy %f\n",lfv->incy);
					lfv->nb_fit_pnts=NUM_FITPNTS;
						fprintf(opf," nbpnts %d\n",lfv->nb_fit_pnts);
					lfv->guess_rng=GUESS_RANGE;
						fprintf(opf," range %d\n",lfv->guess_rng);
					lfv->img_szx=naxis_row;
						fprintf(opf," nax %d\n",lfv->img_szx);
					lfv->img_szy=naxis_col;
						fprintf(opf," nay %d\n",lfv->img_szy);			
					//take R' info as guess or real values (depending on limb_hmi.f)
					lfv->ir=rsun_lf;
					lfv->iy=y0_lf;
					lfv->ix=x0_lf;
					
					retcode=limbfit(lfv,lfr,skipc,debug);
					
					printf("returned code %d\n",retcode);
					fprintf(opf,"returned code %d\n",retcode);
					fprintf(opf,"after limbfit x %6.3f\n",lfr->cenx);
					fprintf(opf,"after limbfit y %6.3f\n",lfr->ceny);
					fprintf(opf,"after limbfit r %6.3f\n",lfr->radius);
					fprintf(opf,"after limbfit err %d\n",lfr->error);
					fprintf(opf,"after limbfit cmean %6.3f\n",lfr->cmean);
					fprintf(opf,"saveldf: %p \n",lfr->fits_ldfs_data); 
					fprintf(opf,"save_as: %p \n",lfr->fits_as); 
					fprintf(opf,"save_es: %p \n",lfr->fits_es); 
					fprintf(opf,"save_r: %p \n",lfr->fits_r); 
					fprintf(opf,"savel_ip: %p \n",lfr->fits_ip); 
					fprintf(opf,"save_a: %p \n",lfr->fits_ab_dataa); 
					fprintf(opf,"save_b: %p \n",lfr->fits_ab_datab); 
					//fprintf(opf,"save_b1: %p \n",lfr->fits_ab_datab1); 
					fprintf(opf,"guess x %6.3f\n",lfr->gxc);
					fprintf(opf,"guess y %6.3f\n",lfr->gyc);
					fprintf(opf,"guess r %6.3f\n",lfr->gr);
					if (retcode >= 0)
					{	
			
						/**********************************************************************
							Write the ldfs FITS file
						**********************************************************************/
							// for test purpose, later will be replaced by new segment written
							filenameout[0] = '\0';
							strcat(filenameout,dirout);
							//sprintf(tmps, "%u", fsnx);
							strcat(filenameout,tmps);
							strcat(filenameout,"_ldfs.fits");			
						
							firstrow=1;
							firstelem=1;
							status=0;
						
							fprintf(opf,"F1: -%s-\n",filenameout);	
						printf("F1: -%s-\n",filenameout);	
						
							remove(filenameout); 
							if (fits_create_file(&outfptr, filenameout, &status)) 
							{
								fits_report_error(stdout, status); 
								fclose(opf);
								//return ERR_FITSIO;
							}
							naxes[0]=lfr->fits_ext0_naxis1;
							naxes[1]=lfr->fits_ext0_naxis2;
							nelements=naxes[0]*naxes[1];
							if ( fits_create_img(outfptr, FLOAT_IMG, lfr->fits_ext0_naxis, naxes, &status) )
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
							if ( fits_update_key(outfptr, TSTRING, "T_OBS", &s_tobs, "Observing date (from original data file)",&status) )
							{
								fits_report_error(stdout, status); 
								//return ERR_FITSIO;
							}
							if (debug)
							{
								if ( fits_update_key(outfptr, TFLOAT, "PERCENTD", &percentd, "Percentage of good data",&status) )
								{
									fits_report_error(stdout, status); 
									//return ERR_FITSIO;
								}
								if ( fits_update_key(outfptr, TDOUBLE, "GUESS_X", &lfr->gxc, "Guess Center Position X",&status) )
								{
									fits_report_error(stdout, status); 
									//return ERR_FITSIO;
								}
								if ( fits_update_key(outfptr, TDOUBLE, "GUESS_Y", &lfr->gyc, "Guess Center Position Y",&status) )
								{
									fits_report_error(stdout, status); 
									//return ERR_FITSIO;
								}
								if ( fits_update_key(outfptr, TDOUBLE, "GUESS_R", &lfr->gr, "Guess Solar Radius",&status) )
								{
									fits_report_error(stdout, status); 
									//return ERR_FITSIO;
								}
							}
							
							if (isnan(x0_lf)) { x0_lf=0.0;}
							if ( fits_update_key(outfptr, TDOUBLE, "X0_LF", &x0_lf, "Center Position X from lev1",&status) )
							{
								fits_report_error(stdout, status); 
								//return ERR_FITSIO;
							}
							if (isnan(y0_lf)) { y0_lf=0.0;}
							if ( fits_update_key(outfptr, TDOUBLE, "Y0_LF", &y0_lf, "Center Position Y from lev1",&status) )
							{
								fits_report_error(stdout, status); 
								//return ERR_FITSIO;
							}
							if (isnan(rsun_lf)) { rsun_lf=0.0;}
							if ( fits_update_key(outfptr, TDOUBLE, "RSUN_LF", &rsun_lf, "Solar Radius from lev1",&status) )
							{
								fits_report_error(stdout, status); 
								//return ERR_FITSIO;
							}
							if ( fits_update_key(outfptr, TINT, "CAMERA", &cam, "Camera ID",&status) )
							{
								fits_report_error(stdout, status); 
								//return ERR_FITSIO;
							}
							if ( fits_update_key(outfptr, TINT, "FID", &fid, "FID from lev1",&status) )
							{
								fits_report_error(stdout, status); 
								//return ERR_FITSIO;
							}
							if ( fits_update_key(outfptr, TINT, "FSN", &fsnx, "FSN from lev1",&status) )
							{
								fits_report_error(stdout, status); 
								//return ERR_FITSIO;
							}
							if ( fits_update_key(outfptr, TINT, "HPLTID", &hpltid, "HPLTID from lev1",&status) )
							{
								fits_report_error(stdout, status); 
								//return ERR_FITSIO;
							}
							if ( fits_update_key(outfptr, TSTRING, "HWLTNSET", &s_hwltnset, "HWLTNSET from lev1",&status) )
							{
								fits_report_error(stdout, status); 
								//return ERR_FITSIO;
							}
							if ( fits_update_key(outfptr, TFLOAT, "SAT_ROT", &sat_rot, "SAT_ROT from lev1",&status) )
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
							if ( fits_update_key(outfptr, TINT, "DATAMAX", &maxval, "Maximum value of all pixels ",&status) )
							{
								fits_report_error(stdout, status); 
								//return ERR_FITSIO;
							}
								if ( fits_write_img(outfptr, TFLOAT, 1, nelements, lfr->fits_ldfs_data, &status) ) 
								{
									fits_report_error(stdout, status); 
									//return ERR_FITSIO;
								}	
							
								status=0;
							
								// As ---EXT#1-------------------------------------------------------------------------
								char extnameA[] = "LDF info / As";
								char *tunitAE[] = { "", "", "", "", "", "", "" };
								char *tformAE[] = { "1D", "1D", "1D", "1D", "1D", "1D" };
								char *ttypeA[]  = { "A0", "A1", "A2", "A3", "A4", "A5"};
							
								if ( fits_create_tbl(outfptr, BINARY_TBL, lfr->fits_ext1_nrows, lfr->fits_ext1_tfields, 
																	ttypeA, tformAE,tunitAE, extnameA, &status) )
								{
										fits_report_error(stdout, status); 
										//return ERR_FITSIO;
								}
								for (i=1;i<=lfr->fits_ext1_tfields;i++)
								{
									if (fits_write_col(outfptr, TDOUBLE, i, firstrow, firstelem, lfr->fits_ext1_nrows, 
																	lfr->fits_as+(i-1)*lfr->fits_ext1_nrows, &status))
									{
										fits_report_error(stdout, status); 
										//return ERR_FITSIO;
									}
								}
						
								// Es ----EXT#2------------------------------------------------------------------------
								char *ttypeE[] = { "E0", "E1", "E2", "E3", "E4", "E5"};
								char extnameE[] = "LDF info / Es";
							
								if ( fits_create_tbl(outfptr, BINARY_TBL, lfr->fits_ext2_nrows, lfr->fits_ext2_tfields, 
																		ttypeE, tformAE, tunitAE, extnameE, &status) )
								   {
										fits_report_error(stdout, status); 
										//return ERR_FITSIO;
									}
								for (i=1;i<=lfr->fits_ext2_tfields;i++)
								{
									if (fits_write_col(outfptr, TDOUBLE, i, firstrow, firstelem, lfr->fits_ext2_nrows, 
																	(lfr->fits_es+(i-1)*lfr->fits_ext2_nrows), &status))
									{
										fits_report_error(stdout, status); 
										//return ERR_FITSIO;
									}
								}
							
								// Radius ---EXT#3---------------------------------------------------------------------
								char *tunitR[] = { "" };
								char *ttypeR[] = { "Radius"};
								char *tformR[] = { "1D" };
								char extnameR[] = "LDF info / Radius";
							
								if ( fits_create_tbl(outfptr, BINARY_TBL, lfr->fits_ext3_nrows, lfr->fits_ext3_tfields, 
																			ttypeR, tformR, tunitR, extnameR, &status) )
								{
									fits_report_error(stdout, status); 
									//return ERR_FITSIO;
								}
								if (fits_write_col(outfptr, TDOUBLE, 1, firstrow, firstelem, lfr->fits_ext3_nrows, 
																								lfr->fits_r, &status))
								{
									fits_report_error(stdout, status); 
									//return ERR_FITSIO;
								}
							
								// AB ---EXT#4---------------------------------------------------------------------
								char *tunitAB[] = { "","" };
								char *ttypeAB[] = { "Alpha", "Beta"};
								char *tformAB[] = { "1E","1E" };
								char extnameAB[]="LDF info / Alpha&Beta";
								if ( fits_create_tbl( outfptr, BINARY_TBL, lfr->fits_ext4_nrows, lfr->fits_ext4_tfields, 
																		ttypeAB, tformAB, tunitAB, extnameAB, &status) )
								{
									fits_report_error(stdout, status); 
									//return ERR_FITSIO;
								}
								if (fits_write_col(outfptr, TFLOAT, 1, firstrow, firstelem, lfr->fits_ext4_nrows, 
																						lfr->fits_ab_dataa, &status))
								{
									fits_report_error(stdout, status); 
									//return ERR_FITSIO;
								}
								if (fits_write_col(outfptr, TFLOAT, 2, firstrow, firstelem, lfr->fits_ext4_nrows,
																						lfr->fits_ab_datab, &status))
								{
									fits_report_error(stdout, status); 
									//return ERR_FITSIO;
								}
		
								// IP ---EXT#5---------------------------------------------------------------------
								char *tunitIP[] = { "" };
								char *ttypeIP[] = { "IP1"};
								char *tformIP[] = { "1D" };
								char extnameIP[]="LDF info / Inflection Point 1";
							
								if ( fits_create_tbl( outfptr, BINARY_TBL, lfr->fits_ext5_nrows,  lfr->fits_ext5_tfields, 
																			ttypeIP, tformIP, tunitIP, extnameIP, &status) )
								{
									fits_report_error(stdout, status); 
									//return ERR_FITSIO;
								}
								if (fits_write_col(outfptr, TDOUBLE, 1, firstrow, firstelem, lfr->fits_ext5_nrows, 
																								lfr->fits_ip, &status))
								{
									fits_report_error(stdout, status); 
									//return ERR_FITSIO;
								}
								//--free local vars
		//						free(extnameA);
		//						free(tunitAE);
		//						free(tformAE);
		//						free(ttypeA);
		//						free(ttypeE);
		//						free(extnameE);
		//						free(tunitR);
		//						free(ttypeR);
		//						free(tformR);
		//						free(extnameR);
		//						free(tunitAB);
		//						free(ttypeAB);
		//						free(tformAB);
		//						free(extnameAB);
		//						free(tunitIP);
		//						free(ttypeIP);
		//						free(tformIP);
		//						free(extnameIP);
								//----------------------
						
							//-----------------------------------------------------------------------
								if ( fits_close_file(outfptr, &status) ) 
								{
									fits_report_error(stdout, status); 
									//return ERR_FITSIO;
								}		
							
							}
							else 
							{
								printf("NO LDFS file created (ifail >0) or retcode <0\n");
								fprintf(opf,"NO LDFS file created (ifail >0) or retcode <0\n");
							}
					//free(hwltnset);
					drms_free_array(img);
				}
				fclose(opf);
			} // if fid ok
			else printf("REJECTED:(1) rec_num:%lld FSN:%u\n", recnum0, fsnx);
		} // if loop ok & no missvals
		else printf("REJECTED:(2) rec_num:%lld FSN:%u\n", recnum0, fsnx);
	} //end process one DRMS record : for(r=0; r < ncnt; r++) 
	drms_free_record(record);
		
/**********************************************************************/
//free(data); 
//free(filenameout);

//+ free all variables coming from limbfit.c

free(lfr->fits_ldfs_data);
free(lfr->fits_as);
free(lfr->fits_es);
free(lfr->fits_r);
free(lfr->fits_ab_dataa);
free(lfr->fits_ab_datab);
free(lfr->fits_ip);

printf("end of main\n");

return(retcode);
}

/*
void get_keys(strkey *gsk)
{

}

void write_fits(STRKEYS *gsk,LIMBFIT_OUTPUT *results)
{

}
*/

