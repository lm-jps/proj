#ifndef LEV0LEV1_INCL
#define LEV0LEV1_INCL 1

//Used by build_lev1.c to pass info to do_flat().
//All the data arrays are filled in except adata1,
//which do_flat() populates.
typedef struct {
  DRMS_Record_t *rs0;	//drms lev0 record
  DRMS_Record_t *rs1;	//drms lev1 record
  DRMS_Record_t *rsff;	//drms flat field record
  short *adata0;	//lev0 segment array data
  int *adata1;		//lev1 segment array data
  float *adataff;	//flat field array data
  float *adatadark;	//bias dark array data
  int *adatabad;	//bad pixel array data
  long long recnum0;	//lev0 record DRMS record number
  long long recnum1;	//lev1 record DRMS record number
  unsigned int fsn;	//fsn of lev0 record
  int himgcfid;		//HMI_SEQ_ID_IMAGE_CNFG
  int datamin;
  int datamax;
  int datamedn;
  double datamean;
  double data_rms;
  double dataskew;
  double datakurt;
  int datavals;
  int missvals;
} LEV0LEV1;

#endif
