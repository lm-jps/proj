/*
 * render_image - make .png or .ppm image from a segment image file in a series.
 */

/**
   @defgroup render_image render_image  Make .png or .ppm image from a segment image file in a series.
   @ingroup su_util

   @brief Make .png or .ppm image from a segment image file in a series.

   @par Synopsis:
   @code
   render_image --help
   render_iamge in=<RecordSet> out=<target> | {outname=<ident_for_filename>}  {scaling=<scaletype>}
                {palette=<color_table>}  {type=<image_protocol>} {outid=<filename_code>}
                {tkey=<time_keyword>} {min=<scale_min>} {max=<scale_max>}
                {scale=<image_size_ratio>} {-c} {-u} {-w} {-x}
   @endcode

   @code
     in=Input data series.
     n=limit limit number of records to n from beginning or -n from end of 'in'
     out=Output data series or full path of directory or pipe via ppm program.
     outname=Output quantity name for filenames.
     scaling=Color table type, minmax, minmax99, minmaxgiven, histeq, ...,, default=MINMAX
     palette=Color table, GREY or filename, default=GREY (bad spelling as pallette still allowed)
     type=Image protocol, png, ppm, ..., default=png
     tkey=keyword name for time, default T_REC
     outid=output identifier, #=record counter, time=time as yyyymmdd_hhmmss, default=#
     min=Min value for scaling, default nan
     max=Max value for scaling, default nan
     scale"Reduction factor list, default 1.
     -c Crop flag, causes crop to RSUN_OBS
     -u use unchanged for rotation and centering
     -w use white for missing pixels, instead of black
     -x High-quality flag, sets bytespercolor to 2 instead of 1,
  @endcode

   @par
   Render_image can make one or many images per record in the input recordset.  It will make one image for each value
   of scale where the values are given as a comma delimited list of divisors for the dimensions of the input image array.
   The image file will be written to a filename created from outid_outname_type where outid is either "time" or "#",
   outname is any literal string, and type is "png" or "ppm" or if out specifies a pipe, then type may be any file type
   that is supported as the output of a pipeline which accepts a PPM file.
   @par
   The file location is spcified with the out parameter.  The target may be a seriesname in which case the image will
   be written into a record directory of that series.  If the target seriesname is the same as the series specified in in, then
   each input record will be cloned.  The image will be written into a segment named "image".
   @par
   If the target location begins with a "/" or "./" then the $target will be used as the path for a directory into which
   the image file will be written.  If the target location beginw with a "|" character, then the target string will be used
   as a pipe into which a .ppm file will be written.  The pipe string will have " > <filename>" appended to it.
   In the case of a pipe, the string given in the 'out' parameter will be scanned
   for replacement tokens delimited by '{...}'.  Inside the '{}' pair the following
   replacement instructions are available:
  @code
   The word 'ID' is replaced by the filename ID as built from outid.
   The word '%' is replaced by that percent of the height of the image rounded
       to the nearest pixel.  If the floating point number after the '%' is itself
       followed by a ':' and an integer, then that integer will be used as the
       smallest number for the replacement token.  So, for instance '{%1.2:5}'
       with the image height of 256 pixels will be processed as: 1.2% of 256 is 3.072
       which rounds to 3 which is smaller than 5 so the result will be 5.
  @endcode
   @par
   The image filename is generated from the outid, outname, and type parameters.  If outid is the string "time" then
   the time of the image as specified by the tkey (defaults to T_REC) will be used with the "." and ":" and "_<zone>" componenets deleted.  I.e.
   the time will be in the form "yyyymmdd_hhmmss".  The outid = time can have a suffix specifying the number of time characters to
   include.  the format of the suffix is ":nn" where nn is from 1 to 15. If outid is a literal "#" then the "%04d" layout will be used to generate
   the image sequence number within the run of the program (i.e. the record number within the current recordset).  the
   type parameter must be "png" or "ppm" unless the out parameter specifies a ppm pipeline in which case type should
   be the suffix for the file type geenrated by the pipe.
   @par
   The scaling of image values to color table indices is controlled by the scaling parameter.  scaling must be one of
   MINMAX, MAXMIN, MINMAXGIVEN, MINMAX99, MAXMIN99, HISTEQ, MAG at present.  The default is MINMAX.  If the scaling is
   MINMAX and both min and max parameters are given, then the MINMAXGIVEN scaling will be used.  MINMAX99 and
   MAXMIN99 mean to eliminate the top and bottom 0.5% of the histogram of the values before computing a linear scaling.
   MAXMIN and MAXMIN99 result in a reversed scaling with the max value mapped to the first color in the table and the
   min value mapped to the max color in the table.  String case is ignored for the scaling parameter.
   @par
   The MAG scaling uses a scaling of 400*M/pow((abs(M)+bias),0.75)  bias=30, from user provided min and max (defaulting to +-1500).
   This scaling emphasizes low absolute field strength but still shows detail up to a max of 3000 gauss.  The default
   of 1500 clips umbral fields with the tradeoff of more contrast for plage and network fields.
   @par
   The color table is specified by the palette parameter.  The table may be a built in table, at persent the only
   build-in is "grey".  ("grey", "gray", "GREY", "GRAY" are all allowed).  If not a build-in table name the palette
   parameter is expected to be a file path to a table of .sao or .lut types.  These table formats are described in the
   ds9(1) documentation.
   @par
   The image will be centered to the nearest pixel, and will be rotated 180 degrees if CROTA2 is near 180.

   @par Flags:
   @c -c:  Crop image.  This flag will casue the image to be cropped at pixel radius RSUN_OBS/CDELT1 about CRPIX1,CRPIX2 center.

   @c -w:  White background, the -w flag will cause missing values to be white.  Else black.

   @c -x:  the -x flag is used to indicate extra resolution in the color table, i.e. 16-bits per color.  

   @par Example to generate a directory full of simple grey scaled png images, perhaps for a movie
   @code
   render_image in='hmi.M_45s_nrt[2010.06.01_TAI/5d@30m]' out=./movie scaling=minmax99 outid=# outname=M 
   @endcode

   @par Example to generate a set of 4096x4096 and 1024x1024, and 256x256 pixel images of magnetograms in the current directory
   using a color table the adds dynamic range to small absolute values.  The times will omit the seconds field in this example.
   @code
   render_image in='hmi.M_45s_nrt'$QRY \
     outname=M \
     palette=/home/phil/apps/mag.lut \
     outid=time:13 \
     -c \
     min=-500 \
     max=500 \
     type=jpg \
     -w \
     scale=1,4,16 \
     out='| ppmlabel -color black -size {%0.75:5} -x 15 -y {%98} -text "SDO/HMI Quick-Look Magnetogram: {ID}" | pnmtojpeg -quality=95'
   @endcode

  @bug
  Using the MAG scaling with a list of scale sizes, any instance of scale=1 must appear last in the list.  This is because the values of the
  image data are replaced in the original array in the scale=1 case.  This will affect any scaling type that modifies the values of the data.

  @bug
  Does not make 16-bit colors correctly, or I can not get ppm tools with correct flags...

*/

#include "mypng.h"
#include "jsoc_main.h"

#define PNGDIE(msg) {fprintf(stderr,"%s\n",msg); \
	if (png_ptr) png_destroy_write_struct(&png_ptr, &info_ptr); \
	if (row_pointers) free(row_pointers); \
	if (fp) fclose(fp); \
	return(1);}

#define DIE(msg) {fflush(stdout);fprintf(stderr,"%s, status=%d\n",msg,status); return(status);}

#define MINMAX	1
#define MAXMIN	2
#define HISTEQ  3
#define MINMAX99  4
#define MAXMIN99  5
#define MINMAXGIVEN  6
#define MAG 7
#define MAGMINMAXGIVEN 8
#define LOG 9
#define SQRT 10

char *module_name = "render_image";

ModuleArgs_t module_args[] =
{
     {ARG_STRING, "in", "NOT SPECIFIED",  "Input data series."},
     {ARG_INT, "n", "0",  "Limit of number of records, 0 for no limin; >0 count from start; <0 count from end"},
     {ARG_STRING, "out", "NOT SPECIFIED",  "Output data series or full path of directory or pipe via ppm program."},
     {ARG_STRING, "outname", "NOT SPECIFIED",  "Output quantity name for filenames"},
     {ARG_STRING, "scaling", "MINMAX", "Color table type, minmax, minmax99, minmaxgiven, histeq, ..."},
     {ARG_STRING, "pallette", "none", "Color table, GREY or filename"},
     {ARG_STRING, "palette", "GREY", "Color table, GREY or filename"},
     {ARG_STRING, "type", "png", "Image protocol, png, ppm, ..."},
     {ARG_STRING, "outid", "#", "output identifier, #=record counter, time=time as yyyymmdd_hhmmss"},
     {ARG_STRING, "tkey", "T_REC", "Time keyword name used if outid==time, defaults to T_REC"},
     {ARG_FLOAT, "min", "DRMS_MISSING_FLOAT", "Min value for scaling"},
     {ARG_FLOAT, "max", "DRMS_MISSING_FLOAT", "Max value for scaling"},
     {ARG_FLOAT, "bias", "30.0", "Max value for mag scaling"},
     {ARG_INTS, "scale", "1", "Reduction factors"},
     {ARG_FLAG, "c", "0", "Crop flag, causes crop to RSUN_OBS"},
     {ARG_FLAG, "u", "0", "use unchanged, is-is for rotation and centering"},
     {ARG_FLAG, "w", "0", "use white for missing pixels, instead of black"},
     {ARG_FLAG, "x", "0", "High-quality flag, sets bytespercolor to 2 instead of 1"},
     {ARG_END}
};

#define     Deg2Rad    (M_PI/180.0)
#define     Rad2arcsec (3600.0/Deg2Rad)
#define     arcsec2Rad (Deg2Rad/3600.0)
#define     Rad2Deg    (180.0/M_PI)

struct ObsInfo_struct
  {
  // from observation info
  TIME  t_obs;
  double rsun_obs, obs_vr, obs_vw, obs_vn;
  double crpix1, crpix2, cdelt1, cdelt2, crota2;
  double crval1, crval2;
  double cosa, sina;
  double obs_b0;
  // observed point
  int i,j;
  // parameters for observed point
  double x,y,r;
  double rho;
  double lon;
    double lat;
  double sinlat, coslat;
  double sig;
  double mu;
  double chi;
  double obs_v;
  };

typedef struct ObsInfo_struct ObsInfo_t;

void rebinArraySF(DRMS_Array_t *out, DRMS_Array_t *in);

int add_png(char *imgPath, DRMS_Array_t *imgArray, int scaletype, int usewhite, char *palette, int colors, int bytespercolor, double *minp, double *maxp);

int add_ppm(char *imgPath, DRMS_Array_t *imgArray, int scaletype, int usewhite, char *palette, int colors, int bytespercolor, double *minp, double *maxp);

int make_png(char *imgPath, unsigned char *data, int height, int width, char *palette, int bytepercolor, int colors);

char *set_scaling(DRMS_Array_t *in, double *minp, double *maxp, int *nmissp, char *palette, int colors, int bytepercolor, int scaling, int usewhite);

int upNcenter(DRMS_Array_t *arr, ObsInfo_t *ObsLoc);
int crop_image(DRMS_Array_t *arr, ObsInfo_t *ObsLoc);

ObsInfo_t *GetObsInfo(DRMS_Segment_t *seg, ObsInfo_t *pObsLoc, int *rstatus);

#define PNG 1
#define PPM 2
#define PPMPIPE 3

int DoIt(void)
  {
  int status = DRMS_SUCCESS;
  DRMS_RecordSet_t *inRS, *outRS = NULL;
  DRMS_Record_t *inRec, *outRec;
  int irec, nrecs;
  int fileonly=0;
  int imgtype;
  int quality;
 
  char *inQuery = (char *)params_get_str(&cmdparams, "in");
  char *outQuery = (char *)params_get_str(&cmdparams, "out");
  char *outName = (char *)params_get_str(&cmdparams, "outname");
  char *outid = (char *)params_get_str(&cmdparams, "outid");
  char *tkey = (char *)params_get_str(&cmdparams, "tkey");
  char *palette = (char *)params_get_str(&cmdparams, "palette");
  if (strcmp((char *)params_get_str(&cmdparams, "pallette"), "none"))
    palette = (char *)params_get_str(&cmdparams, "pallette");
  int limit = params_get_int(&cmdparams, "n");
  int bytespercolor = (params_isflagset(&cmdparams, "x") ? 2 : 1);
  int missingwhite = params_isflagset(&cmdparams, "w");
  int crop = params_isflagset(&cmdparams, "c");
  int asis = params_isflagset(&cmdparams, "u");
  char *type = (char *)params_get_str(&cmdparams, "type");
  char *scaling = (char *)params_get_str(&cmdparams, "scaling");
  double called_min = params_get_double(&cmdparams, "min");
  double called_max = params_get_double(&cmdparams, "max");
  double min, max;
  int *scales = NULL;
  int iscale, nscales;
  int colors;
  int scaletype;
  int srcNx = 0, srcNy = 0;

  if (strcasecmp(palette, "grey") == 0 || strcasecmp(palette, "gray") == 0)
    {
    palette = "GREY";
    colors = 1;
    }
  else
    colors = 3;
    
  if (strcasecmp(type, "png")==0) imgtype = PNG;
  else if (strcasecmp(type, "ppm")==0) imgtype = PPM;
  else imgtype = PPMPIPE;

  nscales = cmdparams_get_intarr(&cmdparams, "scale", &scales, &status);
  if (status)
    DIE("scale parameter not found");

// for (iscale=0; iscale<nscales; iscale++)
// fprintf(stderr,"Scale #%d=%d ", iscale, scales[iscale]);
// fprintf(stderr," using min=%f, max=%f\n", min, max);

  if (strcasecmp(scaling, "minmax") == 0) scaletype = MINMAX;
  else if (strcasecmp(scaling, "maxmin") == 0) scaletype = MAXMIN;
  else if (strcasecmp(scaling, "minmax") == 0) scaletype = MINMAX99;
  else if (strcasecmp(scaling, "maxmin99") == 0) scaletype = MAXMIN99;
  else if (strcasecmp(scaling, "histeq") == 0) scaletype = HISTEQ;
  else if (strcasecmp(scaling, "minmaxgiven") == 0) scaletype = MINMAXGIVEN;
  else if (strncasecmp(scaling, "mag", 3) == 0) scaletype = MAG;
  else if (strncasecmp(scaling, "log", 3) == 0) scaletype = LOG;
  else if (strncasecmp(scaling, "sqrt", 4) == 0) scaletype = SQRT;

  else  scaletype = 0;
  fprintf(stdout, " using scaling %d\n", scaletype);
  if (isnan(called_min) || isnan(called_max))
    {
    if (scaletype == MINMAXGIVEN)
      DIE("min and max must be specified for type MINMAXGIVEN\n");
    }
  else
    {
    if (scaletype == MINMAX)
      scaletype = MINMAXGIVEN;
    else if (scaletype == MAG)
      scaletype = MAGMINMAXGIVEN;
    }

  if (limit == 0)
    inRS = drms_open_records(drms_env, inQuery, &status);
  else
    inRS = drms_open_nrecords(drms_env, inQuery, limit, &status);
  if (status || inRS->n == 0)
       DIE("No input data found");
  drms_stage_records(inRS, 1, 1);
  nrecs = inRS->n;
fprintf(stderr,"nrecs=%d\n",nrecs);

  if (!outQuery[0])
       DIE("No valid output place specified");
  if (outQuery[0] == '/' || (outQuery[0] == '.' && outQuery[1] && outQuery[1] == '/') || outQuery[0] == '|')
    {
    // Output to a directory
    fileonly = 1;
    }
  else
    {
    // output to a record segment
    fileonly = 0;
    if (strcmp(inRS->records[0]->seriesinfo->seriesname, outQuery) == 0)
      outRS = drms_clone_records(inRS, DRMS_PERMANENT, DRMS_SHARE_SEGMENTS, &status);
    else
      outRS = drms_create_records(drms_env, nrecs, outQuery, DRMS_PERMANENT, &status);
    if (status)
      DIE("Output recordset not created");
    }

  for (irec=0; irec<nrecs; irec++)
    {
    int status;
    DRMS_Array_t *srcArray;
    DRMS_Segment_t *srcSeg;
    ObsInfo_t *ObsLoc;
 
    char imageID[100];
    inRec = inRS->records[irec];
    quality = drms_getkey_int(inRec, "QUALITY", &status);
    if (!status && quality < 0)
      {
      fprintf(stderr,"Quality bad for rec %d\n", irec);
      continue;
      }
    else
      {
      srcSeg = drms_segment_lookupnum(inRec, 0);
      srcArray = drms_segment_read(srcSeg, DRMS_TYPE_FLOAT, &status);
      if (!srcArray || status)
        {
         fprintf(stderr," No data file found rec %d, status=%d\n", irec,status);
         if (srcArray)
           {
           drms_free_array(srcArray);
           srcArray = NULL;
           }
         continue;
        }
      ObsLoc = GetObsInfo(srcSeg, NULL, &status);
      if (!asis) upNcenter(srcArray, ObsLoc);
      if (crop) crop_image(srcArray, ObsLoc);
      srcNx = srcArray->axis[0];
      srcNy = srcArray->axis[1];
      }
    if (strcmp(outid,"#") == 0)
      sprintf(imageID, "%04d", irec);
    else if (strncasecmp(outid,"time",4) == 0)
      {
      int timestatus;
      char timebuf[100];
      static int t[] = {0,1,2,3,5,6,8,9,10,11,12,14,15,17,18};
      int i;
      int maxtplace = 15;
      char *tplace = index(outid,':');
      TIME now = drms_getkey_time(inRec, tkey, &timestatus);
      if (timestatus)
        {
        fprintf(stderr,"Time key given in 'tkey' param, %s, is not found in the input series.\n",tkey);
        return(timestatus);
        }
      sprint_time(timebuf, now, "TAI", 0);
      if (tplace)
        sscanf(tplace+1,"%d",&maxtplace);
      for (i=0; i<maxtplace; i++)
        imageID[i] = timebuf[t[i]];
      imageID[i] = '\0';
      }
    else
      {
      fprintf(stderr, "Invalid outid param: |%s|\n", outid);
      return(1);
      }
    for (iscale=0; iscale<nscales; iscale++)
      {
      int dimxy;
      char dimxytxt[50];
      DRMS_Array_t *imageArray = NULL;
      char imgPath[DRMS_MAXPATHLEN];
      int scale = scales[iscale];
      int imageDims[2];
      int tmparray = 0;
      char imageFileName[DRMS_MAXSEGFILENAME];
      min = called_min;
      max = called_max;
      imageDims[0] = srcNx/scale;
      imageDims[1] = srcNy/scale;
      dimxy = imageDims[0];
      sprintf(dimxytxt,"%d",dimxy);
// fprintf(stderr,"scale=%d, imageNy=%d, imageNx=%d\n", scale, imageDims[1], imageDims[0]);

      if (scale != 1)
        {
        imageArray = drms_array_create(DRMS_TYPE_FLOAT, 2, imageDims, NULL, &status);
        rebinArraySF(imageArray, srcArray);
        tmparray = 1;
        }
      else
        {
        imageArray = srcArray;
        tmparray = 0;
        }

      sprintf(imageFileName, "%s%s%s_%s.%s", imageID, (strcmp(outid,"#")==0 ? "." : "_"), outName,
                                                  (dimxy == 4096 ? "4k" :
                                                  (dimxy == 2048 ? "2k" :
                                                  (dimxy == 1024 ? "1k" :
                                                  (dimxy == 512 ? "512" :
                                                  (dimxy == 256 ? "256" :
                                                  (dimxy == 128 ? "128" : dimxytxt)))))), type);
      
      if (fileonly)
        {
        if (outQuery[0] == '|')
          { // process outQuery for '{...}' replacement tokens
          char buf[4096];
          char *pbuf = buf, *c = outQuery;
          while (*c && (pbuf-buf) < sizeof(buf))
            {
            if (*c == '{')
              {
              c++;
              if (*c == '%')
                {
                float x;
                int ix, minx;
                c++;
                x = strtof(c, &c);
                ix = round(0.01*x*imageDims[1]);
                if (*c == ':')
                  {
                  c++;
                  minx = (int)strtol(c, &c, 10);
                  }
                else minx = 0;
                pbuf += sprintf(pbuf,"%d",(ix < minx ? minx : ix));
                }
              else if (strncmp(c, "ID", 2) == 0)
                {
                pbuf += sprintf(pbuf, "%s", imageID);
                c += 2;
                }
             if (*c++ != '}')
               DIE("Unrecognized token in pipe string");
             }
           else
             *pbuf++ = *c++;
           }
           *pbuf = '\0';
          sprintf(imgPath, "%s > %s", buf, imageFileName);
          }
        else
          sprintf(imgPath, "%s/%s", outQuery, imageFileName);
        }
      else
        {
        DRMS_Segment_t *imageSeg;
        char recdir[DRMS_MAXPATHLEN];
        outRec = outRS->records[irec];
        drms_copykeys(outRec, inRec, 0, kDRMS_KeyClass_Explicit);
        drms_record_directory(outRec, recdir, 0);
        imageSeg = drms_segment_lookup(outRec, "image");
        strcpy(imageSeg->filename, imageFileName);
        drms_keyword_setdate(outRec);
        sprintf(imgPath, "%s/%s", recdir, imageFileName);
        }

      if (imgtype == PNG)
        add_png(imgPath, imageArray, scaletype, missingwhite, palette, colors, bytespercolor, &min, &max);
      else if (imgtype == PPM || imgtype == PPMPIPE)
        add_ppm(imgPath, imageArray, scaletype, missingwhite, palette, colors, bytespercolor, &min, &max);

      if (tmparray)
        drms_free_array(imageArray);
      }
    drms_free_array(srcArray);
    }
fprintf(stderr,"render_image done, irec=%d\n",irec);
  drms_close_records(inRS, DRMS_FREE_RECORD);
fprintf(stderr,"render_image close in records done\n");
  if (!fileonly && outRS)
    drms_close_records(outRS, DRMS_INSERT_RECORD);
fprintf(stderr,"leaving module\n");
  return (DRMS_SUCCESS);
  } // end of DoIt

/* ----------------------------------------------------------------------
                    functions
      ------------------------------------------------------------
*/ 

int add_ppm(char *imgPath, DRMS_Array_t *imageArray, int scaletype, int usewhite, char *palette, int colors, int bytespercolor, double *minp, double *maxp)
  {
  int i, Nx, Ny;
  unsigned char *imgData = NULL;
  int nmissing;
  int maxpixval = (bytespercolor == 1 ? 255 : 65535);
  FILE *out;
  char *ppmcode;

  imgData = (unsigned char*)set_scaling(imageArray, minp, maxp, &nmissing, palette, colors, bytespercolor, scaletype, usewhite);
  if (!imgData)
    {
    fprintf(stderr,"scaling failed\n");
    return(1);
    }

  ppmcode = (colors == 1 ? "P5" : "P6");

  Nx = imageArray->axis[0];
  Ny = imageArray->axis[1];

  if (imgPath[0] == '|')
    out = popen(imgPath+1, "w");
  else
    out = fopen(imgPath, "w");

  fprintf(out, "%s\n%d %d %d\n", ppmcode, Nx, Ny, maxpixval); 
  for (i=Ny-1; i >=0; i--)
      fwrite(imgData + i*Nx*bytespercolor*colors, bytespercolor, colors*Nx, out);
  if (imgPath[0] == '|')
    pclose(out);
  else
    fclose(out);

  free(imgData);
  return(0);
  }

int add_png(char *imgPath, DRMS_Array_t *imageArray, int scaletype, int usewhite, char *palette, int colors, int bytespercolor, double *minp, double *maxp)
  {
  int srcNx, srcNy;
  int status=0;
  unsigned char *imgData = NULL;
  int nmissing;

  imgData = (unsigned char*)set_scaling(imageArray, minp, maxp, &nmissing, palette, colors, bytespercolor, scaletype, usewhite);

  if (!imgData)
    {
      fprintf(stderr,"scaling failed, status=%d\n",status);
      return(status);
    }

// should add text info to .png file with min, max, nmissing, etc.

  if (make_png(imgPath, imgData, imageArray->axis[0], imageArray->axis[1], palette, bytespercolor, colors) != 0)
    {
      fprintf(stderr,"png failed, status=%d\n",status);
      free(imgData);
      return(status);
    }
  free(imgData);
  return(0);
  }
// ---------------------------------------------------------------------- 

void rebinArraySF(DRMS_Array_t *out, DRMS_Array_t *in)
  {
    float *inData = in->data;
    float *outData = out->data;
    float *inp, *outp;
    int inNx=in->axis[0], inNy=in->axis[1];
    int outNx=out->axis[0], outNy=out->axis[1];
    int inI, inJ, outI, outJ, i;
    int binsize = inNx/outNx;
    double inRow[inNx], outRow[outNx];
    int inRowN[inNx], outRowN[outNx];

    for (outJ=0; outJ < outNy; outJ++)
      {
        int inRow0 = outJ * binsize;
        for (outI=0; outI<outNx; outI++)
        {
          outRowN[outI] = 0;
          outRow[outI] = 0.0;
        }
        for (inJ = inRow0; inJ < inRow0+binsize; inJ++)
          {
           inp = inData + inJ*inNx;
           for (outI=0; outI<outNx; outI++)
             for (i=0; i<binsize; i++, inp++)
               if (*inp != DRMS_MISSING_FLOAT)
                {
                  outRow[outI] += *inp;
                  outRowN[outI]++;
                }
         }
        for (outI=0; outI < outNx; outI++)
            outData[outJ*outNx + outI] = (outRowN[outI] ? 
                 (outRow[outI]/outRowN[outI]) : DRMS_MISSING_FLOAT);
      }   // end outJ
  }

// ---------------------------------------------------------------------

unsigned short *init_color_table(char * palette, int colors, int bytepercolor)
  {
  int c, i;
  int tablesize = (bytepercolor == 1 ? 256 : 65536);
  int maxcolor = tablesize - 2; // i.e. for 8-bit, max=254, leave top for missing
  int lastcolor = tablesize - 1;
  int mincolor = 1;
  unsigned short *table = (unsigned short *)malloc(colors*tablesize*sizeof(unsigned short));

  if (!table) return NULL;
  if (strcmp(palette,"GREY") == 0)
    for (i=0; i<tablesize; i++)
      for (c = 0; c<colors; c++)
        table[i*colors + c] = i;

  else if (colors == 3 && access(palette, R_OK) == 0)
    {
    FILE *ct = fopen(palette,"r");
    void proc_color(FILE *ct, unsigned short *t, int max);
    char buf[256];
    char *dot = rindex(palette, '.');
    if (dot && strcmp(dot, ".lut")==0)
        {
        float fr, fg, fb;
        table[0] = table[1] = table[2] = 0;
        table[colors*lastcolor] = table[colors*lastcolor+1] = table[colors*lastcolor+2] = lastcolor;
        for (i=1; i<lastcolor; i++)
	    {
	    fscanf(ct, "%f%f%f", &fr, &fg, &fb);
            c = 0;
            table[i*colors + c++] = maxcolor*fr;
            table[i*colors + c++] = maxcolor*fg;
            table[i*colors + c++] = maxcolor*fb;
	    }
        }
    else
        {
        // assume .sao type color table
        unsigned short r[tablesize], g[tablesize], b[tablesize];
        for (i=0; i<=tablesize; i++)
          r[i] = g[i] = b[i] = 0;
        while (fgets(buf,256,ct))
	  {
	  if (strncmp(buf,"RED:",4)==0)
	   proc_color(ct,r,maxcolor);
	  if (strncmp(buf,"GREEN:",6)==0)
	    proc_color(ct,g,maxcolor);
	  if (strncmp(buf,"BLUE:",5)==0)
	   proc_color(ct,b,maxcolor);
	  }
        table[0] = table[1] = table[2] = 0;
        table[colors*lastcolor] = table[colors*lastcolor+1] = table[colors*lastcolor+2] = lastcolor;
        for (i=1; i<lastcolor; i++)
          {
          c = 0;
          table[i*colors + c++] = r[i];
          table[i*colors + c++] = g[i];
          table[i*colors + c++] = b[i];
          }
        }
    fclose(ct);
    }
  return(table);
  }

void proc_color(FILE *ct, unsigned short *t, int max)
  {
  float cut, val;
  int pos = -1;
  double prev = 0.0;
  int prevloc = 0;
  double slope;
  pos = -1;
  prev = 0.0;
  prevloc = 0;
  while (fscanf(ct,"(%f,%f)", &cut, &val) == 2)
        {
        int loc = max*cut + 0.5;
        val *= max;
 // fprintf(stderr,"  got loc=%d, cut=%f %f\n",loc,cut,val);
        if (pos < 0)
                {
                pos = 0;
                prev = val;
                prevloc = 0;
                }
        if (loc == prevloc)
                slope = 0.0;
        else
                slope = (val - prev)/(loc - prevloc);
        while (pos <= loc)
                {
                t[pos] = slope * (pos-prevloc) + prev;
                pos++;
                }
        prev = t[loc];
        prevloc = loc;
        }
  while (pos <= max)
        t[pos++] = prev;
  }


char *set_scaling(DRMS_Array_t *in, double *minp, double *maxp, 
    int *nmissp, char *palette, int colors, int bytepercolor, int scaling, int usewhite)
  {
    int idata, ndata, nmiss = 0;
    float *inData = in->data;
    ndata = in->axis[0] * in->axis[1];
    int tablesize = (bytepercolor == 1 ? 256 : 65536);
    int color;
    int grey = strcmp(palette, "GREY") == 0;
    static unsigned short *table = NULL;
    static unsigned short *greytable = NULL;
    static unsigned short *colortable = NULL;
    char *out;
    out = (char *)malloc(ndata * colors * bytepercolor * sizeof(char));
    int maxcolor = (bytepercolor == 1 ? 254 : 65534);
    int missingcolor = (usewhite ? maxcolor + 1 : 0);

fprintf(stderr,"set_scaling called, scaling=%d, palette=%s, min=%f, max=%f",scaling,palette,*minp,*maxp);
fprintf(stderr,"\n   missingcolor %d maxcolor %d \n", missingcolor, maxcolor);
    if (!out)
        return(NULL);

    // get color table unless palette is GREY and scaling is MINMAX, do grey inline for that case only
    // since that combination is used to call set_scaling recursively for histogram generation
    if (grey && scaling == MINMAX && colors == 1)
      {
      if (!greytable)
        {
        int i;
        greytable = (unsigned short *)malloc(65536 * sizeof(unsigned short));
        for (i=0; i<65536; i++)
          greytable[i] = i;
        colortable = greytable;
        }
      }
    else
      {
      if (!table)
        table = init_color_table(palette, colors, bytepercolor);
      colortable = table;
      }

    switch(scaling)
      {
      case MAG:
      case MAGMINMAXGIVEN:
          {
          float bias = params_get_double(&cmdparams, "bias");
          for (idata=0; idata<ndata; idata++)
             {
             float val = inData[idata];
             if (!isnan(val))
                // inData[idata] = val/sqrt(fabs(val)+bias);
                inData[idata] = 400.0*val/pow((fabs(val)+bias), 0.75);
             }
          if (isnan(*minp)) *minp = -1500.0;
          if (isnan(*maxp)) *maxp = -*minp;
          // no break;
          }
      case MINMAX:
      case MAXMIN:
      case MINMAXGIVEN:
      case LOG:
      case SQRT:
          {
           double min;
           double max;
           if (bytepercolor != 1 && bytepercolor != 2)
                {fprintf(stderr,"wrong bytes per color\n"); return(NULL);}
           if (scaling == MINMAXGIVEN || scaling == MAG || scaling == MAGMINMAXGIVEN)
             {
             min = *minp;
             max = *maxp;
	     /*    fprintf(stderr,"debug pt 1: min=%f max=%f\n scaling=%d",min,max,scaling );  */
             }
           else
             {
             if (isnan(*minp) || isnan(*maxp))
               {
               int n=0;
               min = 1.0e20;
               max = -min;
	       /*     fprintf(stderr,"debug pt 2: min=%f max=%f\n",min,max); */
               for (idata=0; idata<ndata; idata++)
                  {
                   float val = inData[idata];
                   if (!isnan(val))
                     {
                     n++;
                     if (val > max) max = val;
                     if (val < min) min = val;
                     }
                 }
               if (n==0)
                 {
                 fprintf(stderr,"set_scaling, min or max is missing and no data values found.\n");
                 return(NULL);
                 }
               }
             if (isnan(*minp)) *minp = min;
             if (isnan(*maxp)) *maxp = max;
             min = *minp;
             max = *maxp;
             }
           double use_min, use_max;
           if (scaling == LOG)
             {
             if (min <= 0)
               min = 1.0e-10;
             if (max <= 0)
               max = 1.0e-10;
             use_min = log10(min);
             use_max = log10(max);
             }
           else if (scaling == SQRT)
             {
             if (min < 0)
               min = 0;
             if (max < 0)
               max = 0;
             use_min = sqrt(min);
             use_max = sqrt(max);
             }
           for (idata=0; idata<ndata; idata++)
             {
	       /*    fprintf(stderr,"debug pt 4: min=%f max=%f\n",min,max);   */
             float val = inData[idata];
             int newval;
             if (!isnan(val))
                {
                if (val <= min) val = min;
                if (val >= max) val = max;
                if (scaling == LOG)
		  newval = (maxcolor-1)*((log10(val)-use_min)/(use_max-use_min)) + 0.5;
                else if (scaling == SQRT )
                    newval = (maxcolor-1)*((sqrt(val)-use_min)/(use_max-use_min)) + 0.5;                         
                else
                  {
                  newval = (maxcolor-1)*(val-min)/(max-min) + 0.5;                   
                  if (scaling == MAXMIN)
                     newval = maxcolor - newval;
                  }
                if (newval >= maxcolor) newval = maxcolor;
                if (newval < 1) newval = 1;
                }
             else
                {
                newval = missingcolor;
                nmiss += 1;
                }
             for (color=0; color<colors; color++)
                {
                if (bytepercolor == 1)
                   *((unsigned char *)out + colors*idata + color) = colortable[newval*colors + color];
                else
                   *((unsigned short *)out + colors*idata + color) = colortable[newval*colors + color];
                }
             }
        *maxp = max;
        *minp = min;
        *nmissp = nmiss;
fprintf(stderr,"   min %f max %f nmissp %d ndata %d\n", *minp, *maxp, nmiss, ndata);
        break;
        }

    case HISTEQ:
      {
        int ihist, nhist=65536;
	int *hist = (int *)malloc(nhist * sizeof(int));
	int icolor, ncolor = (bytepercolor == 1 ? 255 : 65535);
        int total, want;
        unsigned int newval;
        double valspercolor;
        char *outa;
        outa = set_scaling(in, minp, maxp, nmissp, "GREY", 1, 2, MINMAX, 0);
        valspercolor = (ndata - *nmissp)/(double)ncolor;
        for (ihist = 0; ihist < nhist; ihist++)
          hist[ihist] = 0;
        for (idata=0; idata<ndata; idata++)
          hist[*((unsigned short *)outa+idata)] += 1;
        for (icolor=0, total=0, ihist=0; icolor<=ncolor; icolor++)
          {
          want = icolor * valspercolor;
          for ( ; ihist<nhist; ihist++)
            {
            if (total >= want)
                break;
            total += hist[ihist];
            hist[ihist] = icolor;
            }
          }
        for ( ; ihist<nhist; ihist++)
            hist[ihist] = ncolor;
        for (idata=0; idata<ndata; idata++)
          {
          if (!isnan(inData[idata]))
            {
            newval = hist[*((unsigned short *)outa+idata)];
            if (newval > maxcolor) newval = maxcolor;
            if (newval < 1) newval = 1;
            }
          else
            newval = missingcolor;
          for (color=0; color<colors; color++)
            {
            if (bytepercolor == 1)
               *((unsigned char *)out + colors*idata + color) = colortable[newval*colors + color];
            else
               *((unsigned short *)out + colors*idata + color) = colortable[newval*colors + color];
            }
          }
// fprintf(stderr," HISTEQ scaling %d \n", scaling);
        free(outa);
        free(hist);
        break;
      }

    case MINMAX99:
    case MAXMIN99:
      {
fprintf(stderr,"DOING set_scaling case MAXMIN99 \n");
          int ihist, nhist=65536;
          int *hist = (int *)malloc(nhist * sizeof(int));
          int icolor, ncolor = (bytepercolor == 1 ? 255 : 65535);
          int total, want;
          unsigned int newval;
          double newmin, newmax;
          char *outa;
fprintf(stderr," CALLING set_scaling case MINMAX %d \n", MINMAX);
	*minp = 0;
	*maxp = 0;
        outa = set_scaling(in, minp, maxp, nmissp, "GREY", 1, 2, MINMAX, 0);

fprintf(stderr,"CONTINUING: set_scaling, case MAXMIN99  %f %f \n", *minp, *maxp);

          for (ihist=0; ihist<nhist; ihist++)
             hist[ihist] = 0;
          for (idata=0; idata<ndata; idata++)
            {
              hist[*((unsigned short *)outa+idata)] += 1;
            }

        want = ndata/600;
fprintf(stderr,"ndata %d want %d nhist %d \n", ndata, want, nhist);
        total = 0;
        for (ihist=0; ihist<nhist; total += hist[ihist++])
          if (total >= want)
            break;
        newmin = *minp + ((double)ihist/(double)nhist) * (*maxp - *minp);

        total = 0;
        for (ihist=nhist - 1; ihist>=0; total += hist[ihist--])
         if (total >= want)
            break;
        newmax = *minp + ((double)ihist/(double)nhist) * (*maxp - *minp);
fprintf(stderr,"ihist %d newmax=%f newmin=%f \n", ihist, newmax, newmin);
        if (newmax <= newmin)
          {
          newmax = *maxp;
          newmin = *minp;
          }

fprintf(stderr,"newmax %f minp %f \n", newmax, *minp);
        for (idata=0; idata<ndata; idata++)
          {
          float val = inData[idata];
          unsigned short newval;
          if (!isnan(val)){
            if (val <= newmin) val=newmin;
            if (val >= newmax) val=newmax;
            newval = maxcolor*(val-newmin)/(newmax-newmin);
            if (scaling == MAXMIN99)
                newval = maxcolor - newval;
            if (newval > maxcolor) newval = maxcolor;
            if (newval < 1) newval = 1;
            }
            
          else
            newval = missingcolor;
          for (color=0; color<colors; color++)
            {
            if (bytepercolor == 1)
               *((unsigned char *)out + colors*idata + color) = colortable[newval*colors + color];
            else
               *((unsigned short *)out + colors*idata + color) = colortable[newval*colors + color];
            }
          }
        free(outa);
        free(hist);
        break;
      }
    }
  return(out);
  }

// ----------------------------------------------------------------------
int make_png(char *filename, unsigned char *data, 
           int height, int width, char *palette, int bytepercolor, int colors)
{ 
  int row;
  FILE *fp; // = NULL;
  png_byte **row_pointers = NULL;
  png_structp png_ptr = NULL;
  png_infop info_ptr = NULL;

  fp = fopen(filename, "w");
  if (!fp)
    PNGDIE("cant open output file");

  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, 
          NULL, NULL);
  if (!png_ptr)
    PNGDIE("cant create png_struct");

  info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr)
    PNGDIE("cant create png_infop");

  if (setjmp(png_jmpbuf(png_ptr)))
    PNGDIE("png fatal error");

  png_init_io(png_ptr, fp);
// fprintf(stderr, "ready png, bitdepth=8*bytepercolor=%d, colors=%d \n",8*bytepercolor,colors);
  png_set_IHDR(png_ptr, info_ptr, width, height, 8*bytepercolor,
        (colors == 1 ? PNG_COLOR_TYPE_GRAY : PNG_COLOR_TYPE_RGB), PNG_INTERLACE_ADAM7,
        PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
  png_write_info(png_ptr, info_ptr);
  if (bytepercolor > 1)
    png_set_swap(png_ptr);
  png_set_gAMA(png_ptr, info_ptr, 2.2);
  png_set_sRGB(png_ptr, info_ptr, PNG_sRGB_INTENT_ABSOLUTE);
  row_pointers = (png_byte **)malloc(height * sizeof(png_byte *));

  for (row=0; row<height ; row++)
    {
    if (bytepercolor == 1)
      row_pointers[height - row - 1] = (png_byte *)((char *)data + colors*row*width);
    else
      row_pointers[height - row - 1] = (png_byte *)((short *)data + colors*row*width);
    }
  png_set_rows(png_ptr, info_ptr, row_pointers);
  png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_SWAP_ENDIAN, NULL);
  fclose(fp);

  png_destroy_write_struct(&png_ptr, &info_ptr);
  free(row_pointers);
  return(0);
  }

/* center whith whole pixel shifts and rotate by 180 if needed */
int upNcenter(DRMS_Array_t *arr, ObsInfo_t *ObsLoc)
  {
  int nx, ny, ix, iy, i, j, xoff, yoff;
  double rot, x0, y0, mid;
  float *data;
  if (!arr || !ObsLoc)
    return(1);
  data = arr->data;
  nx = arr->axis[0];
  ny = arr->axis[1];
  x0 = ObsLoc->crpix1 - 1;
  y0 = ObsLoc->crpix2 - 1;
  mid = (nx-1.0)/2.0;
  if ((rot = fabs(ObsLoc->crota2)) > 179 && rot < 181)
    {
    // rotate image by 180 degrees by a flip flip
    float val;
    int half = nx / 2;
    int odd = nx & 1;
    if (odd) half++;
    for (iy=0; iy<half; iy++)
      {
      for (ix=0; ix<nx; ix++)
        {
        i = iy*nx + ix;
        j = (ny - 1 - iy)*nx + (nx - 1 - ix);
        val = data[i];
        data[i] = data[j];
        data[j] = val;
        }
      }
    x0 = nx - x0;
    y0 = ny - y0;
    rot = ObsLoc->crota2 - 180.0;
    if (rot < -90.0) rot += 360.0;
    ObsLoc->crota2 = rot;
    }
  xoff = round(x0 - mid);
  yoff = round(y0 - mid);
  if (abs(xoff) > 1.0)
    {
    for (iy=0; iy<ny; iy++)
      {
      float valarr[nx];
      for (ix=0; ix<nx; ix++)
        {
        int jx = ix - xoff;
        if (jx >= nx) jx -= nx;
        if (jx < 0) jx += nx;
        valarr[jx] = data[iy*nx + ix];
        }
      for (ix=0; ix<nx; ix++)
        data[iy*nx + ix] = valarr[ix];
      }
    x0 -= xoff;
    }
  if (abs(yoff) > 1.0)
    {
    for (ix=0; ix<nx; ix++)
      {
      float valarr[ny];
      for (iy=0; iy<ny; iy++)
        {
        int jy = iy - yoff;
        if (jy >= ny) jy -= ny;
        if (jy < 0) jy += ny;
        valarr[jy] = data[iy*nx + ix];
        }
      for (iy=0; iy<ny; iy++)
        data[iy*nx + ix] = valarr[iy];
      }
    y0 -= yoff;
    }
  ObsLoc->crpix1 = x0 + 1;
  ObsLoc->crpix2 = y0 + 1;
  return(0);
  }

int crop_image(DRMS_Array_t *arr, ObsInfo_t *ObsLoc)
  {
  int nx, ny, ix, iy, i, j, xoff, yoff;
  double x0, y0;
  double rsun = ObsLoc->rsun_obs/ObsLoc->cdelt1;
  double scale, crop_limit;
  float *data;
  if (!arr || !ObsLoc)
    return(1);
  data = arr->data;
  nx = arr->axis[0];
  ny = arr->axis[1];
  x0 = ObsLoc->crpix1 - 1;
  y0 = ObsLoc->crpix2 - 1;
  scale = 1.0/rsun;
  crop_limit = 0.99975; // 1 - 1/4000, 1/2 HMI pixel.
  for (iy=0; iy<ny; iy++)
    for (ix=0; ix<nx; ix++)
      {
      double x, y, R2;
      float *Ip = data + iy*nx + ix;
      if (drms_ismissing_float(*Ip))
        continue;
      x = ((double)ix - x0) * scale; /* x,y in pixel coords */
      y = ((double)iy - y0) * scale;
      R2 = x*x + y*y;
      if (R2 > crop_limit)
        *Ip = DRMS_MISSING_FLOAT;
      }
  return(0);
  }

#define CHECK(keyname) {if (status) {fprintf(stderr,"Keyword failure to find: %s, status=%d\n",keyname,status); *rstatus=status; return(NULL);}}

ObsInfo_t *GetObsInfo(DRMS_Segment_t *seg, ObsInfo_t *pObsLoc, int *rstatus)
  {
  TIME t_prev;
  DRMS_Record_t *rec;
  TIME t_obs;
  double dv;
  ObsInfo_t *ObsLoc;
  int status;

  if (!seg || !(rec = seg->record))
    { *rstatus = 1; return(NULL); }

  ObsLoc = (pObsLoc ? pObsLoc : (ObsInfo_t *)malloc(sizeof(ObsInfo_t)));
  if (!pObsLoc)
    memset(ObsLoc, 0, sizeof(ObsInfo_t));

  t_prev = ObsLoc->t_obs;
  t_obs = drms_getkey_time(rec, "T_OBS", &status); CHECK("T_OBS");

  if (t_obs <= 0.0)
    { *rstatus = 2; return(NULL); }

  if (t_obs != t_prev)
    {
    ObsLoc->crpix1 = drms_getkey_double(rec, "CRPIX1", &status); CHECK("CRPIX1");
    ObsLoc->crpix2 = drms_getkey_double(rec, "CRPIX2", &status); CHECK("CRPIX2");
    ObsLoc->crval1 = drms_getkey_double(rec, "CRVAL1", &status); CHECK("CRVAL1");
    ObsLoc->crval2 = drms_getkey_double(rec, "CRVAL2", &status); CHECK("CRVAL2");
    ObsLoc->cdelt1 = drms_getkey_double(rec, "CDELT1", &status); CHECK("CDELT1");
    ObsLoc->cdelt2 = drms_getkey_double(rec, "CDELT2", &status); CHECK("CDELT1");
    ObsLoc->crota2 = drms_getkey_double(rec, "CROTA2", &status); CHECK("CROTA2");
    ObsLoc->sina = sin(ObsLoc->crota2*Deg2Rad);
    ObsLoc->cosa = sqrt (1.0 - ObsLoc->sina*ObsLoc->sina);
    ObsLoc->rsun_obs = drms_getkey_double(rec, "RSUN_OBS", &status);
    if (status)
      {
      double dsun_obs = drms_getkey_double(rec, "DSUN_OBS", &status); CHECK("DSUN_OBS");
      ObsLoc->rsun_obs = asin(696000000.0/dsun_obs)/arcsec2Rad;
      }
    ObsLoc->obs_vr = drms_getkey_double(rec, "OBS_VR", &status); CHECK("OBS_VR");
    ObsLoc->obs_vw = drms_getkey_double(rec, "OBS_VW", &status); CHECK("OBS_VW");
    ObsLoc->obs_vn = drms_getkey_double(rec, "OBS_VN", &status); CHECK("OBS_VN");
    ObsLoc->obs_b0 = drms_getkey_double(rec, "CRLT_OBS", &status); CHECK("CRLT_OBS");
    ObsLoc->t_obs = t_obs;
    }
  *rstatus = 0;
  return(ObsLoc);
  }


// ----------------------------------------------------------------------
