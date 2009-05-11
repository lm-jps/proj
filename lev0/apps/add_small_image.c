#include "mypng.h"

#define PNGDIE(msg) {fprintf(stderr,"%s\n",msg); \
	if (png_ptr) png_destroy_write_struct(&png_ptr, &info_ptr); \
	if (row_pointers) free(row_pointers); \
	if (fp) fclose(fp); \
	return(1);}

#define GREY	1
#define MINMAX	1
#define MAXMIN	2
#define HISTEQ  3
#define MINMAX99  4
#define MAXMIN99  5

int make_png(char *filename, unsigned char *data, int height, int width, int pallette, int bytepercolor, int colors)
  { // make png
  int row;
  FILE *fp = NULL;
  png_byte **row_pointers = NULL;
  png_structp png_ptr = NULL;
  png_infop info_ptr = NULL;

  fp = fopen(filename, "wb");
  if (!fp)
    PNGDIE("cant open output file");
  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (!png_ptr)
    PNGDIE("cant create png_struct");
  info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr)
    PNGDIE("cant create png_infop");
  if (setjmp(png_jmpbuf(png_ptr)))
    PNGDIE("png fatal error");
  png_init_io(png_ptr, fp);
  png_set_IHDR(png_ptr, info_ptr, width, height, 8*bytepercolor,
	PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE,
	PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
  png_write_info(png_ptr, info_ptr);
  row_pointers = (png_byte **)malloc(height * sizeof(png_byte *));
  for (row=0; row<height ; row++)
    {
    if (bytepercolor == 1)
      row_pointers[height - row - 1] = (png_byte *)((char *)data + row*width);
    else
      row_pointers[height - row - 1] = (png_byte *)((short *)data + row*width);
    }
  png_set_rows(png_ptr, info_ptr, row_pointers);
  png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_SWAP_ENDIAN, NULL);
  fclose(fp);
  png_destroy_write_struct(&png_ptr, &info_ptr);
  free(row_pointers);
  return(0);
  }

char *set_scaling(DRMS_Array_t *in, double *minp, double *maxp, int *nmissp,
        int pallette, int colors, int bytepercolor, int scaling)
  {
  char *out;
  int idata, ndata, nmiss = 0;
  float *inData = in->data;
  ndata = in->axis[0] * in->axis[1];
  out = (char *)malloc(ndata * bytepercolor * sizeof(char));
  int missingcolor = (bytepercolor == 1 ? 255 : 65535);
  int maxcolor = missingcolor - 1;
  if (!out)
	return(NULL);

  switch(scaling)
    {
    case MINMAX:
    case MAXMIN:
	{
	double min = 1.0e20;
	double max = -min;
	if (bytepercolor != 1 && bytepercolor != 2)
		{fprintf(stderr,"wrong bytes per color\n"); return(NULL);}
	for (idata=0; idata<ndata; idata++)
	  {
	  float val = inData[idata];
	  if (!isnan(val))
	    {
            if (val > max) max = val;
            if (val < min) min = val;
            }
          }
        for (idata=0; idata<ndata; idata++)
          {
          float val = inData[idata];
	  unsigned short newval;
          if (!isnan(val))
	    {
            newval = maxcolor*(val-min)/(max-min);
	    if (scaling == MAXMIN)
		newval = maxcolor - newval;
	    }
          else
	    {
            newval = missingcolor;
	    nmiss += 1;
	    }
	  if (bytepercolor == 1)
	    *((unsigned char *)out + idata) = newval;
	  else
            *((unsigned short *)out + idata) = newval;
          }
	*maxp = max;
	*minp = min;
	*nmissp = nmiss;
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
	outa = set_scaling(in, minp, maxp, nmissp, GREY, 1, 2, MINMAX);
	valspercolor = (ndata - *nmissp)/(double)ncolor;
	for (ihist=0; ihist<nhist; ihist++)
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
          if (!isnan(inData[idata]))
	    {
	    if (bytepercolor == 1)
	      *((unsigned char *)out + idata) = hist[*((unsigned short *)outa+idata)];
	    else
              *((unsigned short *)out + idata) = hist[*((unsigned short *)outa+idata)];
	    }
          else
            {
	    if (bytepercolor == 1)
	      *((unsigned char *)out + idata) = missingcolor;
	    else
              *((unsigned short *)out + idata) = missingcolor;
            }
	free(outa);  
	free(hist);
	break;
	}
    case MINMAX99:
    case MAXMIN99:
	{
	int ihist, nhist=65536;
	int *hist = (int *)malloc(nhist * sizeof(int));
	int icolor, ncolor = (bytepercolor == 1 ? 255 : 65535);
	int total, want;
	unsigned int newval;
	double newmin, newmax;
	char *outa;
	outa = set_scaling(in, minp, maxp, nmissp, GREY, 1, 2, MINMAX);

        // create histogram
	for (ihist=0; ihist<nhist; ihist++)
	  hist[ihist] = 0;
	for (idata=0; idata<ndata; idata++)
	  hist[*((unsigned short *)outa+idata)] += 1;

        // find new min and max vals at cut points
        want = ndata/200;
        total = 0;
	for (ihist=0; ihist<nhist; total += hist[ihist++])
	  if (total >= want)
	    break;
        newmin = *minp + ((double)ihist/(double)nhist) * (*maxp - *minp);
        total = 0;
	for (ihist=nhist; ihist>=0; total += hist[ihist--])
	  if (total >= want)
	    break;
        newmax = *minp + ((double)ihist/(double)nhist) * (*maxp - *minp);
        if (newmax <= newmin)
          {
          newmax = *maxp;
          newmin = *minp;
          }

        for (idata=0; idata<ndata; idata++)
          {
          float val = inData[idata];
          unsigned short newval;
          if (!isnan(val))
            {
	    if (val <= newmin) val=newmin;
	    if (val >= newmax) val=newmax;
            newval = maxcolor*(val-newmin)/(newmax-newmin);
            if (scaling == MAXMIN99)
                newval = maxcolor - newval;
            }
          else
            newval = missingcolor;
          if (bytepercolor == 1)
            *((unsigned char *)out + idata) = newval;
          else
            *((unsigned short *)out + idata) = newval;
          }
	free(outa);  
	free(hist);
	break;
	}
    }
  return(out);
  }

// expects short in and returns float DRMS_Array_t's

rebinArraySF(DRMS_Array_t *out, DRMS_Array_t *in)
  {
  short *inData = in->data;
  short *inp;
  float *outData = out->data;
  float *outp;
  int inNx=in->axis[1], inNy=in->axis[0], outNx=out->axis[1], outNy=out->axis[0];
  int inI, inJ, outI, outJ, i;
  int binsize = inNx/outNx;
  double inRow[inNx], outRow[outNx];
  int inRowN[inNx], outRowN[outNx];
  for (outJ=0; outJ<outNy; outJ++)
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
	  if (*inp != DRMS_MISSING_SHORT)
	    {
            outRow[outI] += *inp;
	    outRowN[outI]++;
	    }
      }
    for (outI=0; outI<outNx; outI++)
      outData[outJ*outNx + outI] = (outRowN[outI] ? (outRow[outI]/outRowN[outI]) : DRMS_MISSING_FLOAT);
    }
  }

// Given rec and big, create smaller scale segment reduced by doScaleFits in each
// dimension and create a smaller scale image reduced by doScaleImage.
// Setting either doScale{Fits,Image} to 0 will prevent that component
// from being generated.  Typical use, to convert a 4096x4096 image to
// a 512x512 fits file and 256x256 pixel image file call with 
//  add_small_array(rs, array, 8, 16);
// Scaling for the image will be equal spaced by histogram
// Filenames will be put into segments number 1 and 2 respectively if
// the segments are in the .jsd.
// segment 1 should be a vardim FITS segment and segment 2 should
// be GENERIC.  The
// small FITS and image files will be placed in the record directory
// with filenames derived from segment number 0.
//
// return value is 0 if all OK, else error code.

// Now expects data in "big" to be shorts.

int add_small_array(DRMS_Record_t *rec, DRMS_Array_t *big, int doScaleFits, int doScaleImage)
  {
  DRMS_Segment_t *bigSeg = drms_segment_lookupnum(rec, 0);
  char recdir[DRMS_MAXPATHLEN];
  int bigNx, bigNy;
  int status;

  drms_record_directory(rec, recdir, 1);
  bigNy = big->axis[0];
  bigNx = big->axis[1];
  /* assume second segment, if it is present, is vardims FITS for small image */
  if (doScaleFits && rec->segments.num_total > 1)
    {
    DRMS_Segment_t *smallSeg;
    DRMS_Array_t *small;
    DRMS_SegmentDimInfo_t smallDimInfo;
    int smallNx, smallNy;
    int smallDims[2];
    char smallFileName[DRMS_MAXSEGFILENAME];

    smallNy = bigNy/doScaleFits;
    smallNx = bigNx/doScaleFits;
    smallDims[0] = smallNy;
    smallDims[1] = smallNx;
    small = drms_array_create(DRMS_TYPE_FLOAT, 2, smallDims, NULL, &status);

    rebinArraySF(small, big);
    strcpy(smallFileName, bigSeg->info->name);
    strcat(smallFileName, "_sm.fits");
    smallSeg = drms_segment_lookupnum(rec, 1);
    strcpy(smallSeg->filename, smallFileName);
    smallDimInfo.naxis = 2;
    smallDimInfo.axis[0] = smallDims[0];
    smallDimInfo.axis[1] = smallDims[1];
    drms_segment_setdims(smallSeg, &smallDimInfo);
    // For SDO CCD images of 14 bits use FITS shorts as unsigned and keep full range.
    small->bzero = -(4.0*8192.0);
    small->bscale = 4.0;
    small->israw = 0;
    drms_array_convert_inplace(DRMS_TYPE_SHORT, small->bzero, small->bscale, small);

    small->bzero = 8192.0;
    small->bscale = 0.25;
    small->israw = 1;
    drms_segment_writewithkeys(smallSeg, small, 0);

    drms_free_array(small);
    }

  /* assume third segment, if exists, is for small image, generic segment */
  if (doScaleImage && rec->segments.num_total > 2)
    {
    DRMS_Segment_t *imageSeg;
    DRMS_Array_t *imageArray = NULL;
    char *imgData = NULL;
    int imageDims[2];
    int imageNx, imageNy;
    char imageFileName[DRMS_MAXSEGFILENAME];
    char imagePathName[DRMS_MAXPATHLEN];
    int bytepercolor = 1;
    int nmissing;
    double min, max;
    imageNy = bigNy/doScaleImage;
    imageNx = bigNx/doScaleImage;
    imageDims[0] = imageNy;
    imageDims[1] = imageNx;
    imageArray = drms_array_create(DRMS_TYPE_FLOAT, 2, imageDims, NULL, &status);

    strcpy(imageFileName, bigSeg->info->name);
    strcat(imageFileName, ".png");
    strcpy(imagePathName, recdir);
    strcat(imagePathName, "/");
    strcat(imagePathName, imageFileName);
    imageSeg = drms_segment_lookupnum(rec, 2);
    strcpy(imageSeg->filename, imageFileName);
    rebinArraySF(imageArray, big);
    imgData = (unsigned char *)set_scaling(imageArray, &min, &max, &nmissing, GREY, 1, bytepercolor, MINMAX99);
    drms_free_array(imageArray);
    if (!imgData)
      {
      fprintf(stderr,"scaling failed, status=%d\n",status);
      return(status);
      }
    // should add text info to .png file with min, max, nmissing, etc.
    if (make_png(imagePathName, imgData, imageNy, imageNx, GREY, 1, bytepercolor) != 0)
      {
      fprintf(stderr,"png failed, status=%d\n",status);
      free(imgData);
      return(status);
      }
    free(imgData);
    }
  return(0);
  }

