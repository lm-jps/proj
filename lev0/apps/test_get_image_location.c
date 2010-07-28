#ident "$Header: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/lev0/apps/test_get_image_location.c,v 1.1 2010/07/28 22:14:04 carl Exp $"
/*############################################################################
# Name:        test_get_image_location.c                                     #
# Description:Used to unit test get_image_location.c                         # 
# Author:      Carl                                                          #
# Date:        Move from EGSE to JSOC software environment on July 28, 2010 #
############################################################################*/

/******************** includes ******************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdarg.h>
#include "drms.h"
#include "drms_names.h"
#include "jsoc_main.h"
#include <printk.h>
#include "get_image_location.c"

/************* modules definitions **************************************/
ModuleArgs_t module_args[] =
{
  {ARG_STRING, "src", "Not Specified", "day file source"},
  {ARG_STRING, "in", "Not Specified", "full path to day file"},
  {ARG_STRING, "out", "Not Specified", "Series name"},
  {ARG_FLAG, "p", "0", "print values of keywords to standard out"},
  {ARG_END}
};
ModuleArgs_t   *ggModArgs=module_args;
char* module_name = "test_get_image_location";


/*************************************************************************
 * DoIT                                                                  *
 * FUNCTION: DoIt(void                                                   *
 * DESCRIPTION: Use function to create jsoc module                       *
 *************************************************************************/
int DoIt(void)
{
  /* variables */
  int i=1;
  int ncnt=1;
  int rstatus;
  Image_Location *p_imageloc;
  Image_Location imageloc[10];
  /* set environment variables */
  printf("test_get_image_location() ran!!\n");

  // Obs time TIME tobs;

  // For AIA:1,2,3,4. For HMI:1 or 2 int camera;

  // For AIA:SDO/AIA. For HMI:SDO/HMI char telescope[GMP_MAX_TELESCOPE_STR];

  //wavelength int wavelength;

  //Now fill in info for call to Carl's get_image_location()
  for(i=0; i < ncnt; i++) 
  {
      //imageloc[i].tobs =1648159234; //"2010-03-24T22:00:00Z";
      //time_convert time=2010-04-09T20:00:00Z o=jsoc ->1049918434
      imageloc[i].tobs =1049918434; 
      imageloc[i].camera = 1;
      imageloc[i].wavelength=0;
      strcpy(imageloc[i].telescope,"SDO/HMI");
  }
  p_imageloc = imageloc;
  rstatus = get_image_location(drms_env, ncnt, &p_imageloc);

  printf("rstatus is %d\n",  rstatus);
  printf("x is %f\n",  p_imageloc->x);
  printf("y is %f\n", p_imageloc->y);
  printf("instrot is %f\n", p_imageloc->instrot);
  printf("imscale is %f\n", p_imageloc->imscale);
  printf("yinrtb is %f\n", p_imageloc->yinrtb);
  printf("zinrtb is %f\n", p_imageloc->zinrtb);

  return(0);

}


