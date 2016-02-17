/*
 * Name:		obstime2maskid.c
 *
 * Description:		Convert TAI to noise mask id.
 *
 * Version:		v1.0		Dec 11 2013
 *
 * Issues:
 *			
 */

#include "drms.h"

int obstime2maskid(tobs)
    TIME tobs;
{
    char ttemp[64];
    TIME ChangeTime1, ChangeTime2, ChangeTime3, ChangeTime4, ChangeTime5, ChangeTime6, ChangeTime7, ChangeTime8, ChangeTime9;
    int MaskIndex = 0;

    strcpy(ttemp, "2010.12.13_19:47:00_TAI");
    ChangeTime1 = sscan_time(ttemp);

    strcpy(ttemp, "2011.07.13_18:35:00_TAI");
    ChangeTime2 = sscan_time(ttemp);

    strcpy(ttemp, "2012.01.18_18:15:00_TAI");
    ChangeTime3 = sscan_time(ttemp);

    strcpy(ttemp, "2013.03.14_06:40:00_TAI");
    ChangeTime4 = sscan_time(ttemp);

    strcpy(ttemp, "2014.01.15_19:18:00_TAI");
    ChangeTime5 = sscan_time(ttemp);

    if (tobs<ChangeTime1) MaskIndex = 0;
    if (tobs>ChangeTime1 && tobs<ChangeTime2) MaskIndex = 1;
    if (tobs>ChangeTime2 && tobs<ChangeTime3) MaskIndex = 2;
    if (tobs>ChangeTime3 && tobs<ChangeTime4) MaskIndex = 3;
    if (tobs>ChangeTime4 && tobs<ChangeTime5) MaskIndex = 4;
    if (tobs>ChangeTime5) MaskIndex = 5;

  return MaskIndex;
}

