/*
 * Name:		obstime2maskid.c
 *
 * Description:		Convert TAI to noise mask id.
 *
 * Version:		v1.0		Dec 11 2013
 *				v2.0		Apr 22 2015
 *
 * Note:
 * v2.0
 * 		Move the times into external text file, return -1 if not found
 *		times in the text file must be sorted
 *			
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define TIMEFILE	"proj/mag/ambig/apps/data/masktime.txt"

int obstime2maskid(TIME tobs)
{
    char ttemp[64];
    int MaskIndex = 0;
	FILE *fptr;
	TIME tid;
	
	/* Snippet from Art to obtain data file */
	
	char lpath[256], spath[256], *needle;

	/* Obtain data path relative to binary path. */
	if (readlink("/proc/self/exe", spath, sizeof(spath)) == -1)
	{
		fprintf(stderr, "Cannot locate this binary.\n");
		return -1;
	}
	else
	{
		if ((needle = strstr(spath, DRMS_ARCH)) != NULL)
		{
			needle += strlen(DRMS_ARCH);
			*needle = '\0';
			
			/* We have the path to the make output directory. Now we need
			 * to get the parent and append the data file */
			snprintf(lpath,
					sizeof(lpath),
					"%s/../%s",
					spath,
					TIMEFILE);
		}
		else
		{
			fprintf(stderr, "Cannot find architecture %s subpath.\n", DRMS_ARCH);
			return -1;
		}
	}
	printf("lpath: %s\n", lpath);

	/* Start */
	
	fptr = fopen(lpath, "r");
	if (!fptr) {
		printf("Time file for mask cannot be found, return -1.\n");
		return -1;		// if error, return -1, so no mask will be found
	}
	
	while (fgets(ttemp, 64, fptr) != NULL) {
		tid = sscan_time(ttemp);
//		printf("id=%d; t=%s; tid=%lf\n", MaskIndex, ttemp, tid);
		if (tobs < tid) {
			fclose(fptr);
			return MaskIndex;
		}
		MaskIndex++;
	}

    fclose(fptr);

	return MaskIndex;
	
}

