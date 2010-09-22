#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include "nrutil.h"
#define NR_END 1
#define FREE_ARG char*

float *vector(long nl, long nh, int *status)
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;
	*status=0;
	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) {
		//nrerror("allocation failure in vector()");
		*status=-1;
	}
	return v-nl+NR_END;
}

int *ivector(long nl, long nh, int *status)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;
	*status=0;
	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) {
		//nrerror("allocation failure in ivector()");
		*status=-1;
	}
	return v-nl+NR_END;
}


unsigned long *lvector(long nl, long nh, int *status)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
	unsigned long *v;
	*status=0;
	v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
	if (!v){
		//nrerror("allocation failure in lvector()");
		*status=-1;
	}
	return v-nl+NR_END;
}

void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_lvector(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}



