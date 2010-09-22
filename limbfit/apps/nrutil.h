#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

int sort(unsigned long n, float *arr);
int indexx(unsigned long n, float *arr, unsigned long *indx);
float *vector(long nl, long nh, int *status);
int *ivector(long nl, long nh, int *status);
unsigned long *lvector(long nl, long nh, int *status);
void free_vector(float *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_lvector(unsigned long *v, long nl, long nh);
#endif /* _NR_UTILS_H_ */
