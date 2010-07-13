/*
 * ARGUMENT CHECKING
 * (mexargcheck.c)
 * 
 * Michael Turmon, 2002
 * 
 */

#ifndef _mexargcheck_h_
#define _mexargcheck_h_

#ifdef __cplusplus
extern "C" {
#endif

/* simple size/type checking */
#define IsFullRealArray(pm) (mxIsDouble(pm) && \
                             !mxIsComplex(pm) && \
                             !mxIsSparse(pm))
#define IsFullRealMatrix(pm) (IsFullRealArray(pm) && \
                             (mxGetNumberOfDimensions(pm) <= 2))
#define IsRealVector(pm) (IsFullRealMatrix(pm) && \
                          (mxGetM(pm) == 1 || \
                           mxGetN(pm) == 1))
/* Vectors in broad sense */
#define IsRealLooseVector(pm) (IsFullRealMatrix(pm) && \
                          (mxGetM(pm) <= 1 || \
                           mxGetN(pm) <= 1))
#define IsLooseVector(pm) (mxGetM(pm) <= 1 || \
                           mxGetN(pm) <= 1)
#define IsRealScalar(pm) (IsFullRealMatrix(pm) && \
                           mxGetNumberOfElements(pm) == 1)
/* Scalar in broad sense */
#define IsRealLooseScalar(pm) (IsFullRealMatrix(pm) && \
                           mxGetNumberOfElements(pm) <= 1)
#define IsLooseScalar(pm) (mxGetNumberOfElements(pm) <= 1)
#define IsEmpty(pm)       (mxGetNumberOfElements(pm) == 0)

/*
 *  Function templates 
 */

int
mexargparse(int narg, 
	    const mxArray **args,
	    const char **names, const char **specs, const char **msgs, const char *gen_msg);
int start_sizechecking();
int sizeinit(const mxArray *pm);
int sizecheck(const char *msg, int num);
int sizecheck_msg(const char *msg, const char **argnames, int num);
int sizeagree(const mxArray *pm);
int sizeagreeM(const mxArray *pm);
int sizeagreeN(const mxArray *pm);
int sizeagreeMN(const mxArray *pm);
int sizesare(int *s, int d);
int sizeis(int m, int n);
int sizeis3(int m, int n, int p);
int sizeisM(int m);
int sizeisN(int n);
int sizeisMN(int m, int n);

#ifdef __cplusplus
}	/* extern "C" */
#endif

#endif /* _mexargcheck_h_ */
