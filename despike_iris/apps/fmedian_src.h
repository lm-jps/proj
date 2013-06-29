#include "idl_export.h"

#define SWAP(a,b) (atemp=a,a=b,b=atemp)
#define SWAP2(TYPE,a,b) {TYPE tmp; tmp=a; a=b; b=tmp;}


#define MAX2(a,b) ( (a) > (b) ? (a) : (b) )
#define MIN2(a,b) ( (a) < (b) ? (a) : (b) )

#define QSRT_k_func2_sig(TYPE)\
  void qsrt_k_##TYPE(TYPE a[], long k, long n)

#define QSRT_k_call2(TYPE,a,k,n) qsrt_k_##TYPE(a,k,n)

QSRT_k_func2_sig(UCHAR);
QSRT_k_func2_sig(IDL_INT);
QSRT_k_func2_sig(IDL_UINT);
QSRT_k_func2_sig(IDL_LONG);
QSRT_k_func2_sig(IDL_ULONG);
QSRT_k_func2_sig(IDL_LONG64);
QSRT_k_func2_sig(IDL_ULONG64);
QSRT_k_func2_sig(float);
QSRT_k_func2_sig(double);


#define FMEDIAN_CORE_FUNC2_SIG(TYPE) \
  void fmedian_core_##TYPE						\
  (TYPE *in,TYPE *out,TYPE *work,					\
   int nx, int ny, int xbox,int ybox,					\
   int always,TYPE missing,int find_missing)

#define FMEDIAN_CALL(TYPE,in,out,work,nx,ny,xbox,ybox,			\
		     always, missing, find_missing)			\
  fmedian_core_##TYPE((TYPE*)in,(TYPE*)out,(TYPE*)work,nx,ny,		\
		      xbox,ybox, always, missing, find_missing)

/* Declare signatures for them all */
FMEDIAN_CORE_FUNC2_SIG(UCHAR);
FMEDIAN_CORE_FUNC2_SIG(IDL_INT);
FMEDIAN_CORE_FUNC2_SIG(IDL_UINT);
FMEDIAN_CORE_FUNC2_SIG(IDL_LONG);
FMEDIAN_CORE_FUNC2_SIG(IDL_ULONG);
FMEDIAN_CORE_FUNC2_SIG(IDL_LONG64);
FMEDIAN_CORE_FUNC2_SIG(IDL_ULONG64);
FMEDIAN_CORE_FUNC2_SIG(float);
FMEDIAN_CORE_FUNC2_SIG(double);

