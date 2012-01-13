/*
 * mex_stubs_unimp.c
 * 
 * These are stubs for *unimplemented* mex functions.
 * They are used in both mex_stubs.c and pymex_stubs.c, 
 * which #include this file.
 *
 * It's intended that these are:
 * (a) the exotic mx functions -- cell arrays and the like.
 * (b) the unuseful mx/mex functions -- mexCallMatlab etc.
 *
 * If you want to implement one of these functions in one or the other file,
 * move the code for "unimplemented" to the appropriate stubs file and 
 * put the real implementation in the other file.
 *
 * Michael Turmon, 2008
 */

/* no #includes needed */

/* remember, the internals of the mxArray structure will be different
 * here depending on compiler options in the parent file.
 */

/* 
 * un-useful mexFOO functions when outside the Matlab environment
 */

/* this is used up to r13 */
int
mexPutArray(mxArray *pa, const char *ws)
{
  return(0);
}

/* this is used in r13 and later */
int 
mexPutVariable(const char *workspace, const char *name, const mxArray *parray)
{
  return(0);
}

int
mexEvalString(const char *string)
{
   return(0);
}

int
mexCallMATLAB(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[], 
	      const char *name)
{
   return(0);
}

void mexLock(void) {return;}
void mexUnlock(void) {return;}
bool mexIsLocked(void) {return true;}

int  mexAtExit(void (*exit_fcn)(void)) {
  mex2c_unimplemented("mxAtExit");return(0);
}



/* 
 * creation of sparses
 */
mxArray *
mxCreateSparse(mwSize m, mwSize n, mwSize nzmax, mxComplexity flag) {
  mex2c_unimplemented("mxCreateSparse"); return 0;}
mxArray *
mxCreateSparseLogicalMatrix(mwSize m, mwSize n, mwSize nzmax) {
  mex2c_unimplemented("mxCreateSparseLogicalMatrix"); return 0;}

/* 
 * creation of structs/cells
 */
mxArray *
mxCreateCellMatrix(mwSize m, mwSize n) {
  mex2c_unimplemented("mxCreateCellMatrix"); return 0;}
mxArray *
mxCreateCellArray(mwSize ndim, const mwSize *dims) {
  mex2c_unimplemented("mxCreateCellArray"); return 0;}
mxArray *
mxCreateStructMatrix(mwSize m, mwSize n, int nf, const char **names) {
  mex2c_unimplemented("mxCreateStructMatrix"); return 0;}
mxArray *
mxCreateStructArray(mwSize ndim, const mwSize *dims, int nf, const char **names) {
  mex2c_unimplemented("mxCreateStructArray"); return 0;}

/* get/set sparse */
mwIndex *mxGetIr(const mxArray *pa) {
  mex2c_unimplemented("mxGetIr"); return 0;}
void     mxSetIr(mxArray *pa, mwIndex *newir) {
  mex2c_unimplemented("mxSetIr");}
mwIndex *mxGetJc(const mxArray *pa) {
  mex2c_unimplemented("mxGetJc"); return 0;}
void     mxSetJc(mxArray *pa, mwIndex *newjc) {
  mex2c_unimplemented("mxSetJc");}
mwSize   mxGetNzmax(const mxArray *pa) {
  mex2c_unimplemented("mxGetNzmax"); return 0;}
void     mxSetNzmax(mxArray *pa, mwSize nzmax) {
  mex2c_unimplemented("mxSetNzmax");}

/* get/set structs and cells */
int      mxGetNumberOfFields(const mxArray *pa) {
  mex2c_unimplemented("mxGetNumberOfFields"); return 0;}
mxArray *mxGetCell(const mxArray *pa, mwIndex i) {
  mex2c_unimplemented("mxGetCell"); return NULL;}
void     mxSetCell(mxArray *pa, mwIndex i, mxArray *value) {
  mex2c_unimplemented("mxSetCell");}
int      mxGetFieldNumber(const mxArray *pa, const char *name) {
  mex2c_unimplemented("mxGetFieldNumber"); return 0;}
mxArray *mxGetFieldByNumber(const mxArray *pa, mwIndex i, int fieldnum) {
  mex2c_unimplemented("mxGetFieldByNumber"); return NULL;}
void     mxSetFieldByNumber(mxArray *pa, mwIndex i, int num, mxArray *value) {
  mex2c_unimplemented("mxSetFieldByNumber");}
mxArray *mxGetField(const mxArray *pa, mwIndex i, const char *fieldname) {
    mex2c_unimplemented("mxGetField"); return NULL;}
void     mxSetField(mxArray *pa, mwIndex i, const char *name, mxArray *val) {
  mex2c_unimplemented("mxSetField");}
const char *mxGetFieldNameByNumber(const mxArray *pa, int n) {
  mex2c_unimplemented("mxGetFieldNameByNumber"); return NULL;}
int      mxSetClassName(mxArray *pa, const char *classname) {
  mex2c_unimplemented("mxSetClassName"); return 0;}
int      mxAddField(mxArray *pa, const char *fieldname) {
  mex2c_unimplemented("mxAddField"); return 0;}
void     mxRemoveField(mxArray *pa, int field) {
  mex2c_unimplemented("mxRemoveField");}

/* get/set objects */
mxArray *mxGetProperty(const mxArray *pa, mwIndex i, const char *propname) {
  mex2c_unimplemented("mxGetProperty"); return NULL;}
void     mxSetProperty(mxArray *pa, mwIndex i, const char *propname, const mxArray *value) {
  mex2c_unimplemented("mxSetProperty");}

