/*
 * mex2c_dispatcher
 * 
 * This is some common code that wraps a mexFunction in a dispatcher
 * which handles error recovery via setjmp/longjmp.
 *
 * It is included in two other files here, one for making a C library
 * from a mexFunction, and one for making a command-line executable
 * from a mexFunction.
 *
 * Michael Turmon, 2010
 */

/* no system #includes because this file is itself #included */

/* assumes the following global variables to hold program state:
 *
 * char *mex2c_progname; 
 * char mex2c_phase[...]; 
 *
 */

/* hook for within-mexFunction exception handling via setjmp/longjmp */
static jmp_buf MXT_exception_env;
static int     MXT_exception_valid = 0; // is the above valid
static char    MXT_errstr[256]; // transmit error message from longjmp


/************************************************************************
 * mxt_ErrHandler
 * 
 * the mex2c versions of mexErrMsgTxt and friends route back thru this point.
 * 
 * This hook can be defined by each interface mode.  
 *
 ************************************************************************/
void
mxt_ErrHandler(const char *msg)
{
#ifdef DEBUG
  printf("  in mxt_ErrHandler, returning to %s.\n", 
	 MXT_exception_valid ? "caller" : "system");
#endif
  // jump back to the matching setjmp() where we called mexFunction()
  if (MXT_exception_valid) {
    // pass error message back, but keep quiet
    if (msg) {
      strncpy(MXT_errstr, msg, sizeof(MXT_errstr));
      MXT_errstr[sizeof(MXT_errstr)-1] = '\0'; // ensure nul-terminated
    }
    longjmp(MXT_exception_env, 1);
  } else {
    // if reached, this is a code error
    fprintf(stderr, "mxt_ErrHandler: no matching setjmp: should not happen.\n");
    if (msg)
      fprintf(stderr, "\tError message is: %s.\n", msg);
    fprintf(stderr, "\tReturning to system.\n");
    exit(1); // nowhere else to go
  }
}


/************************************************************************
 * 
 * mxt_mexDispatcher: wrapper for mexfunction execution from C
 * 
 * set jump marker, call mexFunction, with return on error via longjmp
 *
 ***********************************************************************/

char *
mxt_mexDispatcher(mexfn_t *mexfn, 
		  const char *progname,
		  int verbosity,
		  int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs)
{
  char *retval;
  jmp_buf save_MXT_exception_env;
  int save_MXT_exception_valid;
  int i;

  strcpy(mex2c_phase, "mexFunction dispatcher");
  mex2c_progname = strdup(progname);

  // save setjmp state (it may be invalid, that's OK)
  save_MXT_exception_valid = MXT_exception_valid;
  memcpy(save_MXT_exception_env, MXT_exception_env, sizeof(MXT_exception_env));

  /* call mexfunction */
  strcpy(mex2c_phase, "mexFunction internal");
  MXT_exception_valid = 1; // jmp_buf for longjmp is set
  if (setjmp(MXT_exception_env) == 0) {
#ifdef DEBUG
    printf("  calling mexfunction %s...\n", progname);
#endif
    (*mexfn)(nlhs, plhs, nrhs, prhs);
#ifdef DEBUG
    printf("  returned normally from mexfunction %s.\n", progname);
#endif
    // indicates OK
    retval = NULL;
  } else {
    // mexErrMsgTxt within mexFunction() will route thru this clause
    if (verbosity) {
      fprintf(stderr, "dispatcher: encountered error within mexfunction `%s'.\n", 
	      progname);
      fprintf(stderr, "\tError message returned by `%s' follows:\n\t%s\n", 
	      progname, MXT_errstr);
    }
    // null out plhs to ensure it is not used
    for (i = 0; i < nlhs; i++)
      plhs[i] = NULL;
    // pass back the error message; the pointer is non-NULL
    retval = MXT_errstr; 
  }
  // restore prior setjmp state
  MXT_exception_valid = save_MXT_exception_valid;
  memcpy(MXT_exception_env, save_MXT_exception_env, sizeof(MXT_exception_env));
  return retval;
}


