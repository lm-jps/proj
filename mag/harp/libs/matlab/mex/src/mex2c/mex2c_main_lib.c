/*
 * mex2c_main_lib
 * 
 * This is a shell around the matlab/mex "mexFunction" conventional
 * function entry point.  This shell allows
 * functions written for compilation and execution
 * under the mex environment to be used standalone from C, 
 * without the rest of the matlab environment, hence the name
 * mex2c.  The cli part indicates the arguments to the mexfunction
 * are passed in by a command-line interface.
 * This works by using the command line to instantiate the Matrix
 * data structures that are used by the mexFunction entry point.
 *
 * version 1, Michael Turmon, 2010
 */

/*LINTLIBRARY*/

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <setjmp.h> /* for mexFunction exit via mexErrMsgTxt */

/*
 * mex includes
 */
#include "mex.h"        /* mex data structures */
#include "mexhead.h"    /* mex utilities */


/* global variables to hold program state */
char *mex2c_progname;  /* name of this program, for error messages;
			  is also used for mexFunctionName() in
			  mex_stubs.c */
char mex2c_phase[256]; /* for error messages: where are we */

// the error handler and dispatcher wrapper for CLI and mex-as-library
#include "mex2c_dispatcher.c"



