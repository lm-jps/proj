/*
 *  donought.c							~rick/src/cb
 *
 *  An extremely simple module, does (almost) nothing at all.
 *
 *  Responsible:  Rick Bogart			RBogart@solar.stanford.edu
 *
 *  Usage:
 *    donought [-avVH]
 *
 *  Bugs:
 *    The module is of no particular use, and exists merely for heuristic
 *      and testing purposes.
 *
 *  Revision history is at end of file.
 */
						  /*  required .h inclusion  */
#include <jsoc_main.h>
#include <stdio.h>
					  /*  optional(?) module identifier  */
char *module_name = "donought";
char *version_id = "0.9";
							 /*  arguments list  */
ModuleArgs_t module_args[] = {
  {ARG_FLAG,    "a",    "", "force an abort"},
  {ARG_FLAG,    "v",    "", "run in verbose mode"},
  {}
};

int DoIt (void) {
  int status;

  if (params_isflagset (&cmdparams, "v"))
    printf ("running module %s version %s\n", module_name, version_id);
  if (params_isflagset (&cmdparams, "a")) {
    printf ("aborting\n");
    return 1;
  }
  printf ("done\n");
  return 0;
}

/*
 *  Revision History
 *
 *  05.11.14	Rick Bogart	created this file
 */
