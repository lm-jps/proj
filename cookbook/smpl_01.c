/*
 *  smpl_01.c						$DRMS/proj/cookbook/
 *
 *  An extremely simple module that does (almost) nothing at all.
 *    It illustrates the structure of a DRMS module and exhibits
 *    return status behavior with use of the module verbose flags
 *
 *  Usage:
 *    smpl_01 [-avVH]
 *
 *  Bugs:
 *    The module is of no particular use, and exists merely for heuristic
 *      and testing purposes.
 *
 *  Revision history is at end of file.
 */
						  /*  required .h inclusion  */
#include <jsoc_main.h>
					     /*  required module identifier  */
char *module_name = "CookbookRecipe:01";
					    /*  optional version identifier
		      (recommended, and required for this particular module  */
char *version_id = "1.0";
				  /*  required module arguments declaration  */
ModuleArgs_t module_args[] = {
				   /*  module-specific argument declarators  */
  {ARG_FLAG,    "a",    "", "force an abort"},
  {ARG_FLAG,    "v",    "", "run in verbose mode"},
  {ARG_STRING,	"print", "done", "message to print on successful completion"},
				       /*  required end (or blank) argument  */
  {ARG_END}
};
					    /*  required module declaration  */
int DoIt (void) {
  int status;

  char *msg = strdup (cmdparams_get_str (&cmdparams, "print", &status));

  if (params_isflagset (&cmdparams, "v"))
    printf ("running module %s version %s\n", module_name, version_id);
  if (params_isflagset (&cmdparams, "a")) {
    printf ("aborting\n");
    return 1;
  }

  printf ("%s\n", msg);
				    /*  required status return to jsoc_main  */
  return 0;
}

/*
 *  Revision History
 *
 *  09.04.13	file created by R Bogart
 */
