/*
 *  smpl_00.c						$DRMS/proj/cookbook/
 *
 *  An extremely simple program that does (almost) nothing at all.
 *    It illustrates how a module is called from a driver program, in this
 *    case a locally-provided main program rather than jsoc_main, and the
 *    use of the command line parsing features in (and outside of) the
 *    module.
 *
 *  Usage:
 *    smpl_00 [print= ...]
 *
 *  Bugs:
 *    The program is of no particular use, and exists merely for heuristic
 *      and testing purposes.
 *
 *  Revision history is at end of file.
 */
#include <cmdparams.h>
#include <timeio.h>
				    /*  sample module arguments declaration  */
ModuleArgs_t module_args[] = {
				   /*  module-specific argument declarators  */
  {ARG_STRING,	"print", "done", "message to print on successful completion"},
				       /*  required end (or blank) argument  */
  {ARG_END}
};
       /*  required global declaration (normally included in module driver)  */
CmdParams_t cmdparams;

int DoIt (void) {
  int status;
  char *msg = cmdparams_get_str (&cmdparams, "print", &status);
  printf ("%s\n", msg);
  return 0;
}

int main (int argc, char **argv) {
  int status = cmdparams_parse (&cmdparams, argc, argv);
  if (status >= 0) return DoIt ();
}

/*
 *  Revision History
 *
 *  09.07.27	file created by R Bogart
 */
