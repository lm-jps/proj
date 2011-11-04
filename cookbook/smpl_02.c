/*
 *  smpl_02.c						$DRMS/proj/cookbook/
 *
 *  Annother simple module that does nothing but echo back the values of
 *    arguments. It illustrates the use of command line parsing and the
 *    variety of argument types; also introduces time string functions
 *
 *  Usage:
 *    smpl_02 [-avVH]
 *
 *  Bugs:
 *    The module is of no particular use, and exists merely for heuristic
 *      and testing purposes.
 *
 *  Revision history is at end of file.
 */

#include <jsoc_main.h>
char *module_name = "CookbookRecipe:02";
char *version_id = "1.0";

ModuleArgs_t module_args[] = {
  {ARG_STRING,	"name",	"Not Specified", "a string"}, 
  {ARG_INT,	"ival",	"1", "a positive integer", "[1,)"},
  {ARG_INTS,	"iptr", "[0]", "array of integers"},
  {ARG_FLOAT,	"fval",	"0.0", "a real number"},
  {ARG_FLOAT,	"fvalns", "Unspecified", "a real number"},
  {ARG_FLOATS,	"fptr", "{2.71828, 3.14159, -1}", "array of real numbers"},
  {ARG_TIME,	"time", "1582.10.5_00:00:00",
      "a time, in standard date_time format"},
  {ARG_NUME, "colour", "", "enumerated choice without a default",
      "red, orange, yellow, green, blue, indigo, violet"},
  {ARG_NUME, "mois", "Brumaire", "enumerated choice with a default",
      "Vendémiaire, Brumaire, Frimaire, Nivôse, Pluviôse, Ventôse, Germinal, \
      Florial, Prairial, Messidor, Thermidor, Fructidor"},
  {ARG_FLAG,	"e",	"", "a flag value"}, 
  {ARG_END}
};

int DoIt (void) {
  CmdParams_t *params = &cmdparams;
  double *fpval;
  int *ipval;
  int i, fpvals, ipvals;
  char *colours[] = {"red", "orange", "yellow", "green", "blue", "violet"};
  char *moiss[] = {"Vendémiaire", "Brumaire", "Frimaire", "Nivôse", "Pluviôse",
      "Ventôse", "Germinal", "Florial", "Prairial", "Messidor", "Thermidor",
      "Fructidor"};
  char key[64], tbuf[64];

  TIME tval = params_get_time (params, "time");
  double fval = params_get_double (params, "fval");
  int ival = params_get_int (params, "ival");
  int colour = params_get_int (params, "colour");
  int mois = params_get_int (params, "mois");
  char *name = strdup (params_get_str (params, "name"));
  int flagset = params_isflagset (params, "e");

  printf ("name = %s\n", name);
  printf ("fval = %g\n", fval);
  printf ("ival = %d\n", ival);
  sprint_time (tbuf, tval, "UT", 3);
  printf ("time = %s\n", tbuf);
  printf ("colour = %d (%s)\n", colour, colours[colour]);
  printf ("mois = %d (%s)\n", mois, moiss[mois]);
  printf ("-e? : %d\n", flagset);

  printf ("fpvals = %d\n", fpvals = params_get_int (params, "fptr_nvals"));
  fpval = (double *)malloc (fpvals * sizeof (double));
  for (i = 0; i < fpvals; i++) {
    sprintf (key, "fptr_%d_value", i);
    fpval[i] = params_get_double (params, key);
    printf ("fpval[%d] = %g\n", i, fpval[i]);
  }
  printf ("ipvals = %d\n", ipvals = params_get_int (params, "iptr_nvals"));
  ipval = (int *)malloc (ipvals * sizeof (int));
  for (i = 0; i < ipvals; i++) {
    sprintf (key, "iptr_%d_value", i);
    ipval[i] = params_get_int (params, key);
    printf ("ipval[%d] = %d\n", i, ipval[i]);
  }

  return (0);
}

/*
 *  Revision History
 *
 *  07.02.27	created by RSB, based on original using SSSC API
 *  09.04.14	minor mods to illustrate additional features of arg parsing
 */
