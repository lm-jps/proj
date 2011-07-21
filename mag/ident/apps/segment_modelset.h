/*
 * Name: segment_modelset.h
 *
 * Include file for model sets used for image segmentation
 *
 * No code changes are needed to define a new modelset.  You can do so
 * either in the .c file, to statically allocate it, or in a JSON format file 
 * that you name.  See the .c file for more on how this works.
 *
 * Michael Turmon, JPL
 *  created 6/2011
 *
 */

#ifndef _segment_modelset_h_
#define _segment_modelset_h_

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <sys/stat.h>
// JSON structure parsing: see http://sourceforge.net/projects/cjson/
// #include "cJSON.h"

// name field starts this way for a builtin model
#define SEG_BUILTIN_PREFIX "/builtin"

/***********************************************************************
 *
 * MODEL STRUCTURE DECLARATIONS
 *
 ***********************************************************************/

// a model for one region type
typedef struct {
  char *name;      // name for this class, e.g. quiet Sun
  int k;           // number of gaussian components making the class
  double *params;  // list of k, six-tuples of gaussian parameters
} seg_onemodel_t;

// a model-set for an entire image-classification method
// Note: The code currently assumes one component of the modelset
// is for quiet sun, and 
// identifies it by the presence of "quiet" somewhere in the name.
// The code also assumes one component is for active regions, and 
// identifies it by the presence of "active" somewhere in the name.
// Case is significant.
// So, it is best if some models[r].name contains "quiet" and another
// contains "active".
typedef struct {
  char *name;      // identifying name for the whole model
  char *descrip;   // readable description for the whole model
  char *version;   // version identifier
  int nclass;      // number of classes, e.g. 2 for quiet Sun vs. active region
  int ndim;        // number of image observations (always ndim = 2 = M + Ic)
  char *format;    // ordering of params (always = "var,chol,diag-then-upper")
  double *alpha;   // length-nclass vector of per-class probabilities
  seg_onemodel_t models[]; // length-nclass list of models, one for each class
} seg_modelset_t;

/***********************************************************************
 *
 * RUDIMENTARY MODEL API
 *
 ***********************************************************************/

// parse json text into modelset object
seg_modelset_t *seg_text2modelset(char *text, char **errmsg);
// Read a file, return a modelset for that file
seg_modelset_t *seg_modelfile2modelset(const char *filename, char **errmsg);
// frees allocated seg_modelset_t (not a statically-allocated modelset!)
void seg_modelset_destroy(seg_modelset_t *m);
// get a modelset, built-in or external
seg_modelset_t *seg_modelsetname2modelset(const char *name, char **errmsg, int enable_files);
// return the i'th built-in modelset, or NULL
seg_modelset_t *seg_modelset_builtin(int i);
// get the number, within a modelset, of a certain class (e.g., quiet)
int seg_modelname2modelnum(seg_modelset_t *modelset, char *modelname);
// print out a modelset
void seg_modelset_fprintf(FILE *fp, int verbose, seg_modelset_t *m);
// print out all builtin modelsets
void seg_modelsets_fprintf(FILE *fp, int verbose);

#endif // _hmi_segment_models_h_

