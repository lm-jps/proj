/*
 * Name: segment_modelset.c
 *
 * Source file for model sets used for image segmentation
 *
 * No code changes are needed to define a new modelset.  You can do so
 * either here, if you want to statically allocate it as a built-in 
 * data structure, or in a JSON format file that you name.
 *
 * To use a JSON file, you create a file with the attributes here, and refer 
 * to it using its filename.  The file looks like this:
 *
 * {"modelset": { 
 *   "name":     "hmi.M_Ic_720s.test",
 *   "descrip":  "Test model",
 *   "version":  "1 from 2010.10.10",
 *   "nclass":   2,
 *   "ndim":     2,
 *   "format":   "var,diag-then-upper",
 *   "alpha":    [0.0, -4.0],
 *   "models": [
 * 	{ "name": "quiet Sun", "k": 3, "params": [ (list of k*6 doubles) ] }, 
 *	{ "name": "active",    "k": 5, "params": [ (list of k*6 doubles) ] }
 *    ]}}
 *
 * Or, to use a built-in model, name a pseudo-file referring to a statically-
 * allocated structure like the ones in this file.
 *
 * The naming scheme we have used for the name (in the modelset structure
 * below) uses the form:
 *   /builtin/<instrument>.<observables>.<type>
 * where instrument is hmi, observables is M_Ic_720s (for example), and
 * type is "production" or "test" depending on the status of the model.
 *
 * The /builtin prefix just means that the model is statically-allocated
 * in the source code.  Its symbolic value is in SEG_BUILTIN_PREFIX
 *
 * To add a model having K classes, define K static arrays, one for each class,
 * with separate names (K arrays -> K names).
 * Be sure one model has "quiet" in it (quiet sun is special)
 * The number of Gaussian components for class k (1 thru K) is, say, R_k.
 * The array for class k is then of length 6 * R_k.
 * Lastly, you need a prior probability array of length K (alpha).
 *
 * Then, insert these arrays into a seg_modelset_t structure (following
 * the examples here).  Finally you add your new seg_modelset_t structure
 * to the list of models near the end of the declarations here.
 *
 * Michael Turmon, JPL
 *  created 6/2011
 *
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <sys/stat.h>
// JSON structure parsing: see http://sourceforge.net/projects/cjson/
#include "cJSON.h"
// current interface
#include "segment_modelset.h"


/***********************************************************************
 *
 * BUILT-IN MODELS
 *
 * These are statically-allocated structures that we build up in a 
 * series of definitions here, from simple vectors to nested structures.
 *
 ***********************************************************************/

/*
 *
 * MODEL PARAMETER DEFINITIONS
 *
 * 6-by-R matrices of gaussian parameters, for each class and each model set.
 *
 */

// valid for Yang's test images 
// only for test, never used in production
static double seg_model_qs_test_yang[] = {
  0.204635504622,  1.20531725063, 1.00777710322,   19.346252856, 0.00419579448962,  0.0564904666445, 
  0.204635504622, -1.20531725063, 1.00777710322,   19.346252856, 0.00419579448962, -0.0564904666445, 
  0.232398997138,  0,             1.02955618645,   2.7404836868, 0.00427555589387,  0, 
  0.145253155252,  0,             0.992142854282, 247.248235524, 0.00504149434045,  0, 
  0.139157664971,  0,             0.950374859774, 74.2478459635, 0.000555807820391, 0, 
  0.0617507192555, 0,             0.986125619938, 2854.35303874, 0.00151332362565,  0, 
  0.01216845414,   0,             1.01070360771,  89583.9686521, 0.0145007628344,   0};
static double seg_model_ar_test_yang[] = {
  0.192105993124,  -846.854676681, 0.743104620169, 118603.827293, 0.00499678354147,   10.6891118764, 
  0.192105993124,   846.854676681, 0.743104620169, 118603.827293, 0.00499678354147,  -10.6891118764, 
  0.126980271901,   311.120029981, 0.773004292296, 18546.8591144, 0.00277679162354,  0.270699599935, 
  0.126980271901,  -311.120029981, 0.773004292296, 18546.8591144, 0.00277679162354, -0.270699599935, 
  0.111840866217,   1725.55076542, 0.432208758685,  148513.00487,  0.0244679802074,  -43.2940104149, 
  0.111840866217,  -1725.55076542, 0.432208758685,  148513.00487,  0.0244679802074,   43.2940104149, 
  0.0690728687575, -2565.81843681, 0.13790051683,  124931.169245,  0.0014162131028,   10.5502256214, 
  0.0690728687575,  2565.81843681, 0.13790051683,  124931.169245,  0.0014162131028,  -10.5502256214};
static double seg_model_alpha_test_yang[] = {0, -4.0};


// valid for M_720s and Ic_720s images beginning late march 2010
// used in production during 2010-2011
static double seg_model_qs_M_Ic_720s[] = {
  0.676700379049,  0, 13321.3804543, 137.665818814, 133580569.057, 0, 
  0.30516524961,   0, 13833.1964163, 214.193316404, 127827406.447, 0, 
  0.0181343713409, 0, 13300.8328273, 4715.87050328, 147742870.796, 0};
static double seg_model_ar_M_Ic_720s[] = {
  0.373249775824, -534.166827748, 13470.5783991, 87955.4275489, 140934500.036, -21833.1675272, 
  0.373249775824,  534.166827748, 13470.5783991, 87955.4275489, 140934500.036,  21833.1675272, 
  0.118699362858,  194.068579955, 13621.338073,  6919.60349716, 150211728.679, -33384.7350668, 
  0.118699362858, -194.068579955, 13621.338073,  6919.60349716, 150211728.679,  33384.7350668, 
  0.0161017226367, 0,             13038.102521,  165305.443351, 73124251.7542,  0};
static double seg_model_alpha_M_Ic_720s[] = {0, -4.0};

// in development for Ic_noLimbDark_720s, June 2011
// (AR class proved to be too permissive, discontinued November 2011)
// scaled from M_Ic_720s models above by 13500 along Ic coordinate
static double seg_model_qs_M_Ic_noLD_V1_720s[] = {
  0.6767003790, 0.0000000000, 0.9867689225,   137.6658188140, 0.7329523679, 0.0000000000, 
  0.3051652496, 0.0000000000, 1.0246812160,   214.1933164040, 0.7013849462, 0.0000000000, 
  0.0181343713, 0.0000000000, 0.9852468761,  4715.8705032800, 0.8106604708, 0.0000000000};
static double seg_model_ar_M_Ic_noLD_V1_720s[] = {
  0.3732497758, -534.1668277480, 0.9978206222, 87955.4275489000, 0.7733031552, -1.6172716687, 
  0.3732497758,  534.1668277480, 0.9978206222, 87955.4275489000, 0.7733031552,  1.6172716687, 
  0.1186993629,  194.0685799550, 1.0089880054,  6919.6034971600, 0.8242070161, -2.4729433383, 
  0.1186993629, -194.0685799550, 1.0089880054,  6919.6034971600, 0.8242070161,  2.4729433383, 
  0.0161017226,    0.0000000000, 0.9657853719, 165305.443351000, 0.4012304623,  0.0000000000};
static double seg_model_alpha_M_Ic_noLD_V1_720s[] = {0, -4.0};

// production model for Ic_noLimbDark_720s, October 2011
// scaled from original M_IC_noLD_720s model by 1.3 along M coordinate
// latter had been already, in turn:
// scaled from M_Ic_720s models above by 13500 along Ic coordinate
static double seg_model_qs_M_Ic_noLD_720s[] = {
  0.6767003790, 0.0000000000, 0.9867689225,  232.65523379566, 0.7329523679, 0.0000000000,
  0.3051652496, 0.0000000000, 1.0246812160,  361.98670472276, 0.7013849462, 0.0000000000,
  0.0181343713, 0.0000000000, 0.9852468761, 7969.82115054320, 0.8106604708, 0.0000000000};
static double seg_model_ar_M_Ic_noLD_720s[] = {
  0.3732497758, -694.4168760724, 0.9978206222, 148644.672557641, 0.7733031552, -2.1024531693,
  0.3732497758,  694.4168760724, 0.9978206222, 148644.672557641, 0.7733031552,  2.1024531693,
  0.1186993629,  252.2891539415, 1.0089880054,  11694.129910200, 0.8242070161, -3.2148263398,
  0.1186993629, -252.2891539415, 1.0089880054,  11694.129910200, 0.8242070161,  3.2148263398,
  0.0161017226,    0.0000000000, 0.9657853719, 279366.199263190, 0.4012304623,  0.0000000000};
static double seg_model_alpha_M_Ic_noLD_720s[] = {0, -4.0};

// in development for masks for the mdi resident archive, August 2013.
// scaled from hmi-production M_IC_noLD_720s model (just above) by 1.4 along M coordinate,
//   in accordance with Table 1 in the Liu et al., 2012, HMI/MDI calibration paper.  
//   (This 1.4x is in addition to the 1.3x used by the model we used as a starting point.)
static double seg_model_qs_mdi_fdM_II_96m[] = {
  0.6767003790, 0.0000000000, 0.9867689225,  456.00425823949, 0.7329523679, 0.0000000000,
  0.3051652496, 0.0000000000, 1.0246812160,  709.49394125661, 0.7013849462, 0.0000000000,
  0.0181343713, 0.0000000000, 0.9852468761,15620.84945506467, 0.8106604708, 0.0000000000};
static double seg_model_ar_mdi_fdM_II_96m[] = {
  0.3732497758, -972.18362650136, 0.9978206222, 291343.55821297632, 0.7733031552, -2.94343443702,
  0.3732497758,  972.18362650136, 0.9978206222, 291343.55821297632, 0.7733031552,  2.94343443702,
  0.1186993629,  353.20481551810, 1.0089880054,  22920.49462399200, 0.8242070161, -4.50075687572,
  0.1186993629, -353.20481551810, 1.0089880054,  22920.49462399200, 0.8242070161,  4.50075687572,
  0.0161017226,    0.00000000000, 0.9657853719, 547557.75055585231, 0.4012304623,  0.00000000000};
static double seg_model_alpha_mdi_fdM_II_96m[] = {0, -4.0};


/*
 *
 * STATIC MODEL STRUCTURE DEFINTIONS
 *
 * Fill in modelset structures with above parameter blocks
 * Note: all static models should have a name beginning with 
 * the string in `SEG_BUILTIN'.
 *
 */

// one row of gaussian parameters, for convenience
#define SEG_MODEL_ONEROW (6 * sizeof(double))

static seg_modelset_t 
seg_model_test_yang = {
  SEG_BUILTIN_PREFIX "/hmi.M_Ic_proxy.test",
  "Pre-launch model for proxy images made by Yang Liu",
  "1 from 2010.03.01",
  2, 
  2, "var,diag-then-upper", // always like this
  seg_model_alpha_test_yang, 
  { 
    // this is a list of seg_onemodel_t structures
    {"quiet Sun",     
      sizeof(seg_model_qs_test_yang)/SEG_MODEL_ONEROW, 
      seg_model_qs_test_yang}, 
    {"active region", 
     sizeof(seg_model_ar_test_yang)/SEG_MODEL_ONEROW, 
     seg_model_ar_test_yang}
  }
};

static seg_modelset_t 
seg_model_M_Ic_720s = {
  SEG_BUILTIN_PREFIX "/hmi.M_Ic_720s.production",
  "First post-launch model using hmi.Ic_720s and hmi.M_720s",
  "1 from 2010.04.15",
  2, 
  2, "var,diag-then-upper", // always like this
  seg_model_alpha_M_Ic_720s, 
  {
    // this is a list of seg_onemodel_t structures
    {"quiet Sun",     
     sizeof(seg_model_qs_M_Ic_720s)/SEG_MODEL_ONEROW, 
     seg_model_qs_M_Ic_720s}, 
    {"active region", 
     sizeof(seg_model_ar_M_Ic_720s)/SEG_MODEL_ONEROW, 
     seg_model_ar_M_Ic_720s}
  }
};

static seg_modelset_t 
seg_model_M_Ic_noLD_V1_720s = {
  SEG_BUILTIN_PREFIX "/hmi.M_Ic_noLimbDark_720s.version1",
  "Model using hmi.Ic_noLimbDark_720s and hmi.M_720s, June 2011, obsolete",
  "1 from 2011.06.20",
  2, 
  2, "var,diag-then-upper", // always like this
  seg_model_alpha_M_Ic_noLD_V1_720s,
  {
    // this is a list of seg_onemodel_t structures
    {"quiet Sun",     
     sizeof(seg_model_qs_M_Ic_noLD_V1_720s)/SEG_MODEL_ONEROW, 
     seg_model_qs_M_Ic_noLD_V1_720s}, 
    {"active region", 
     sizeof(seg_model_ar_M_Ic_noLD_V1_720s)/SEG_MODEL_ONEROW, 
     seg_model_ar_M_Ic_noLD_V1_720s}
  }
};

static seg_modelset_t 
seg_model_M_Ic_noLD_720s = {
  SEG_BUILTIN_PREFIX "/hmi.M_Ic_noLimbDark_720s.production",
  "Model using hmi.Ic_noLimbDark_720s and hmi.M_720s",
  "1 from 2011.10.10",
  2, 
  2, "var,diag-then-upper", // always like this
  seg_model_alpha_M_Ic_noLD_720s,
  {
    // this is a list of seg_onemodel_t structures
    {"quiet Sun",     
     sizeof(seg_model_qs_M_Ic_noLD_720s)/SEG_MODEL_ONEROW, 
     seg_model_qs_M_Ic_noLD_720s}, 
    {"active region", 
     sizeof(seg_model_ar_M_Ic_noLD_720s)/SEG_MODEL_ONEROW, 
     seg_model_ar_M_Ic_noLD_720s}
  }
};

static seg_modelset_t 
seg_model_mdi_fd_M_II_96m = {
  SEG_BUILTIN_PREFIX "/mdi.fd_M_II_96m.production",
  "Model using su_turmon.mdi_II_96m and mdi.fd_M_96m_lev182",
  "1 from 2013.09.01",
  2, 
  2, "var,diag-then-upper", // always like this
  seg_model_alpha_mdi_fdM_II_96m,
  {
    // this is a list of seg_onemodel_t structures
    {"quiet Sun", 
     sizeof(seg_model_qs_mdi_fdM_II_96m)/SEG_MODEL_ONEROW, 
     seg_model_qs_mdi_fdM_II_96m}, 
    {"active region", 
     sizeof(seg_model_ar_mdi_fdM_II_96m)/SEG_MODEL_ONEROW, 
     seg_model_ar_mdi_fdM_II_96m}
  }
};


// keep the namespace clean
#undef SEG_MODEL_ONEROW

// Add new static models to this list. 
// No further code modifications are needed.
static seg_modelset_t *
seg_models_known[] = {
  &seg_model_test_yang, 
  &seg_model_M_Ic_720s, 
  &seg_model_M_Ic_noLD_V1_720s,  // unused in production
  &seg_model_M_Ic_noLD_720s, 
  &seg_model_mdi_fd_M_II_96m,    // draft model for mdi
  NULL}; // MUST be null-terminated


/***********************************************************************
 *
 * JSON MODEL INTERFACE
 *
 ***********************************************************************/

/*
 * Functions for retrieving objects from cJSON structures
 *
 * All these functions have these semantics upon error:
 *    (1) sets *status argument to nonzero, and
 *    (2) returns NULL or 0, and
 *    (3) puts a message into a global error buffer
 *
 * Also, for all these routines, a nonzero *status upon entry 
 * means there was already an error, so we immediately and 
 * unconditionally return NULL or 0.  
 * The error buffer is not disturbed, so the message is preserved.
 *
 */

/*
 * Using arguments identical to printf, puts an error message
 * into a global char array.
 *
 * Except: if fmt is given as NULL, returns the current buffer.
 *
 * You must free this string when you're done with it.
 */

static char seg_err_buf[200]; // current error message

static
char *
seg_err_log(const char *fmt, ...)
{
  va_list ap;
  char *rv;   // return value

  va_start(ap, fmt);
  if (fmt) {
    vsnprintf(seg_err_buf, sizeof(seg_err_buf), fmt, ap);
    rv = NULL;
  }
  else {
    rv = seg_err_buf;
  }
  va_end(ap);
  return rv;
}

/*
 * Return a string in the field `name' from the object `obj'
 *
 * You must free this string when you're done with it.
 */
static
char *
retrieve_str(cJSON *obj, char *name, char *parent, int *status)
{
  cJSON *field;
  char *str;

  if (*status) return NULL; // existing error condition
  if ((field = cJSON_GetObjectItem(obj, name)) == NULL) {
    *status = 1;
    seg_err_log("string: no field `%s' in `%s'.", name, parent);
    return NULL;
  }
  if (field->type != cJSON_String) {
    *status = 1;
    seg_err_log("string: field `%s' in `%s' was not a string.", name, parent);
    return NULL;
  }
  // we must retain our own copy
  if ((str = strdup(field->valuestring)) == NULL) {
    *status = 1;
    seg_err_log("string: out of memory for field `%s' in `%s'.", name, parent);
    return NULL;
  }
  return str;
}

/*
 * Return an integer in the field `name' from the object `obj'
 *
 * If error, returns 0.
 */
static
int
retrieve_int(cJSON *obj, char *name, char *parent, int *status)
{
  cJSON *field;

  if (*status) return 0; // existing error condition
  if ((field = cJSON_GetObjectItem(obj, name)) == NULL) {
    *status = 1;
    seg_err_log("int: no field `%s' in `%s'.", name, parent);
    return 0;
  }
  if (field->type != cJSON_Number) {
    *status = 1;
    seg_err_log("int: field `%s' in `%s' was not an int.", name, parent);
    return 0;
  }
  return field->valueint;
}

/*
 * Return an object in the field `name' from the object `obj'
 *
 * It returns a reference, not an independent copy, so no need to free.
 *
 */
static
cJSON *
retrieve_obj(cJSON *obj, char *name, char *parent, int *status)
{
  cJSON *obj2;

  if (*status) return NULL; // existing error condition
  if ((obj2 = cJSON_GetObjectItem(obj, name)) == NULL) {
    *status = 1;
    seg_err_log("obj: no field `%s' in `%s'.", name, parent);
    return NULL;
  }
  return obj2;
}

/*
 * Return the i'th item in the array `obj'
 *
 * It returns a reference, not an independent copy, so no need to free.
 *
 */
static
cJSON *
retrieve_itm(cJSON *obj, int i, char *parent, int *status)
{
  cJSON *obj2;

  if (*status) return NULL; // existing error condition
  if (obj->type != cJSON_Array) {
    *status = 1;
    seg_err_log("item: parent object `%s' was not an array.", parent);
    return NULL;
  }
  if ((obj2 = cJSON_GetArrayItem(obj, i)) == NULL) {
    *status = 1;
    seg_err_log("item: no item #%d in parent `%s'.", i+1, parent);
    return NULL;
  }
  return obj2;
}

/*
 * Return a double array in the field `name' from the object `obj'
 * The length should be len, but if len < 0, the check is suppressed.
 *
 * You must free this array when you're done with it.
 */
static
double *
retrieve_arr(cJSON *obj, char *name, char *parent, int len, int *status)
{
  int len1, i;
  double *arr;
  cJSON *field, *entry;

  if (*status) return NULL; // existing error condition
  if ((field = cJSON_GetObjectItem(obj, name)) == NULL) {
    *status = 1;
    seg_err_log("array: no field `%s' in `%s'.", name, parent);
    return NULL;
  }
  if (field->type != cJSON_Array) {
    *status = 1;
    seg_err_log("array: field `%s' in `%s' was not an array.", name, parent);
    return NULL;
  }
  len1 = cJSON_GetArraySize(field);
  if (len < 0) {
    len = len1; // let the predicted length be the observed length
  } else if (len1 != len) {
    // not the expected length
    *status = 1;
    seg_err_log("array: field `%s' in `%s' had length %d, expected %d.", 
		name, parent, len1, len);
    return NULL;
  }
  if ((arr = calloc(len, sizeof(double))) == NULL) {
    *status = 1;
    seg_err_log("array: field `%s' in `%s': failed calloc (%d doubles).", 
		name, parent, len);
    return NULL;
  }
  // loop to get values
  for (i = 0; i < len; i++) {
    entry = cJSON_GetArrayItem(field, i);
    if (!entry) {
      *status = 1;
      seg_err_log("array: field `%s' in `%s': could not find item #%d of %d.", 
		  name, parent, i+1, len);
      return NULL;
    }
    if (entry->type != cJSON_Number) {
      *status = 1;
      seg_err_log("array: field `%s' in `%s': expected item #%d of %d to be numeric", 
		  name, parent, i+1, len);
      return NULL;
    }
    arr[i] = entry->valuedouble;
  }
  return arr;
}

/*
 * extract modelset structure from json parser structure
 *
 * returns allocated modelset structure, or NULL.
 * if modelset is NULL, an error message is in *errmsg.
 */
static
seg_modelset_t *
extract_modelset(cJSON *root)
{
  int i, Nc;
  seg_modelset_t *m = NULL;
  seg_onemodel_t *m1;
  cJSON *JS_set;         // modelset object at top of object
  cJSON *JS_models;      // model array within modelset
  cJSON *JS_m1;          // one model
  int status = 0;        // OK
  char str_m1[40];

  // ensure an empty error message to start
  seg_err_log("");
  // need number of classes to make modelset structure
  JS_set      = retrieve_obj(root,   "modelset","root",     &status);
  Nc          = retrieve_int(JS_set, "nclass",  "modelset", &status);
  // check status now, to ensure a valid Nc before calloc
  if (status) return NULL; // error message is already set up
  // modelset struct, plus one struct for each model in models[]
  if ((m = calloc(1, 
		  sizeof(seg_modelset_t) + 
		  sizeof(seg_onemodel_t)*Nc)) == NULL) {
    seg_err_log("extract: failed calloc for modelset");
    return NULL;
  }
  m->nclass   = Nc; // got this above
  m->name     = retrieve_str(JS_set, "name",    "modelset", &status);
  m->descrip  = retrieve_str(JS_set, "descrip", "modelset", &status);
  m->version  = retrieve_str(JS_set, "version", "modelset", &status);
  m->ndim     = retrieve_int(JS_set, "ndim",    "modelset", &status);
  m->format   = retrieve_str(JS_set, "format",  "modelset", &status);
  m->alpha    = retrieve_arr(JS_set, "alpha",   "modelset", Nc, &status);
  JS_models   = retrieve_obj(JS_set, "models",  "modelset", &status);
  for (i = 0; i < Nc; i++) {
    // a descriptor for the containing object
    snprintf(str_m1, sizeof(str_m1), "models(%d) in modelset", i+1);
    // the model we're setting
    m1 = m->models+i;
    JS_m1     = retrieve_itm(JS_models, i,    str_m1, &status);
    m1->name  = retrieve_str(JS_m1, "name",   str_m1, &status);
    m1->k     = retrieve_int(JS_m1, "k",      str_m1, &status);
    m1->params= retrieve_arr(JS_m1, "params", str_m1, (m1->k)*6, &status);
  }
  if (status) {
    seg_modelset_destroy(m);
    return NULL;
  } else {
    return m;
  }
}


/***********************************************************************
 *
 * RUDIMENTARY MODEL API
 *
 ***********************************************************************/

/* 
 * parse json text into modelset object
 *
 * If error, NULL is returned, and *errmsg contains a descriptive
 * message.
 *
 */
seg_modelset_t *
seg_text2modelset(char *text, char **errmsg)
{
  const int errmsg_max = 140; // longest explanation we want
  const char *msg_low; // low-level error message
  int low_len;         // its length
  const char *msg_fmt; // format of *errmsg
  char *msg_ez;        // last-ditch message
  char *suffix;
  cJSON *json;
  seg_modelset_t *m;
  
  json = cJSON_Parse(text);
  if (json == NULL) {
    m = NULL; // signal failure later
    // capture the low-level error message, reformat later
    msg_low = cJSON_GetErrorPtr();
    msg_fmt = "JSON parse error before the string:\n%.96s%s";
    msg_ez  = "JSON parse error. Could not allocate memory for description";
  } else {
    m = extract_modelset(json);
    cJSON_Delete(json);
    if (!m) {
      msg_low = seg_err_log(NULL); // as above, capture local-eror buffer
      msg_fmt = "text2modelset: Could not retrieve %.96s%s";
      msg_ez  =  "text2modelset: Could not retrieve an element";
    }
  }
  // set up error message
  if (m) {
    *errmsg = NULL;
  } else {
    // reformat the low-level error message
    low_len = strlen(msg_low);
    if (low_len > 96) suffix = " [...]"; else suffix = "";
    if ((*errmsg = calloc(errmsg_max, 1)) == NULL) {
      *errmsg = msg_ez; // not even room for the calloc
    } else {
      snprintf(*errmsg, errmsg_max, msg_fmt, msg_low, suffix);
    }
  } 
  return m;
}

/* 
 * Read a file, return a modelset for that file.
 * 
 * Error message put into *errmsg.
 */
seg_modelset_t *
seg_modelfile2modelset(const char *filename, char **errmsg)
{
  long len;
  char *text;
  seg_modelset_t *m;
  FILE *fp;

  *errmsg = NULL; // assume no error
  if ((fp = fopen(filename,"r")) == NULL) {
    *errmsg = "Could not open file for reading";
    return NULL;
  }
  fseek(fp, 0, SEEK_END);
  len = ftell(fp);
  rewind(fp);
  if ((text = calloc(len+1, 1)) == NULL) {
    *errmsg = "Insufficient memory to read named file";
    fclose(fp);
    return NULL;
  }
  if (fread(text, 1, len, fp) != len) {
    *errmsg = "Failed to read all of named file";
    fclose(fp);
    return NULL;
  }
  fclose(fp);
  // done with the file, convert the text
  m = seg_text2modelset(text, errmsg);
  free(text);
  return m;
}

/*
 * Given a modelset name, look it up in the list of builtin sets, or
 * optionally allow the name to refer to a JSON modelset file.
 *
 * Text error message in *errmsg contains explanatory message if 
 * returned model is NULL, which signals error.
 */
seg_modelset_t *
seg_modelsetname2modelset(const char *name, char **errmsg, int enable_files)
{
  seg_modelset_t *m;
  int i;
  struct stat BUF;

  *errmsg = "";  // clear the error
  if (strstr(name, SEG_BUILTIN_PREFIX) == name) {
    // look thru builtin files
    for (i = 0; (m = seg_models_known[i]) != NULL; i++) {
      // printf("%d: %s\n", i, m->name);
      if (strcmp(name, m->name) == 0)
	return m;
    }
    *errmsg = "modelset given as builtin, but not found as builtin";
    return NULL;
  }
  // load a file, if allowed
  if (!enable_files) {
    m = NULL;
    *errmsg = "modelset not found as builtin";
  } else {
    // Try to read a file, and return a modelset for that file
    m = seg_modelfile2modelset(name, errmsg);
  }
  return m;
}

/*
 * frees allocated seg_modelset_t
 *
 * the reverse of: modelsetname2modelset
 *
 * Can be used with a statically-allocated modelset, as long as the
 * convention that the name of such a modelset begins with 
 * SEG_BUILTIN_PREFIX.
 */
void
seg_modelset_destroy(seg_modelset_t *m)
{
  seg_onemodel_t *m1;
  int i;

  if (!m) return;
  // don't destroy a builtin
  if (strstr(m->name, SEG_BUILTIN_PREFIX) == m->name) return;
  // alpha is a double *
  if (m->alpha)
    free(m->alpha);
  // models -- array of single model structures
  for (i = 0; i < m->nclass; i++) {
    m1 = m->models+i; // a single model
    if (m1 && m1->params)
      free(m1->params); // it's a double *
  }
  free(m);
}

/*
 * returns i'th seg_modelset_t
 *
 * returns NULL once the last valid modelset is passed;
 * it is illegal to proceed past the first NULL.
 *
 */
seg_modelset_t *
seg_modelset_builtin(int i)
{
  return seg_models_known[i];
}

/*
 * Given a single modelset and a class name, returns the index of the
 * class name in the modelset.  Used for, e.g., finding the number of
 * the quiet sun class.
 *
 * Return value is in 0..(modelset->nclass-1),
 * or -1 for "not found".
 * With 0 for off-disk, the class number recorded in the 
 * image would be one larger than the returned number.
 */
int 
seg_modelname2modelnum(seg_modelset_t *modelset, char *modelname) 
{
  int class;
  seg_onemodel_t *m1;

  for (class = 0; class < modelset->nclass; class++) {
    m1 = &(modelset->models[class]);
    if (strstr(m1->name, modelname) != 0)
      return class;
  }
  return -1;
}

// output format, for convenience
#define SEG_FLFMT "%.4g "

/*
 * Print a summary of a modelset structure to a file pointer, with
 * verbosity of 0 or 1.
 */

void
seg_modelset_fprintf(FILE *fp, int verbose, seg_modelset_t *m)
{
  int class, inx;
  int nbump, bump;
  double *params;
  seg_onemodel_t *m1;

  fprintf(fp, "Name = <%s>.\n\tVersion = %s.\n\t%s\n", 
	  m->name, m->version, m->descrip);
  fprintf(fp, "\tModel set has %d components in %d dimensions.\n", 
	  m->nclass, m->ndim);
  if (verbose)
    fprintf(fp, "\tModel set uses vectorization format %s.\n", m->format);
  if (verbose) {
    // if nonstandard m->mode, this order will vary
    fprintf(fp, "\tParams listed in this order: "
	    "wt, mu(B), mu(Ic), var(B), var(Ic), cov(B,Ic)\n");
    for (class = 0; class < m->nclass; class++) {
      m1 = &(m->models[class]);
      nbump = m1->k;
      params = m1->params;
      fprintf(fp, "\n\tFor %s:\n", m1->name);
      for (inx = bump = 0; bump < nbump; bump++, inx += 6) {
	fprintf(fp, 
		"\t"
		SEG_FLFMT SEG_FLFMT SEG_FLFMT 
		SEG_FLFMT SEG_FLFMT SEG_FLFMT
		"\n", 
		params[inx+0], 
		params[inx+1],
		params[inx+2],
		params[inx+3],
		params[inx+4],
		params[inx+5]);
      } // for inx
    } // for class
  } // if verbose
}

/*
 * Print a summary of all built-in modelset structures 
 * to a file pointer, with verbosity of 0 or 1.
 */

void
seg_modelsets_fprintf(FILE *fp, int verbose)
{
  int i;
  seg_modelset_t *m;

  for (i = 0; seg_models_known[i] != NULL; i++) {
    m = seg_models_known[i];
    fprintf(fp, "\n");
    seg_modelset_fprintf(fp, verbose, m);
  }
  fprintf(fp, "\n");
}

// clean namespace
#undef SEG_FLFMT


