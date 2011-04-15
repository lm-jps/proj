
#include <stdlib.h>
#include <string.h>
#include "cmdparams.h"
#define CPSAVE_SUCCESS			(0)
#define CPSAVE_UNKNOWN_PARAM		(1)
#define CPSAVE_INVALID_CONVERSION	(2)
#define CPSAVE_OUTOFMEMORY		(4)
#define CPSAVE_UNKNOWN_ERROR		(8)

#define SAVESTRUNIT 256
#define PARMSEPRTR "\n"
char *savestr=NULL;
int savestrmax=0;
int savestrlen=0;

float 	cmdparams_save_float (CmdParams_t *parms, char *name, int *status);
double 	cmdparams_save_double (CmdParams_t *parms, char *name, int *status);
int 	cmdparams_save_int (CmdParams_t *parms, char *name, int *status);
char * 	cmdparams_save_str (CmdParams_t *parms, char *name, int *status);
double 	cmdparams_save_time (CmdParams_t *parms, char *name, int *status);
int	cmdparams_save_flag (CmdParams_t *parms, char *name, int *status);
char *	cmdparams_save_arg (CmdParams_t *parms, int num, int *status);

void cpsave_decode_error(int status) {
  if (status == 0) return;
  if (status & CPSAVE_UNKNOWN_PARAM) fprintf(stderr, "CPSAVE: unknown parameter.\n");
  if (status & CPSAVE_INVALID_CONVERSION) fprintf(stderr, "CPSAVE: invalid conversion.\n");
  if (status & CPSAVE_OUTOFMEMORY) fprintf(stderr, "CPSAVE: out of memory.\n");
  if (status & CPSAVE_UNKNOWN_ERROR) fprintf(stderr, "CPSAVE: unknown error.\n");
}


float cmdparams_save_float (CmdParams_t *parms, char *name, int *status) {

  float retval;
  char *strval, *buf;
  int nadd, stat;
  int newstat=0;

  retval=cmdparams_get_float (parms, name, &stat);
  switch (stat) {
  case CMDPARAMS_SUCCESS:
    newstat=CPSAVE_SUCCESS;
    break;
  case CMDPARAMS_UNKNOWN_PARAM:
    newstat=CPSAVE_UNKNOWN_PARAM;
    break;
  case CMDPARAMS_INVALID_CONVERSION:
    newstat=CPSAVE_INVALID_CONVERSION;
    break;
  case CMDPARAMS_OUTOFMEMORY:
    newstat=CPSAVE_OUTOFMEMORY;
    break;
  default:
    newstat=CPSAVE_UNKNOWN_ERROR;
  }

  if (stat == CMDPARAMS_UNKNOWN_PARAM) {  // can't find the parameter so give up
    *status = *status | newstat;
    return retval;
  }

  if (savestrmax == 0) {
    savestr=calloc(SAVESTRUNIT,sizeof(char));
    if (savestr == NULL) {
      *status = *status | newstat | CPSAVE_OUTOFMEMORY;
      return retval;
    }
    savestrmax=SAVESTRUNIT;
  }

  strval=(char *)cmdparams_get_str (parms, name, NULL);  // already know the status is success
  nadd=strlen(strval)+strlen(name)+strlen(PARMSEPRTR)+2; 
  if (savestrlen+nadd > savestrmax) {
    savestrmax += (nadd > SAVESTRUNIT) ? nadd : SAVESTRUNIT;
    buf=malloc(savestrmax*sizeof(char));
    if (buf == NULL) {
      *status = *status | newstat | CPSAVE_OUTOFMEMORY;
      return retval;
    }
    else {
      strcpy(buf,savestr);
      free(savestr);
      savestr=buf;
    }
  }

  strcat(savestr,name);
  strcat(savestr,"=");
  strcat(savestr,strval);
  strcat(savestr,PARMSEPRTR);
  savestrlen+=nadd-1;

  *status = *status | newstat;
  return retval;
}


double cmdparams_save_double (CmdParams_t *parms, char *name, int *status) {

  double retval;
  char *strval, *buf;
  int nadd, stat;
  int newstat=0;

  retval=cmdparams_get_double (parms, name, &stat);
  switch (stat) {
  case CMDPARAMS_SUCCESS:
    newstat=CPSAVE_SUCCESS;
    break;
  case CMDPARAMS_UNKNOWN_PARAM:
    newstat=CPSAVE_UNKNOWN_PARAM;
    break;
  case CMDPARAMS_INVALID_CONVERSION:
    newstat=CPSAVE_INVALID_CONVERSION;
    break;
  case CMDPARAMS_OUTOFMEMORY:
    newstat=CPSAVE_OUTOFMEMORY;
    break;
  default:
    newstat=CPSAVE_UNKNOWN_ERROR;
  }

  if (stat == CMDPARAMS_UNKNOWN_PARAM) {
    *status = *status | newstat;
    return retval;
  }

  if (savestrmax == 0) {
    savestr=calloc(SAVESTRUNIT,sizeof(char));
    if (savestr == NULL) {
      *status = *status | newstat | CPSAVE_OUTOFMEMORY;
      return retval;
    }
    savestrmax=SAVESTRUNIT;
  }

  strval=(char *)cmdparams_get_str (parms, name, NULL);
  nadd=strlen(strval)+strlen(name)+strlen(PARMSEPRTR)+2;
  if (savestrlen+nadd > savestrmax) {
    savestrmax += (nadd > SAVESTRUNIT) ? nadd : SAVESTRUNIT;
    buf=malloc(savestrmax*sizeof(char));
    if (buf == NULL) {
      *status = *status | newstat | CPSAVE_OUTOFMEMORY;
      return retval;
    }
    else {
      strcpy(buf,savestr);
      free(savestr);
      savestr=buf;
    }
  }

  strcat(savestr,name);
  strcat(savestr,"=");
  strcat(savestr,strval);
  strcat(savestr,PARMSEPRTR);
  savestrlen+=nadd-1;

  *status = *status | newstat;
  return retval;
}


int cmdparams_save_int (CmdParams_t *parms, char *name, int *status) {

  int retval;
  char *strval, *buf;
  int nadd, stat;
  int newstat=0;

  retval=cmdparams_get_int (parms, name, &stat);
  switch (stat) {
  case CMDPARAMS_SUCCESS:
    newstat=CPSAVE_SUCCESS;
    break;
  case CMDPARAMS_UNKNOWN_PARAM:
    newstat=CPSAVE_UNKNOWN_PARAM;
    break;
  case CMDPARAMS_INVALID_CONVERSION:
    newstat=CPSAVE_INVALID_CONVERSION;
    break;
  case CMDPARAMS_OUTOFMEMORY:
    newstat=CPSAVE_OUTOFMEMORY;
    break;
  default:
    newstat=CPSAVE_UNKNOWN_ERROR;
  }

  if (stat == CMDPARAMS_UNKNOWN_PARAM) {
    *status = *status | newstat;
    return retval;
  }

  if (savestrmax == 0) {
    savestr=calloc(SAVESTRUNIT,sizeof(char));
    if (savestr == NULL) {
      *status = *status | newstat | CPSAVE_OUTOFMEMORY;
      return retval;
    }
    savestrmax=SAVESTRUNIT;
  }

  strval=(char *)cmdparams_get_str (parms, name, NULL);
  nadd=strlen(strval)+strlen(name)+strlen(PARMSEPRTR)+2;
  if (savestrlen+nadd > savestrmax) {
    savestrmax += (nadd > SAVESTRUNIT) ? nadd : SAVESTRUNIT;
    buf=malloc(savestrmax*sizeof(char));
    if (buf == NULL) {
      *status = *status | newstat | CPSAVE_OUTOFMEMORY;
      return retval;
    }
    else {
      strcpy(buf,savestr);
      free(savestr);
      savestr=buf;
    }
  }

  strcat(savestr,name);
  strcat(savestr,"=");
  strcat(savestr,strval);
  strcat(savestr,PARMSEPRTR);
  savestrlen+=nadd-1;

  *status = *status | newstat;
  return retval;
}


double cmdparams_save_time (CmdParams_t *parms, char *name, int *status) {

  double retval;
  char *strval, *buf;
  int nadd, stat;
  int newstat=0;

  retval=cmdparams_get_time (parms, name, &stat);
  switch (stat) {
  case CMDPARAMS_SUCCESS:
    newstat=CPSAVE_SUCCESS;
    break;
  case CMDPARAMS_UNKNOWN_PARAM:
    newstat=CPSAVE_UNKNOWN_PARAM;
    break;
  case CMDPARAMS_INVALID_CONVERSION:
    newstat=CPSAVE_INVALID_CONVERSION;
    break;
  case CMDPARAMS_OUTOFMEMORY:
    newstat=CPSAVE_OUTOFMEMORY;
    break;
  default:
    newstat=CPSAVE_UNKNOWN_ERROR;
  }

  if (stat == CMDPARAMS_UNKNOWN_PARAM) {
    *status = *status | newstat;
    return retval;
  }

  if (savestrmax == 0) {
    savestr=calloc(SAVESTRUNIT,sizeof(char));
    if (savestr == NULL) {
      *status = *status | newstat | CPSAVE_OUTOFMEMORY;
      return retval;
    }
    savestrmax=SAVESTRUNIT;
  }

  strval=(char *)cmdparams_get_str (parms, name, NULL);
  nadd=strlen(strval)+strlen(name)+strlen(PARMSEPRTR)+2;
  if (savestrlen+nadd > savestrmax) {
    savestrmax += (nadd > SAVESTRUNIT) ? nadd : SAVESTRUNIT;
    buf=malloc(savestrmax*sizeof(char));
    if (buf == NULL) {
      *status = *status | newstat | CPSAVE_OUTOFMEMORY;
      return retval;
    }
    else {
      strcpy(buf,savestr);
      free(savestr);
      savestr=buf;
    }
  }

  strcat(savestr,name);
  strcat(savestr,"=");
  strcat(savestr,strval);
  strcat(savestr,PARMSEPRTR);
  savestrlen+=nadd-1;

  *status = *status | newstat;
  return retval;
}


int cmdparams_save_flag (CmdParams_t *parms, char *name, int *status) {

  int retval;
  char *strval, *buf;
  int nadd, stat;
  int newstat=0;

  if (cmdparams_exists (parms, name)) {
     retval = cmdparams_get_int (parms, name, &stat);
  } 
  else {
    stat=CMDPARAMS_SUCCESS;
    retval=0;
  }

  switch (stat) {
  case CMDPARAMS_SUCCESS:
    newstat=CPSAVE_SUCCESS;
    break;
  case CMDPARAMS_UNKNOWN_PARAM:
    newstat=CPSAVE_UNKNOWN_PARAM;
    break;
  case CMDPARAMS_INVALID_CONVERSION:
    newstat=CPSAVE_INVALID_CONVERSION;
    break;
  case CMDPARAMS_OUTOFMEMORY:
    newstat=CPSAVE_OUTOFMEMORY;
    break;
  default:
    newstat=CPSAVE_UNKNOWN_ERROR;
  }

  if (savestrmax == 0) {
    savestr=calloc(SAVESTRUNIT,sizeof(char));
    if (savestr == NULL) {
      *status = *status | newstat | CPSAVE_OUTOFMEMORY;
      return retval;
    }
    savestrmax=SAVESTRUNIT;
  }

  strval=(char *)cmdparams_get_str (parms, name, NULL);
  nadd=strlen(strval)+strlen(name)+strlen(PARMSEPRTR)+2;
  if (savestrlen+nadd > savestrmax) {
    savestrmax += (nadd > SAVESTRUNIT) ? nadd : SAVESTRUNIT;
    buf=malloc(savestrmax*sizeof(char));
    if (buf == NULL) {
      *status = *status | newstat | CPSAVE_OUTOFMEMORY;
      return retval;
    }
    else {
      strcpy(buf,savestr);
      free(savestr);
      savestr=buf;
    }
  }

  strcat(savestr,name);
  strcat(savestr,"=");
  strcat(savestr,strval);
  strcat(savestr,PARMSEPRTR);
  savestrlen+=nadd-1;

  *status = *status | newstat;
  return retval;
}


char *	cmdparams_save_arg (CmdParams_t *parms, int num, int *status) {

  char *strval, *istr, *buf;
  int nadd, i;
  static int savenum;

  strval=(char *)cmdparams_getarg(parms,num);
  if (strval == NULL || num <= savenum) return strval;

  if (savestrmax == 0) {
    savestr=calloc(SAVESTRUNIT,sizeof(char));
    if (savestr == NULL) {
      *status = *status | CPSAVE_OUTOFMEMORY;
      return strval;
    }
    savestrmax=SAVESTRUNIT;
  }

  for (i=savenum+1; i<=num; i++) {
    istr=(char *)cmdparams_getarg(parms,i);
    nadd=strlen(istr)+strlen(PARMSEPRTR)+1;
    if (savestrlen+nadd > savestrmax) {
      savestrmax += (nadd > SAVESTRUNIT) ? nadd : SAVESTRUNIT;
      buf=malloc(savestrmax*sizeof(char));
      if (buf == NULL) {
        *status = *status | CPSAVE_OUTOFMEMORY;
        break;
      }
      else {
        strcpy(buf,savestr);
        free(savestr);
        savestr=buf;
      }
    }

    strcat(savestr,istr);
    strcat(savestr,PARMSEPRTR);
    savestrlen+=nadd-1;
  }
  savenum=num;
  return strval;
}


char * cmdparams_save_str (CmdParams_t *parms, char *name, int *status) {

  char *strval;
  int nadd, stat;
  int newstat=0;
  char *buf;

  strval=(char *)cmdparams_get_str (parms, name, &stat);
  switch (stat) {
  case CMDPARAMS_SUCCESS:
    newstat=CPSAVE_SUCCESS;
    break;
  case CMDPARAMS_UNKNOWN_PARAM:
    newstat=CPSAVE_UNKNOWN_PARAM;
    break;
  case CMDPARAMS_INVALID_CONVERSION:
    newstat=CPSAVE_INVALID_CONVERSION;
    break;
  case CMDPARAMS_OUTOFMEMORY:
    newstat=CPSAVE_OUTOFMEMORY;
    break;
  default:
    newstat=CPSAVE_UNKNOWN_ERROR;
  }

  if (stat == CMDPARAMS_UNKNOWN_PARAM) {
    *status = *status | newstat;
    return strval;
  }

  if (savestrmax == 0) {
    savestr=calloc(SAVESTRUNIT,sizeof(char));
    if (savestr == NULL) {
      *status = *status | newstat | CPSAVE_OUTOFMEMORY;
      return strval;
    }
    savestrmax=SAVESTRUNIT;
  }

  nadd=strlen(strval)+strlen(name)+strlen(PARMSEPRTR)+4;
  if (savestrlen+nadd > savestrmax) {
    savestrmax += (nadd > SAVESTRUNIT) ? nadd : SAVESTRUNIT;
    buf=malloc(savestrmax*sizeof(char));
    if (buf == NULL) {
      *status = *status | newstat | CPSAVE_OUTOFMEMORY;
      return strval;
    }
    else {
      strcpy(buf,savestr);
      free(savestr);
      savestr=buf;
    }
  }

  strcat(savestr,name);
  strcat(savestr,"=\"");
  strcat(savestr,strval);
  strcat(savestr,"\"");
  strcat(savestr,PARMSEPRTR);
  savestrlen+=nadd-1;

  *status = *status | newstat;
  return strval;
}
