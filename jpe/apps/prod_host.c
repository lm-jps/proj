static char rcsid[] = "$Header: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/jpe/apps/prod_host.c,v 1.1 2009/02/03 22:18:47 production Exp $";
/* prod_host.c
 *
 * This will set the static variables according to the info
 * in the file /home/soi/CM/tables/prod_host.cfg
 * This is called once at the initialization of a program like
 * ampex_svc or dsds_svc, etc..
 *
*/

#include <stdlib.h>
#include <stdio.h>
#include "pe.h"

/* hide these vrbls */
static char prod_prime[256];
static char prod_second[256];
static char prod_prime_short[256];
static char prod_second_short[256];

/* Set the variables from the file. Return 1 on success, else 0
 * The file typically looks like this:
 * PROD_PRIME=sonar
 * PROD_SECOND=tarax
*/
int prod_host_set() {
  FILE *fp;
  char line[256], ltmp[256];
  char *cptr;

  if(!(fp=fopen(PHNAME, "r"))) {
    return(0);
  }
  while(fgets(line, 256, fp)) {
    if(line[0] == '#' || line[0] == '\n') continue;
    if(cptr = (char *)rindex(line, '\n')) *cptr = NULL; /* elim term. CR */
    if(strstr(line, "PROD_PRIME")) {
      if(!(cptr = (char *)index(line, '='))) {
        fclose(fp);
        return(0);			/* file has bad format */
      }
      else {
        strcpy(prod_prime, cptr+1);
        strcpy(ltmp, prod_prime);
        cptr = strchr(ltmp, '.');
        if(cptr) { *cptr = NULL; };
        strcpy(prod_prime_short, ltmp);
      }
    }
    if(strstr(line, "PROD_SECOND")) {
      if(!(cptr = (char *)index(line, '='))) {
        fclose(fp);
        return(0);			/* file has bad format */
      }
      else {
        strcpy(prod_second, cptr+1);
        strcpy(ltmp, prod_second);
        cptr = strchr(ltmp, '.');
        if(cptr) { *cptr = NULL; };
        strcpy(prod_second_short, ltmp);
      }
    }
  }
  fclose(fp);
  return(1);
}

char *prod_host_prime() {
  return(prod_prime);
}
char *prod_host_second() {
  return(prod_second);
}
char *prod_host_prime_short() {
  return(prod_prime_short);
}
char *prod_host_second_short() {
  return(prod_second_short);
}

/*
$Id: prod_host.c,v 1.1 2009/02/03 22:18:47 production Exp $
$Source: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/jpe/apps/prod_host.c,v $
$Author: production $
*/
/* $Log: prod_host.c,v $
 * Revision 1.1  2009/02/03 22:18:47  production
 * initial
 *
 * Revision 1.3  2001/10/22  21:59:29  jim
 * add prod_host_prime_short and prod_host_second_short
 *
 * Revision 1.2  2001/08/24  22:03:51  jim
 * initial
 * */
