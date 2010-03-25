#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <strings.h>
#include <sys/types.h>
#include <sys/time.h>


static void sprint_time_ISO (char *tstring, double t)
{
sprint_time(tstring,t,"UTC",0);
tstring[4] = tstring[7] = '-';
tstring[10] = 'T';
tstring[19] = '\0';
}

int main(int argc, char *argv[]) {
  char string[128];
  double t;

  t = strtod(argv[1], NULL);
  sprint_time(string, t, "UTC", 0);
  printf("%s\n", string);
}
