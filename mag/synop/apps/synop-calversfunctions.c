char *cvsinfo = "cvsinfo: $Header: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/mag/synop/apps/synop-calversfunctions.c,v 1.1 2016/05/05 19:49:28 arta Exp $";

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "jsoc_main.h"


/* 
the following function is taken directly from chapter 2 of Kernighan and Ritchie
getbits: get n bits from position p
*/
unsigned long long getbits(unsigned long long x, int p, int n)
{
  return (x >> (p+1-n)) & ~(~0ull << n);
}

/*
Richard Heathfield's solution to exercise 2-6 of K&R 
*/
unsigned long long setbits(unsigned long long x, int p, int n, unsigned long long y)
{
  return (x & ((~0ull << (p + 1)) | (~(~0ull << (p + 1 - n))))) | ((y & ~(~0ull << n)) << (p + 1 - n));
// another solution (untested):
//return (x & ~(~(~0<<n)<<(p+1)) | (y&~(~0<<n)) << (p+1));
}

unsigned long long fixcalver64(unsigned long long x)
{
  int firstnibble = getbits(x,3,4);
  if (firstnibble == 0)
    return setbits(x,3,4,2);
  else
    return x;
}


