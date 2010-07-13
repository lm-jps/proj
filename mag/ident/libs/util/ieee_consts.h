/*
 * GENERATING IEEE FLOATING-POINT CONSTANTS
 * (ieee_consts.c)
 *
 * Michael Turmon, 2002
 * 
 */

#ifndef _ieee_consts_h_
#define _ieee_consts_h_

#ifdef __cplusplus
extern "C" {
#endif

float  mxt_getnanf(void);
double mxt_getnand(void);
double mxt_getinfd(void);

#ifdef __cplusplus
}	/* extern "C" */
#endif

#endif /* _ieee_consts_h_ */
