/*
 *  soi_missing.h				~soi/CM/include/soi_missing.h
 *
 *  Definitions of standard missing or fill values for SOI datasets
 *  and programs.
 *  Additional information is in the following man pages:
 *	soi_NaN.h
 *	SOI TechNote 115
 *
 *
 *  Bugs:
 *
 *  Revision history is at the end of the file.
 */
#ifndef SOI_MISSING_INCL
/****************************************************************************/
/**************************  INCLUDE STATEMENTS  ****************************/
/****************************************************************************/

#ifndef SOI_VERSION_INCL
#include <soi_version.h>
#endif

/****************************************************************************/
/****************************  DEFINE STATEMENTS  ***************************/
/****************************************************************************/

#define SOI_MISSING_VERSION_NUM	(0.9)
#define SOI_MISSING_INCL	1

#define B_MISSING	(-128)
#define S_MISSING	(-32768)
#define I_MISSING	(-2147483647-1)
#define UB_MISSING	(255)
#define US_MISSING	(65535)
#define UI_MISSING	(2147483647)
#define F_MISSING	(A_Quiet_fNaN())
#define D_MISSING	(A_Quiet_dNaN())

#define T_MISSING	(-211087684800.0)
#define T_MISSING_STR	("-4712.01.01_12:00:00.000_UT")

#ifndef MISSING
#define MISSING		(-8388608.0e10)
#endif

/****************************************************************************/
/*******************************  TYPEDEFS  *********************************/
/****************************************************************************/
  
/****************************************************************************/
/*********************  GLOBAL & EXTERN DECLARATIONS  ***********************/
/****************************************************************************/

/****************************************************************************/
/****************************  MACRO DEFINITIONS  ***************************/
/****************************************************************************/

#define is_B_MISSING(v)	(v == B_MISSING)
#define is_S_MISSING(v)	(v == S_MISSING)
#define is_I_MISSING(v) (v == I_MISSING)
#define is_UB_MISSING(v) (v == UB_MISSING)
#define is_US_MISSING(v) (v == US_MISSING)
#define is_UI_MISSING(v) (v == UI_MISSING)
#define is_F_MISSING(v)	(IsfNaN(v))
#define is_D_MISSING(v) (IsdNaN(v))
#define is_C_MISSING(v) (IsfNaN((v).r) || IsfNaN((v).i))

#define is_T_MISSING(v) (v == T_MISSING)

/****************************************************************************/
/***************************  FUNCTION PROTOTYPES  **************************/
/****************************************************************************/

#endif
/*
 *  Revision History
 */

/*
 * Revision 1.1  1994/02/16  23:21:14  CM
 * Initial revision
 * */

/*
$Id: soi_missing.h,v 1.1 2009/04/24 21:53:02 production Exp $
$Source: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/jpe/apps/soi_missing.h,v $
$Author: production $
*/
/* $Log: soi_missing.h,v $
 * Revision 1.1  2009/04/24 21:53:02  production
 * *** empty log message ***
 *
 * Revision 1.7  1997/12/15  18:18:12  phil
 * added is_T_MISSING
 *
 * Revision 1.6  1997/04/16  21:54:02  kehcheng
 * added #include <soi_version.h>
 *
 * Revision 1.5  1995/08/29  18:15:21  CM
 * auto rcsfix by CM
 * */
