/*
 *  module.h					~soi/(version)/include/module.h
 *
 *  Data structures used in 4-level model for SSSC compute modules, described
 *    in "Programming in the SOI Analysis Environment", SOI TN-93-107.
 *  Levels are:  use (interface), strategy, function, and library (support).
 *
 *  Responsible:  Kay Leibrand			KLeibrand@solar.Stanford.EDU
 *
 *  Bugs:
 *    There are no unsigned fixed-point types for arguments
 *
 *  Revision history is at the end of the file.
 */
#ifndef __MODULE_TYPES__

/****************************************************************************/
/**************************  INCLUDE STATEMENTS  ****************************/
/****************************************************************************/

#include "soi_version.h"
#include "soi_args.h"
#include "soi_sds.h"
//#include "soi_fun.h"
#include "soi_names.h"
#include "soi_error.h"
#include "timeio.h"
#include "soi_vds.h"

/****************************************************************************/
/****************************  DEFINE STATEMENTS  ***************************/
/****************************************************************************/

#define MODULE_ERROR	(-1)
#define MODULE_OK	(0)

#define MODULES_VERSION_NUM	(4.8)
#define __MODULE_TYPES__

/****************************************************************************/
/*********************  GLOBAL & EXTERN DECLARATIONS  ***********************/
/****************************************************************************/

//extern int (*DoIt)();  //old. now in jsoc_main.h
		  /*  table of arguments needed by strategy level function  */
extern  argument arguments[];
		 /*  table of arguments created by strategy level function  */
extern  argument outarguments[];

#endif
/*
 *  Revision History
 *  V 0.0   93.02.26	Kay Leibrand	created this file			
 *	    93.04.17	Rick Bogart	include soi_sds and soi_fun 
 *	    93.06.07	K Leibrand	added outarguments
 *  V 0.7   93.08.13	R Bogart	include soi_names &
 *		soi_error; soi_key included by soi_names
 *	    93.09.10	R Bogart	defined version number
 *	    94.02.07	K Leibrand	modified argument struct to
 *		put kind first, added DATASET types for IN and OUT,
 *  V 0.8   94.02.08    K Leibrand	defined version number
 *	    94.09.02	Jim Aloise	Replace ARG_DATASET with
 *		ARG_DATA_IN and ARG_DATA_OUT
 *  v 0.9   94.11.14	R Bogart	V 0.9
 *  v 1.0   95.01.30	R Bogart	include soi_time.h; updated version #
 *	    95.09.25	R Bogart	include soi_vds.h
 *  v 4.5   99.11.11	R Bogart	added arg types for pointers
 *	    99.12.13	R Bogart	changed ARG struct member names from
 *		minvalid & maxvalid to range and description
 *  v 4.8   00.06.16	R Bogart	added arg type nume for menus (enum)
 *	    01.07.24	R Bogart	removed argument declarations to new
 *		include files soi_args.h	
 */

/*
$Id: module.h,v 1.1 2009/04/24 21:50:44 production Exp $
$Source: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/jpe/apps/module.h,v $
$Author: production $
 * $Log: module.h,v $
 * Revision 1.1  2009/04/24 21:50:44  production
 * *** empty log message ***
 *
 * Revision 1.13  2007/05/04  21:13:19  rick
 * changed soi_time.h to timeio.h
 *
 * Revision 1.12  2001/07/24 23:15:57  rick
 * see above
 *
 * Revision 1.11  2001/04/11  19:02:41  rick
 * see above
 *
 * Revision 1.10  1999/12/14  00:12:22  rick
 * see above
 *
 * Revision 1.9  1999/11/18  01:03:43  rick
 * see above
 *
 * Revision 1.8  1997/04/16  21:47:07  kehcheng
 * added #include <soi_version.h>
 *
 * Revision 1.7  1995/09/25  15:38:03  rick
 * include soi_vds.h
 *
 * Revision 1.6  1995/08/28  18:38:08  kay
 * no change
 *
 * Revision 1.5  1995/01/30  22:41:31  rick
 * included soi_time.h
 *
 * Revision 1.4  1994/11/14  23:12:02  rick
 * version 0.9
 *
 * Revision 1.3  1994/09/12  22:59:27  jim
 * added ARG_DATA_IN and OUT
 *
 * Revision 1.2  1994/09/02  17:50:39  jim
 * Replace ARG_DATASET with ARG_DATA_IN and ARG_DATA_OUT
 *
 * Revision 1.1  1994/02/16  23:21:14  CM
 * Initial revision
 * */
