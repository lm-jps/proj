/*
 *  soi_args.h						~CM/include/soi_args.h
 *
 *  Argument declarations for modules
 *
 *  Responsible:  Rick Bogart			RBogart@solar.Stanford.EDU
 *
 *  Bugs:
 *    There are no unsigned fixed-point types for arguments
 *
 *  Revision history is at the end of the file.
 */

#ifndef __ARGS_DECLARED__

/****************************************************************************/
/****************************  DEFINE STATEMENTS  ***************************/
/****************************************************************************/

#define __ARGS_DECLARED__
					     /*  kinds of module arguments  */
#define ARG_END		(0)
     /*  1 formerly assigned to ARG_DATASET; replaced with ARG_DATA_IN/OUT  */
#define ARG_OBSOLETE	(1)
#define ARG_DATASET     (1)     /* replaced with ARG_DATA_IN/OUT */
#define ARG_FLAG	(2)
#define ARG_TIME	(3)
#define ARG_INT		(4)
#define ARG_FLOAT	(5)
#define ARG_STRING	(6)
#define ARG_FILEPTR	(7)
#define ARG_DATA_IN	(8)
#define ARG_DATA_OUT	(9)
#define ARG_INTS	(14)
#define ARG_FLOATS	(15)
#define ARG_NUME	(24)

/****************************************************************************/
/*******************************  TYPEDEFS  *********************************/
/****************************************************************************/
  
 /*  form of entry in table of arguments needed by strategy level function  */
typedef struct ARG {
   int   kind;
   char *key;
   char *default_value;
   char *range;
   char *description;
} argument;

#endif

/*
 *  Revision History
 *  01.07.24	Rick Bogart	created this file; removed the declarations
 *			from module.h in order to make them accessible to
 *			the sds library
 */
