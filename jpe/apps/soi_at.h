#ifndef AT_OK
#include <timeio.h>

#ifndef SOI_VERSION_INCL
#include "soi_version.h"
#endif

typedef struct at_struct {
	char *pathname;		/* current file name */
	char *delim;		/* field delimiter, default '\t' */
	char *tokdelim;		/* strtok match delim */
	int rows;		/* number of rows of data present */
	int row;		/* current row number for reading */
	int maxrows;		/* max rows with space allocated */
	int cols;		/* number of cols of data present */
	int maxcols;		/* max cols with space allocated */
	int header_done;	/* flag showing if header is read */
	char **col_name;	/* pointer to array of strings */
	char ***col_vals;	/* pointer to array of arrays of strings */
	int *col_width;		/* pointer to array of field widths */
	int last_error;		/* error kind */
	char *last_called;	/* name of last user level function called */
} AT;
								   /*  at.c  */
/* function prototypes, I/O */

AT *at_read (char *pathname);
int at_write (AT *at, char *pathname);

/* function prototypes, create and insert info */

AT *at_new ();
AT *at_create (char *header_info);
AT *at_define_table (int cols, int rows);
int at_put_info (AT *at, char *header_info);
int at_put_value (AT *at, int col, int row, double value);
int at_put_time (AT *at, int col, int row, double value);
int at_put_string (AT *at, int col, int row, char *text);
int at_put_name (AT *at, int col, char *name);
char *at_error_msg (AT *at);

/* function prototypes, retrieve info */

int at_cols (AT *at);
int at_rows (AT *at);
int at_col (AT *at, char *data_name);
char *at_get_name (AT *at, int col);
double at_get_value (AT *at, int col, int row);
double at_get_time (AT *at, int col, int row);
char *at_get_string (AT *at, int col, int row);
char *at_get_ascii (AT *at, int col, int row);
int at_in_bounds (AT *at, int col, int row);
int at_check_n_grow (AT *at, int col, int row);
void at_free (AT *at);
							    /*  at_setkey.c  */
/* getkey definitions for AT library */

char *at_rec_getkey_str (AT *at, int rec, char *key);
char *AT_rec_getkey_str (AT *at, int rec, char *key);
int at_rec_getkey_int (AT *at, int rec, char *key);
unsigned int at_rec_getkey_uint (AT *at, int rec, char *key);
double at_rec_getkey_double (AT *at, int rec, char *key);
TIME at_rec_getkey_time (AT *at, int rec, char *key);
TIME at_rec_getkey_time_interval (AT *at, int rec, char *key);

/* setkey definitions for AT library */

int at_rec_setkey_str (AT *at, int rec, char *key, char *val);
int at_rec_setkey_int (AT *at, int rec, char *key, int ival);
int at_rec_setkey_uint (AT *at, int rec, char *key, unsigned int ival);
int at_rec_setkey_double (AT *at, int rec, char *key, double dval);
int at_rec_setkey_time (AT *at, int rec, char *key, TIME tval);
int at_rec_setkey_time_interval (AT *at, int rec, char *key, TIME tival);
							      /*  at_util.c  */
int at_copy_row (AT *at, int new, int old);
int at_remove_row (AT *at, int row);
int at_remove_rows (AT *at, int row0, int row1);

#define AT_EMPTY_LINE		(-7)
#define AT_BAD_NAMELIST		(-6)
#define AT_MALLOC_ERROR		(-5)
#define AT_BAD_DATA_VALUE	(-4)
#define AT_BAD_TYPE		(-3)
#define AT_BAD_REC_NUMBER 	(-2)
#define AT_BAD_COL_COUNT 	(-1)
#define AT_BAD			(0)
#define AT_OK 			(1)
#define AT_MIN_ROW_ALLOC	(256)
#define AT_MIN_COL_ALLOC	(16)
#define AT_MAX_LINE		(10000)
#define AT_DELIM		"\t"
#define AT_TOKDELIM		"\n\t"

#endif

/*
$Id: soi_at.h,v 1.1 2009/04/24 21:52:33 production Exp $
$Source: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/jpe/apps/soi_at.h,v $
$Author: production $
 * $Log: soi_at.h,v $
 * Revision 1.1  2009/04/24 21:52:33  production
 * *** empty log message ***
 *
 * Revision 1.16  2007/05/04  21:12:11  rick
 * changed soi_time.h to timeio.h
 *
 * Revision 1.15  2006/06/14 00:21:21  rick
 * added declaration for at_new()
 *
 * Revision 1.14  2001/07/06 21:25:55  rick
 * added prototypes from at_util.c
 *
 * Revision 1.13  1998/04/13  23:52:42  kay
 * added at_rec_setkey_uint and at_rec_getkey_uint
 *
 * Revision 1.12  1997/10/08  01:06:36  phil
 * *** empty log message ***
 *
 * Revision 1.11  1997/04/16  21:50:32  kehcheng
 * added #include <soi_version.h>
 *
 * Revision 1.10  1996/10/15  05:45:56  phil
 * added at_rec setkey/getkey functions
 *
 * Revision 1.9  1996/05/10  18:47:23  phil
 * changed AT_MIN_ALLOC to AT_MIN_COL_ALLOC and AT_MIN_ROW_ALLOC
 *
 * Revision 1.8  1996/03/11  21:43:17  jim
 * noop
 *
 * Revision 1.7  1996/01/10  02:44:08  phil
 * Removed AT_MISSING and MISSING definitions since the only place
 * AT_MISSING was used was changed to D_MISSING.
 *
 * Revision 1.6  1995/03/15  02:09:07  phil
 * increased AT_MAXLINE
 *
 * Revision 1.5  1995/02/13  20:17:48  phil
 * *** empty log message ***
 *
 * Revision 1.4  1995/02/01  21:47:41  phil
 * *** empty log message ***
 *
 * Revision 1.3  1994/09/28  16:09:22  jim
 * added at_in_bounds() prototype
 *
 * Revision 1.2  1994/09/27  18:23:05  CM
 * fixed at_create() prototype
 *
 * Revision 1.1  1994/09/27  18:01:51  CM
 * Initial revision
 * */
