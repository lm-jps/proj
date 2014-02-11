/*
$Header: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/mag/synop/apps/dataset.h,v 1.1 2014/02/11 18:30:19 arta Exp $
$Log: dataset.h,v $
Revision 1.1  2014/02/11 18:30:19  arta
Finish removing dependencies on Yangs home directory. Added several headers needed for the synop project, but that resided in Yangs myproj directory.

 * Revision 1.3  90/12/18  15:52:07  phil
 * local.h included
 * 
 * Revision 1.2  90/10/05  12:30:22  phil
 * remove kw macros
 * 
 * Revision 1.1  90/09/17  14:32:02  phil
 * Initial revision
 * 
*/
#ifndef DATASETINCLUDED
#define DATASETINCLUDED

#include <local.h>
#include <stdtime.h>
#include <stdio.h>
#include <constant.h>
#include <strings.h>
#include <math.h>
#include <param.h>
#include <keywords.h>
#include <nametables.h>

#include <sys/types.h>
#include <unistd.h>
#include <errno.h>

#ifndef __linux__
extern int errno;
#endif

#define ds_is(such,ds) (ds->ds_access&such)
#define ds_return(r) \
		{\
		ds->ds_blocktrap=FALSE;\
		if(ds->ds_pendingtrap)\
			ds_die();\
		return(r);\
		}

#define TRUE	1
#define FALSE	0
#define MAXT	16
#define MAXPATH	256

	/* ds_rbuf flags */
#define ds_USED 	'\01'
#define ds_EMPTY	'\0'

#define ds_allmissing(ds) (		\
	ds->ds_now >= ds->ds_last  ||	\
	ds->ds_now <  ds->ds_first ||	\
	(ds->ds_gaps && *ds->ds_rbuf == ds_EMPTY) )

	/* ds_state flags */
#define AFTEROPEN	1
#define IO_PENDING	2
#define IO_OK		3
#define SEEKERROR	4
#define AFTERCLOSE	5
#define CLOSED		6
	/* ds_access flags */
#define RD	001
#define WRT	002
#define RAND	004
#define APND	010
#define PIPE	020
#define BLK	040
#define FOURK	0100
	/* ds_kind flags */
#define ds_LIST		1
#define ds_TIME		2
#define ds_CARR		3
#define ds_FREQ		4

	/* ds_erract flags */
#define ds_DIE		1
 
	/* keyword list support */
#define	ds_KEYBEGIN	"### KEYWORD BEGIN ###\n"
#define	ds_KEYEND	"### KEYWORD END ###\n"

	/* d_type flags */
#define d_CHAR		7
#define	d_BYTE		6
#define	d_SHORT		1
#define	d_LONG		2
#define	d_FLOAT		3
#define	d_DOUBLE	4
#define	d_TIME		5

extern short  d_tsize[8];

	/* Machine number representation flags */
#define	ds_f_VAX	0
#define	ds_f_IEEE	16

	/* Machine architecture flag words */
#define	ds_VAX		( ds_f_VAX )
#define	ds_DECMIPS	( ds_f_IEEE )
#define	ds_SGIMIPS	( ds_f_IEEE )
#define	ds_NOMACHINE	(0177777)

	/* set machine type */

/* ALL datasets will be written little-endian so conversions need to be done on
 * input and output to maintain external little-endian standard. Also, no checking
 * needs to be done on the input side. */

#ifdef	vax
#define	ds_MACHINE	ds_VAX

#define ds_INDS_ID(id)	/* no-op */
#define ds_INDOUBLE(d)	{ if (ds->ds_id & ds_f_IEEE) IEEE_to_VAX_d(d); }
#define ds_INFLOAT(d)	{ if (ds->ds_id & ds_f_IEEE) IEEE_to_VAX_f(d); }
#define ds_INLONG(i)	/* no-op */
#define ds_INSHORT(i)	/* no-op */
#define ds_OUTDOUBLE(d)	{ if (ds->ds_id & ds_f_IEEE) VAX_to_IEEE_d(d); }
#define ds_OUTFLOAT(d)	{ if (ds->ds_id & ds_f_IEEE) VAX_to_IEEE_f(d); }
#define ds_OUTLONG(i)	/* no-op */
#define ds_OUTSHORT(i)	/* no-op */
#endif

/* HACK */
#ifdef __linux
#define MIPSEL
#endif

#ifdef	MIPSEL
#define ds_MACHINE	ds_DECMIPS
#define ds_INDS_ID(id)	/* no -o */
#define ds_INDOUBLE(d)	{ if (!(ds->ds_id & ds_f_IEEE)) VAX_to_IEEE_d(d); }
#define ds_INFLOAT(d)	{ if (!(ds->ds_id & ds_f_IEEE)) VAX_to_IEEE_f(d); }
#define ds_INLONG(i)	/* no -o */
#define ds_INSHORT(i)	/* no -o */
#define ds_OUTDOUBLE(d)	{ if (!(ds->ds_id & ds_f_IEEE)) IEEE_to_VAX_d(d); }
#define ds_OUTFLOAT(d)	{ if (!(ds->ds_id & ds_f_IEEE)) IEEE_to_VAX_f(d); }
#define ds_OUTLONG(i)	/* no -o */
#define ds_OUTSHORT(i)	/* no -o */
#endif

/* HACK */
#ifdef __sun
#define MIPSEB
#endif

#ifdef	MIPSEB
#define ds_MACHINE	ds_SGIMIPS
#define ds_INDS_ID(id)	{ gen_swap(id,1,2); }
#define ds_INDOUBLE(d)	{						\
			if (ds->ds_id & ds_f_IEEE)			\
				gen_swap(d,1,sizeof(double));		\
			else						\
				{					\
				VAX_to_IEEE_d(d);			\
				}					\
			}
#define ds_INFLOAT(d)	{						\
			if (ds->ds_id & ds_f_IEEE)			\
				gen_swap(d,1,sizeof(float));		\
			else						\
				{					\
				VAX_to_IEEE_f(d);			\
				}					\
			}
#define ds_INLONG(i)	{ gen_swap(i,1,sizeof(long)); }
#define ds_INSHORT(i)	{ gen_swap(i,1,sizeof(short)); }
#define ds_OUTDOUBLE(d)	{						\
			if (ds->ds_id & ds_f_IEEE)			\
				gen_swap(d,1,sizeof(double));		\
			else						\
				{					\
				IEEE_to_VAX_d(d);			\
				}					\
			}
#define ds_OUTFLOAT(d)	{						\
			if (ds->ds_id & ds_f_IEEE)			\
				gen_swap(d,1,sizeof(float));		\
			else						\
				{					\
				IEEE_to_VAX_f(d);			\
				}					\
			}
#define ds_OUTLONG(i)	{ gen_swap(i,1,sizeof(long)); }
#define ds_OUTSHORT(i)	{ gen_swap(i,1,sizeof(short)); }
#endif

#ifndef	ds_MACHINE
#define	ds_MACHINE	ds_NOMACHINE
#define ds_INDOUBLE(d)	/* no-op */
#define ds_INFLOAT(d)	/* no-op */
#define ds_INLONG(i)	/* no-op */
#define ds_INSHORT(i)	/* no-op */
#define ds_OUTDOUBLE(d)	/* no-op */
#define ds_OUTFLOAT(d)	/* no-op */
#define ds_OUTLONG(i)	/* no-op */
#define ds_OUTSHORT(i)	/* no-op */
#endif

typedef struct 
	{
	short	d_type;		/* data type				*/
	short	d_xwidth;	/* x dimension of data array		*/
	short	d_ywidth;	/* y dimension of data array		*/
	short	d_off_old;	/* old byte offset within group 2/5/87	*/
	char	d_layout[MAXT];	/* format for printf			*/
	char	d_name[MAXT];	/* descriptive name			*/
	char	d_unit[MAXT];	/* unit of measure i.e. "m/s"		*/
	short	d_chars;	/* number of character places in format */
	short	d_valid;	/* true if stats accurately reflect ds	*/
	double	d_scale;
	double	d_base;		/* data = base + scale*(stored data)	*/
	long	d_n;		/* number of non-missing values		*/
	double	d_min;		/* minimum value of data present	*/
	double	d_max;		/* maximum value			*/
	double	d_mean;		/* mean value				*/
	double	d_stddev;	/* standard deviation of data		*/
	double	d_miss;		/* value for missing data		*/
	long	d_off;		/* byte offset within group 2/5/87	*/
	long	d_l1; 		/* extra				*/
	double	d_d2;		/* extra				*/
	} DS_DESC;
#define ds_desc DS_DESC

#define ds_HEADSIZE	84	/*size of recorded part of dataset head	*/

typedef struct
	{
	char	ds_name[MAXT];	/* descriptive name			*/
	short	ds_numitem;	/* number of items per group		*/
	short	ds_gaps;	/* gaps in file flag. (see lseek)	*/
	char	ds_tname[MAXT];	/* "time" inc descriptive name		*/
	TIME	ds_first;	/* first time of data			*/
	TIME	ds_last;	/* last time in dataset + 1		*/
	TIME	ds_inc;		/* time increment between groups	*/
	short	ds_kind;	/* "time" kind				*/
	short	ds_id;		/* machine architecture			*/
	long	ds_txtsz;	/* size of text				*/
	double	ds_d1;		/* extra				*/
	double	ds_d2;		/* extra				*/
/* preceding strcture is equivalent to a dataset header			*/
	ds_desc	**ds_item;	/* array of pointers to descriptors	*/
	TIME	ds_now;		/* user present time			*/
	long	ds_ngrps;	/* number of groups in dataset		*/
	FILE 	*ds_fp;		/* stream pointer			*/
	long 	ds_fd;		/* file descriptor for ds_fp		*/
	short	ds_access;	/* opendata access flag			*/
	char    *ds_fullpath;   /* pathname used to open dataset	*/
	short	ds_state;	/* dataset state			*/
	TIME	ds_clk;		/* internal present time		*/
	long	ds_offset;	/* present location within dataset	*/
	long	ds_begin;	/* offset of first datagroup		*/
	long	ds_end;		/* offset of "last" datagroup		*/
	long	ds_sizeof;	/* group size + 2 if gaps == TRUE	*/
	char	*ds_rbuf;	/* data group i/o buffer		*/
	char	*ds_group;	/* data group pointer			*/
	char	*ds_mgrp;	/* pointer to all missing group		*/
	short	ds_erract;	/* code indicates action on error	*/
	FILE	*ds_txtfp;	/* file pointer for descriptive text	*/
	char	ds_txtnm[20];	/* file name for descriptive text	*/
	KEYWORDS *ds_keywords;	/* keyword list pointer			*/
	FILE	*ds_tmpfp;	/* file descriptor for plumbing		*/
	char	ds_tmpnm[20];	/* file name for plumbing		*/
	char	ds_firststr[80];/* working string for first time	*/
	char	ds_laststr[80];	/* working string for last time		*/
	long	ds_thisds;	/* index of all open datasets		*/
	long	ds_blocktrap;	/* flag to block user traps		*/
	long	ds_pendingtrap;	/* flag that trap pending by user	*/
	} DATASET;

#define dataset DATASET 

typedef struct
	{
	long item;
	long takeall;
	long xind;
	long yind;
	} SELECTLIST;

#define dsdefine ds_def
double	dget(DATASET *ds, int sel, int ix, int iy);
int	dget_i(DATASET *ds, int sel, int ix, int iy);
double	dsgetind();
TIME	dstimeparam(DATASET *ds, char *parm, TIME dflt);
TIME	dsincparam(DATASET *ds, char *parm, TIME dflt);
double	*dgetitem(DATASET *ds, int sel, double *buf);
double	*dgetitem_d(DATASET *ds, int sel, double *buf);
float	*dgetitem_f(DATASET *ds, int sel, float *buf);
long	*dgetitem_l(DATASET *ds, int sel, long *buf);
short	*dgetitem_s(DATASET *ds, int sel, short *buf);
char	*dgetitem_c(DATASET *ds, int sel, char *buf);

int	ds_ismissing(DATASET *ds);
DATASET	*ds_def(int numitems);
char	*dslayout(DS_DESC *ditem);
char    *ditemalloc();
DATASET	*dsmake(FILE *fp);
DS_DESC	*d_deflt();
DATASET	*dsdup();
DATASET	*scanopen();
SELECTLIST *dsselectlist(DATASET *dsout, DATASET *dsin);
SELECTLIST *dsselparam();
unsigned int ds_minalloc(DATASET *ds);
void	*malloc();
TIME	dstimeok(DATASET *ds, TIME want);
int	dsOKstats(DATASET *ds, char *dsname, int sel);
int	ds_lookup(char *n, nametable *t);
char	*ds_pukool(int a, nametable *t);
void	dsset (DATASET *ds, char *sds);
TIME	atodsindex(int kind, char *stringtime);
int	dput(DATASET *ds, int sel, int ix, int iy, double data);
int	dput_i(DATASET *ds, int sel, int ix, int iy, int data_i);
void	dputitem(DATASET *ds, int sel, double *data);
void	dputitem_d(DATASET *ds, int sel, double *data);
void	dputitem_f(DATASET *ds, int sel, float *data);
void	dputitem_l(DATASET *ds, int sel, long *data);
void	dputitem_s(DATASET *ds, int sel, short *data);
void	dputitem_c(DATASET *ds, int sel, char *data);
void	ds_cleanlbl(DATASET *ds);
void	ds_cpgrp(DATASET *dsin, DATASET *dsout);
void	ds_cpmgrp(DATASET *dsin, DATASET *dsout);
DATASET *dsdup(DATASET *ds);
void	ds_err(DATASET *ds);
void	ds_fillin(DATASET *ds, char *key, char *tok);
int	ds_geth(DATASET *ds);
void	ds_gtxt(DATASET *ds);
int	ds_scantxt(DATASET *ds, FILE *fp);
void	ds_help(DATASET *ds, char *str);
int	ds_init(DATASET *ds);
void	ds_ptxt(DATASET *ds);
void	ds_putkeys(DATASET *ds);
void	ds_fprintkeys(DATASET *ds, FILE *fp);
int	ds_sizetext(DATASET *ds);
int	ds_puth(DATASET *ds);
DS_DESC	*dsadditem(DATASET *ds, DS_DESC *desc);
DATASET *dsgetalld(char *fn, double **dp, int *np);
DATASET *dsgetall(char *fn, float **dp, int *np);
int	dsgrp(DATASET *ds);
void	ds_flush(DATASET *ds);
void	ds_clean(DATASET *ds);
void	dsclose(DATASET *ds);
char	*ds_info(DATASET *ds, char *key);
char	*ds_info_value(DATASET *ds, char *key);
char	*ds_info_comment(DATASET *ds, char *key);
void	ds_info_new(DATASET *ds, char *key, char *value);
void	ds_info_newval(DATASET *ds, char *key, char *val, char *comment);
void	dscopy(DATASET *ds, char *key, DATASET *oldds);
int	wdparam(DATASET *ds, double *offset);
void	dsmk (char *filename,char *name,char *type,TIME start,TIME inc,float *darray,int n);
void	dsmkd (char *filename,char *name,char *type,TIME start,TIME inc,double *darray,int n);
int	dsnew(DATASET *ds, char *fname, char *mode);
DATASET *dsopen(char *str, char *mode);
int	dsparam(DATASET *ds, char *par, int def);
void	dsplotmap(DATASET *ds,TIME tt,int sel, double *levels);
void	printtds(DATASET *ds, TIME B);
void	fprinttds(DATASET *ds, FILE *A, TIME B);
void	dsfprint(DATASET *ds,FILE *fp);
void	dfprint(DATASET *ds, FILE *fp,int sel);
void	dfprintq(DATASET *ds, FILE *fp, int sel);
void	dsplot1D(DATASET *ds, TIME start, TIME stop, int sel, double *levels);
void	dprint(DATASET *ds, int sel, int ix, int iy, FILE *fp);
void	ditemprint(DATASET *ds, int sel, FILE *fp);
void	dsprintgrp(DATASET *ds, FILE *fp);
int	dscan(DATASET *ds, int sel, int ix, int iy, FILE *fp);
int	ditemscan(DATASET *ds, int sel, FILE *fp);
int	dsscangrp(DATASET *ds, FILE *fp);
int	dssametime(DATASET *ds, TIME a, TIME b);
int	dsseek(DATASET *ds, TIME tm);
SELECTLIST *dsselparam(char *sel, DATASET *dsin);
void	dstitle(DATASET *ds, char *progname, char *params);
int	nmtosel(DATASET *ds, char *nm);


#endif
