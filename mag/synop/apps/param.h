/*
$Header: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/mag/synop/apps/param.h,v 1.1 2014/02/11 18:30:19 arta Exp $
$Log: param.h,v $
Revision 1.1  2014/02/11 18:30:19  arta
Finish removing dependencies on Yangs home directory. Added several headers needed for the synop project, but that resided in Yangs myproj directory.

 * Revision 1.6  1992/03/06  00:55:11  todd
 * Fixed def of strdup to work on Banneker
 *
 * Revision 1.5  1991/12/17  20:48:18  phil
 * added strdup()
 *
 * Revision 1.4  1991/12/11  17:27:54  phil
 * *** empty log message ***
 *
 * Revision 1.3  90/10/22  14:45:21  phil
 * extract timeparam()
 * delete s*param()
 * 
 * Revision 1.2  90/10/12  11:58:15  phil
 * add align_2,4,8
 * 
 * Revision 1.1  90/09/17  16:40:29  phil
 * Initial revision
 * 
*/

/* external definitions for routines in ./src/libdate.d */
extern double parameter();
extern char *nameparam();
extern int is_param();
extern int openparam();
extern FILE *fopenparam();
extern char *which();
extern char *open_nameparam();
extern double nthparam();
extern double *listparameter();
extern double *listparam();
extern char _allargs[];
extern char *string(char *format, ... );
extern int Strncmp();
extern int Strcmp();
extern int nl();
extern int cscanf();
extern void sprintfreq();
extern void fprintfreq();
extern void tofrom();
extern void init_zero();
#ifndef __linux
extern char *strdup();
#endif
extern char *nthname();

extern char *mytty();

extern short	align_s();
extern long	align_l();
extern float	align_f();
extern double	align_d();
extern char	*align_2();
extern char	*align_4();
extern char	*align_8();

int is_flag();
void getargs();
void setargs();
void addparam();


