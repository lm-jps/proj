/*
$Header: /home/akoufos/Development/Testing/jsoc-4-repos-0914/JSOC-mirror/JSOC/proj/mag/synop/apps/plt.h,v 1.1 2014/02/11 18:30:19 arta Exp $
$Log: plt.h,v $
Revision 1.1  2014/02/11 18:30:19  arta
Finish removing dependencies on Yangs home directory. Added several headers needed for the synop project, but that resided in Yangs myproj directory.

 * Revision 1.1  90/09/17  16:40:38  phil
 * Initial revision
 * 
*/

#include "/home/wso/include/plot.h"
#include <stdtime.h>

PFILE *openplot();
PFILE *initplot();
PFILE *openinitplot();
int addplot();
/*
#define pointplot(a,b,c,d) lineplot(a,b,c,b,c,d)
*/
#define beginp(a) beginplot(stdplt,a)
#define defp(a,b,c,d) defplot(stdplt,a,b,c,d)
#define drawp(a,b,c) drawplot(stdplt,a,b,c)
#define endp() endplot(stdplt)
#define gridp(a,b,c,d) gridplot(stdplt,a,b,c,d)
#define pointp(a,b,c) lineplot(stdplt,a,b,a,b,c)
#define linep(a,b,c,d,e) lineplot(stdplt,a,b,c,d,e)
#define modep(a) modeplot(stdplt,a)
#define movep(a,b,c) moveplot(stdplt,a,b,c)
#define setp(a,b,c,d) setplot(stdplt,a,b,c,d)
#define getp(a,b,c) getplot(stdplt,a,b,c)
#define addp(a) addplot(stdplt,a)
#define setstd(P) setplot(P,P->p_xsize/8.0,P->p_ysize/8.0,\
P->p_xsize*0.75,P->p_ysize*0.75)
#define marginplot(P,L,B,R,T) setplot(P,L,B,P->p_xsize-R,P->p_ysize-T)
