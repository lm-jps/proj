/*these are structures used in many of the ANA modules */
#define ABS(x) ((x)>=0?(x):-(x))
#define	MIN(a,b) (((a)<(b))?(a):(b))
#define	MAX(a,b) (((a)>(b))?(a):(b))
/*plus some parameters we need in some */
/* note that NSCRAT is a count of int's in scrat */
#define NSCRAT 8196
#define	MAX_USER_SUBRS	900
#define	MAX_USER_FUNCS	600
#define	MAX_USER_BLOCKS	600
#define MXNEST 5
#define ANA_SUBSC_FUN 1
#define ANA_POW_FUN 3
typedef unsigned char byte;
struct	sdesc	{ int n;	byte *p; };
/*symbol descriptor structure, 4 types are included - scalar, array (strings
	have same symbol structure with a different class), evb, and edb
	the xx field is used for a hash # for named symbols */
struct	sym_desc { byte class,type; short xx;
	union { union { byte b; short w; int l; float f; double d;} scalar;
		struct { int *ptr; int bstore; } array;
		struct { short args[4]; } evb;
		struct { short args[4]; } edb; } spec;};
/*note, each desc. header is 12 bytes currently for 32 bit machines
and 16 byte for 64 bit machines (because of the 8 byte ptr's) */
struct ana_subr_struct {
	char	*name;		/* ptr to name of subroutine */
	int	minargs;	/* minimum # of args */
	int	maxargs;	/* maximum ... */
	int	(*ptr)();	/* ptr to the subroutine code */
	};
/*ahead is a structure used for array headers */
struct	ahead {	byte ndim, c1, c2, c3;	/*c's are for characteristics
						not yet defined */
		int dims[8];		/*space for 8 dimensions */
		struct flist *facts;};	/*a pointer to things known about
					the array */
struct	flist {	int type;
/*for linked lists of "facts" such as max, min, mean, whatever */
		union { byte b; short w; int l; float f; double d;} value;
		struct flist *next; };
/*scalar union definition */
union scalar { byte b; short w; int l; float f; double d; };
/*linked structure for user subroutine names and first line */
struct user_subr_table {
	char *name; char *line; int num; struct user_subr_table *next; };
/*symbol list structure for subrs, funcs, etc. */
struct sym_list { int num; struct sym_list *next; };
