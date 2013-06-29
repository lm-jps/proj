#define SDC_DESPIKE_SIG \
  int sdc_despike(\
    int *image,         /* To ble cleaned */				\
    char *mask,         /* Ignore these - guaranteed unchanged */	\
    int nx, int ny,     /* Size (nx,ny) of image */			\
    int neighbour,      /* Number of neighbour iterations */		\
    char *kernel,       /* Neighbour convolution. Generated if needed */ \
    int kernel_size,    /* Note: kernel is square */			\
    int xbox, int ybox, /* Size of median box for detection */		\
    int max_var_low,    /* Max variation above median in "low" areas */	\
    float max_factor_hi,/* Max factor above median in "high" areas */	\
    int limit,          /* Value separating high from low areas */	\
    int *badblobs,      /* Blobs of "bad pixels", convert to missing */ \
    int sizeofbads,     /* Number of pre-flagged bad pixels */		\
    int *nspikes,       /* Number of pixels flagged */			\
    int **oldvalues,    /* Old values of flagged pixels */		\
    int **spikelocs,    /* Their one-dimensional index location */	\
    int **newvalues)    /* And their *new* values */

/* Additional notes:
   
   - Bad blobs are filled in, using the same box median as for spiked pixels,
   but do not appear in the output oldvalues/spikelocs/newvalues.
   
   - The value of BLANK/MISSING is hardcoded at 0x80000000
   
   - Existing BLANK/MISSING pixels *and* masked pixels are *unchanged* upon
   exit - this also goes for pixels flagged as bad blobs - i.e.
   masking/pre-existing missing pixels take precedence.

   - As with the AIA routine, a fake mask is generated if needed

   - The kernel is a *square* array used for flagging neighbours of "primary"
   spike pixels, through convolution. If kernel_size or kernel is zero, a
   default 3x3 crosshair kernel is used.

   - As with the AIA routine, *oldvalues, *spikelocs, and *newvalues should be
   free()'d by the caller. All other malloc'ed pointers are freed internally.
   
*/

   
