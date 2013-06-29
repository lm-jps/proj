#include <stdio.h>
#include <stdlib.h>

#include "sdc_despike_src.h"
#include "fmedian_src.h"
#include "idl_export.h"

#define DEBUG 0
#define BLANK 0x80000000   /* From aia_despike behaviour */
#define MISSING BLANK      /* What I'm used to           */

#define DPRINTF(a) {if (DEBUG) printf a;}
#define PRINTI(var) DPRINTF(("%20s= %ld\n",#var,(long)var));
#define PRINTF(var) DPRINTF(("%20s= %f\n",#var,var));
#define PRINTX(var) DPRINTF(("%20s= x%llx\n",#var,(unsigned long long)var))

/* "Safe" allocation macro */
#define MALLOC(dest,num)						\
  {dest = malloc((num) * sizeof(*dest));				\
    if ( ! dest ) {							\
      printf("malloc error for " #dest " in sdc_despike\n"); return 1; } \
  }

/* When applying kernel, consider only pre-recorded MISSING pixels */
int pre_spike(int loc,int *locs, int last)
{
  int spikeno;
  for (spikeno=0; spikeno<last; spikeno++) if (locs[spikeno] == loc) return 1;
  return 0;
}

SDC_DESPIKE_SIG
{
  int xi,yi;   /* Image indices */
  int loc;     /* Local shorthand for image location counter / xi + nx*yi */
  int badblobix;
  int blobbed = 0;
  
  int free_mask =   0; /* Is mask faked? */
  int free_kernel = 0; /* Ditto for kernel */

  int *locs; /* Shorthands for **spikelocs, **oldvalues, **newvalues */
  int *olds;
  int *news;

  int  ksz = kernel_size; /* Short */

  int *median;
  int *median_workspace;

  int primary_spike;     /* spike count before neighbours are flagged */
  int spike;             /* Spikes registered so far... */

  int fill;              /* Whether to fill one more time or not */

  int  *restore_val;  /* Copy of image for restoration of MISSING/masked */
  char *restore_mask; /* Ditto mask */

  DPRINTF(("sdc_despike_src() called\n"));
  PRINTX(image);
  PRINTX(mask);
  PRINTI(nx);
  PRINTI(ny);
  PRINTI(neighbour);
  PRINTX(kernel);
  PRINTI(kernel_size);
  PRINTI(xbox);
  PRINTI(ybox);
  PRINTI(max_var_low);
  PRINTF(max_factor_hi);
  PRINTI(limit);          /* Value separating high from low areas */
  PRINTX(badblobs);        /* Pre-flagged "bad pixels" - convert to MISSING */
  PRINTI(sizeofbads);     /* Number of pre-flagged bad pixels */
  PRINTX(nspikes);
  PRINTX(*oldvalues);
  PRINTX(*spikelocs);
  PRINTX(*newvalues);

  /* Fake mask if not supplied */
  if (mask==0) {
    free_mask = 1;
    MALLOC(mask,nx * ny); /* to-be-freed - conditionally */
    for (loc=0; loc< nx*ny; loc++) mask[loc] = 1;
    DPRINTF(("Generated fake mask\n"));
  }
  
  /* Copy image for restoring pre-existing MISSING and/or */
  /* masked pixels, set restore mask for those pixels */
  MALLOC(restore_val, nx * ny);   /* to-be-freed */
  MALLOC(restore_mask, nx * ny);  /* to-be-freed */
  for (loc=0; loc < nx * ny; loc++) {
    restore_val[loc] = image[loc];
    restore_mask[loc] = ( image[loc]==MISSING  ||  !mask[loc] );
  }
    
  /* Deal w/kernel if not supplied; ksz and kern already init'ed */
  if ( !ksz || !kernel) {
    int i;
    
    ksz = 3;
    MALLOC(kernel, ksz * ksz); /* to-be-freed - conditionally */
    free_kernel=1;
    for (i=0; i < ksz*ksz; kernel[i++] = 0);

    kernel[1 + 0*3]=1; /* Cross, vertical line */
    kernel[1 + 2*3]=1;
    kernel[0 + 1*3]=1; /* Cross, horizontal line */
    kernel[2 + 1*3]=1;
  }

  /* Definitely enough space for spikelocs/olvalues/newvalues, shorthands */
  MALLOC(*spikelocs, nx * ny); locs = *spikelocs;  /* pass-back */
  MALLOC(*oldvalues, nx * ny); olds = *oldvalues;  /* pass-back */
  MALLOC(*newvalues, nx * ny); news = *newvalues;  /* pass-back */

  /* Allocate output & workspace for median routine                */
  MALLOC(median, nx * ny);                            /* to-be-freed */
  MALLOC(median_workspace, xbox * ybox);              /* to-be-freed */


  /* Doit ***********************************************/

  DPRINTF(("Blobbing - flag as MISSING but don't record\n"));
  /* [nperim, perim_offsets, npix, pix_offsets] */
  badblobix=0;
  while (badblobix<sizeofbads) {
    int nperim = badblobs[ badblobix++ ];
    int npix;
    badblobix += nperim;
    npix = badblobs[ badblobix++ ];

    while (npix--) {
      int loc = badblobs[ badblobix++ ];
      image[loc]=MISSING; /* DO NOT USE HANDLE_SPIKE! */
      blobbed++;
    }
  }
  DPRINTF(("Blobbed: %d\n",blobbed));
  loc = nx*ny - 1;

  FMEDIAN_CALL(IDL_LONG,image,median,median_workspace,nx,ny,xbox,ybox,
	       1/*always*/, MISSING, 0/*don't-find_MISSING*/);

  /**** Calculate median array */
  
  /* MAIN ALGORITHM - detect pixels, flag as MISSING, count in spike */

#define HANDLE_SPIKE(spike,loc,locs,olds,image)			\
  {									\
    locs[spike] = loc;						\
    olds[spike] = image[loc];    /* Store org. value */		\
    image[loc] = MISSING;           /* Flag */				\
    spike++;                        /* Detected one - count it! */	\
  }								       

  spike=0;

  for (yi=0; yi<ny; yi++) { /* Should we do a single loop w/loc<nlocs? */
    for (xi=0; xi<nx; xi++) {
      loc = xi + yi*nx;
      if (mask[loc] && image[loc] != MISSING) {
	if (image[loc] < limit) {                          /* Low data */
	  if ( (image[loc]-median[loc]) > max_var_low ) {
	    HANDLE_SPIKE(spike,loc,locs,olds,image);
	    DPRINTF(("HANDLED (%d,%d)\n",xi,yi));
	  }
	} else {                                           /* High data */
	  float med = median[loc] ? median[loc] : 1.0;
	  if ( image[loc]/med > max_factor_hi) {
	    HANDLE_SPIKE(spike,loc,locs,olds,image);
	    DPRINTF(("HANDLED (%d,%d)\n",xi,yi));
	  }
	}
      }
    }
  }
  DPRINTF(("Main algorithm done\n"));

  DPRINTF(("Primary spikes (no neighbour): %d\n",spike));
  primary_spike = spike;

  /* When applying kernel, we should *NOT* flag pixels as MISSING */
  /* during each iteration - since that would "spread" the newly  */
  /* flagged pixels! */

  while (neighbour--) {           /* Allow for multiple neighbour iterations */
    int midk = ksz/2;

    /* Pixels flagged *during* this run should not be regarded as MISSING */
    /* during this run - so verify that any MISSING neighbour was already */
    /* recorded in locs[ first_spike...last_spike-1 ] */

    int last_spike = spike; /* Keep for run-through flagging */


    /* Convolve w/kernel & flag */
    for (yi=0; yi<ny; yi++) {     
      for (xi=0; xi<nx; xi++) {
	/* Calc. bounds of kernel indices for this image pixel */

	/* Ideal stop at < xi-midk+ksz, but actual stop is bounded */
	int max_xi = MIN2(nx, xi-midk+ksz);

	/* Ideal max_kern_xi: ksz. But subtract diff between ideal image stop */
	/* and actual image stop */
	int max_kern_xi = ksz - ((xi-midk+ksz) - max_xi);

	/* Ideal start at xi-midk, but actual stop is bounded */
	int min_xi = MAX2( 0, xi-midk);

	/* Ideal min_kern_xi: 0 but *add* diff between ideal start & actual */
	/* start (reverse of above) */
	int min_kern_xi = min_xi - (xi-midk);

	/* Identical arguments for y direction: */
	
	int max_yi = MIN2(ny, yi-midk+ksz);
	int max_kern_yi = ksz - ((yi-midk+ksz) - max_yi);
	int min_yi = MAX2( 0, yi-midk);
	int min_kern_yi = min_yi - (yi-midk);

	int kern_xi,kern_yi; 	  /* Kernel indices */
		
	loc = xi + nx*yi;
	
	/* Don't add/touch already MISSING or masked pxls */
	if (image[loc] == MISSING || !mask[loc]) {
	  char *why = mask[loc] ? "MISSING" : "MASKED";
	  DPRINTF(("(%d,%d; %d) already %s - skipping\n",xi,yi,loc,why));
	  continue;
	}

	/* Loop over all allowed kernel pixels */
	for (kern_yi=min_kern_yi; kern_yi<max_kern_yi; kern_yi++) {
	  for (kern_xi=min_kern_xi; kern_xi<max_kern_xi; kern_xi++) {
	    
	    if (kernel[kern_xi + ksz*kern_yi]) {
	      
	      /* kernel[0,0] in 3x3 kernel refers */
	      /* to "is neighbour one-down, one-left flagged as spike" */
	      
	      int neighb_xi = (xi + kern_xi - midk);
	      int neighb_yi = (yi + kern_yi - midk);
	      int neighb_loc = neighb_xi + nx * neighb_yi;
	      
	      if (image[neighb_loc] == MISSING) {
		char *what = "but not prev-round spike";
		DPRINTF(("(%d,%d; %d) neighb. (%d,%d; %d) is MISSING ",
			 xi,yi,loc,neighb_xi,neighb_yi,neighb_loc)); 
		if (pre_spike(neighb_loc,locs,last_spike)) {
		  what = "and a last-round spike";
		  HANDLE_SPIKE(spike,loc,locs,olds,image);
		  goto next_image_pixel; /* Don't do this one again */
		}
		DPRINTF((" - %s\n",what));
	      }
	    }
	  }
	}
      next_image_pixel: ;
      }
    }
    DPRINTF(("Neighbour iteration done, flagged %d\n",spike-last_spike));
  }
  DPRINTF(("Neighbouring done, total flagged %d\n",spike-primary_spike));
  
  *nspikes = spike; /* Record the final number */
  
  /* Fill all spikes, if requested */

  fill=1;
  while (fill) {
    fill=0; /* Assume we'll fix all MISSING pixels */
    
    /* Calc. new median (w/only-MISSING set) */
    FMEDIAN_CALL(IDL_LONG,image,median,median_workspace,nx,ny,xbox,ybox,
		 1/*always*/,MISSING,0/*don't-find_MISSING*/);
    
    /* Fill in the ones where we now have a median */
    for (loc=0; loc<nx*ny; loc++) {
      /* Needs to be filled? */
      if (image[loc] == MISSING) {                 
	if (median[loc] != MISSING) image[loc] =  median[loc];    /* Fill */
	else fill=1;                             /* Do it next iteration! */
      }
    }
  }
  DPRINTF(("\nFilled\n"));

  /* Put back masked & originally-flagged-as-MISSING pixels */
  for (loc=0; loc < nx*ny; loc++) {
    if (restore_mask[loc]) image[loc] = restore_val[loc];
  }
  
  /* Freeing MALLOC'ed stuff that's not passed back to caller */
  if (free_mask) free(mask);
  if (free_kernel) free(kernel);
  free(restore_val);
  free(restore_mask);
  free(median);
  free(median_workspace);
  
  /* Record new values in newvalues - bizarre */
  for (spike=0; spike<*nspikes; spike++) {
    loc = locs[spike];
    news[spike] = image[loc]; /* Bizarre ;-) */
  }
  /* Done it(?) *******************************************/
  
  return 0;
}
