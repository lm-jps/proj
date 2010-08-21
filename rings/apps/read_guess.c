#include <stdio.h>
#ifndef max
   #define max(a,b) ((a) > (b) ? (a) : (b))
#endif
#ifndef min
   #define min(a,b) ((a) < (b) ? (a) : (b))
#endif

/* program to evaluate a polynomial with coefficients c at x
    where the polynomial is of the form c(n)*x^n + c(n-1)*x^(n-1) + ... c(0)
    the n passed to this function needs to be = len-1 of the coef array */

/*   double polyval(int n, double **c, double x) {
     double poly;
     poly = (*c)[n];
     while (n) poly = poly*x + (*c)[--n];
     return poly;
   }  */
   double polyval(int n, double *c, double x) {
     double poly;
     poly = c[n];
     while (n) poly = poly*x + c[--n];
     return poly;
   }  

   int max_arr(int *arr, int n) {
     int i, tmp;
     tmp = arr[0];
     for (i=1; i < n; i++) { 
       tmp = max(tmp, arr[i]);
     }
     return tmp;
   }
   int min_arr(int *arr, int n) {
     int i, tmp;
     tmp = arr[0];
     for (i=1; i < n; i++) { 
       tmp = min(tmp, arr[i]);
     }
     return tmp;
   }
    

 int read_guess(char *guessfile, int nrdtot, int *m, double **fnua, 
                double **fwa, double **pa, double **bga, int *ma, int verbose) {

      FILE *gfile;
      int ma1, i, n, status, ip; 

      gfile = fopen (guessfile, "r"); 
      if (!gfile) { 
        fprintf (stderr, "Error: unable to open guess table %s\n", guessfile);
        return 1;
      }
      fscanf(gfile, "%d %d ", &nrdtot, &ma1);
      for (i=0; i < 4; i++) fscanf(gfile,"%d", &(m[i])); 
      if (verbose) {
        printf("guessfile = %s \n", guessfile);
        printf("nrdtot = %d ma = %d \n",nrdtot,ma1);
        printf("m \n");
        for (i=0; i < 4; i++) printf("%d", (m[i]));
        printf("\n\n");
      }
      *fnua = (double*) malloc(ma1 * nrdtot * sizeof(double));
      *fwa = (double*) malloc(ma1 * nrdtot * sizeof(double));
      *pa = (double*) malloc(ma1 * nrdtot * sizeof(double));
      *bga = (double*) malloc(ma1 * nrdtot * sizeof(double));

      for (i=0; i < nrdtot; i++) {
         for (n=0; n < ma1; n++) fscanf(gfile,"%lf", &(*fnua)[i*ma1 + n]);  
            if (fscanf == NULL) { 
              fprintf (stderr, "Problem reading in fnua\n");
              return 1;
         }
      }
      for (i=0; i < nrdtot; i++) {
         for (n=0; n < ma1; n++) fscanf(gfile,"%lf", &(*fwa)[i*ma1 + n]);  
         if (fscanf == NULL) { 
           fprintf (stderr, "Problem reading in fwa\n");
           return 1;
         }
      }
      for (i=0; i < nrdtot; i++) {
         for (n=0; n < ma1; n++) fscanf(gfile,"%lf", &(*pa)[i*ma1 + n]);  
         if (fscanf == NULL) { 
           fprintf (stderr, "Problem reading in pa\n");
           return 1;
         }
      }
      for (i=0; i < nrdtot; i++) {
         for (n=0; n < ma1; n++) fscanf(gfile,"%lf", &(*bga)[i*ma1 + n]);  
         if (fscanf == NULL) { 
           fprintf (stderr, "Problem reading in bga\n");
           return 1;
         }
      }
      if(verbose) { 
        printf("fnua\n");
        for (i=0; i < nrdtot; i++) {
         for (n=0; n < ma1; n++) printf("%12.5f", (*fnua)[i*ma1 + n]);  
         printf("\n");
        }
        printf("\n, \n");
        printf("fwa\n");
        for (i=0; i < nrdtot; i++) {
         for (n=0; n < ma1; n++) printf("%12.5f", (*fwa)[i*ma1 + n]);  
         printf("\n");
        }
        printf("\n, \n");
        printf("pa\n");
        for (i=0; i < nrdtot; i++) {
         for (n=0; n < ma1; n++) printf("%12.5f", (*pa)[i*ma1 + n]);  
         printf("\n");
        }
        printf("\n, \n");
        printf("bga\n");
        for (i=0; i < nrdtot; i++) {
         for (n=0; n < ma1; n++) printf("%12.5f", (*bga)[i*ma1 + n]);  
         printf("\n");
        }
      }
      *ma = ma1;
      return 0;
 }
 int make_table(double *fnua, double *fwa, double *pa, double *bga,
                double *fnu, double *width, double *amp, double *bkg,
                int *m, int nk, double dk, int nrdtot, int kbmin,
                int kbmax, int ma, int verbose)  {
      int i, j, m1, m2, m3, m4, n, ma1;
      double logk;
      ma1 = ma;
      m1 = m[0] - 1;
      m2 = m[1] - 1;
      m3 = m[2] - 1;
      m4 = m[3] - 1;

      for (i=0; i < nrdtot; i++) {
         for (j = kbmin-1; j <= kbmax-1; j++) {
           logk= log(dk * (j+1));
           fnu[i * nk + j] = polyval(m1, &(fnua[i*ma1]), logk);
         }
 
         for (j = kbmin-1; j <= kbmax-1; j++) {
           logk= log(dk * (j+1));
           width[i * nk + j] = polyval(m2, &(fwa[i*ma1]), logk);
         }

         for (j = kbmin-1; j <= kbmax-1; j++) {
           logk= log(dk * (j+1));
           amp[i * nk + j] = polyval(m3, &(pa[i*ma1]), logk);
         }

         for (j = kbmin-1; j <= kbmax-1; j++) {
           logk= log(dk * (j+1));
           bkg[i * nk + j] = polyval(m4, &(bga[i*ma1]), logk);
         }
      }
     /* exponentiate fnu but not the other parameters  */ 

      for (i=0; i < nrdtot; i++) {
        for (j = kbmin-1; j <= kbmax-1; j++)  {
          fnu[i*nk + j] = exp(fnu[i*nk + j]);
        }
      }
/*      if (verbose) {
       for (i=0; i < nrdtot; i++) {
        printf("fnu   n = %d\n", i);
        for (j = kbmin-1; j <= kbmax-1; j++)  {
          printf(" %lf  ", fnu[i*nk + j]);
        }
        printf("\n\n");    
       }
      
       for (i=0; i < nrdtot; i++) {
        printf("width   n = %d\n", i);
        for (j = kbmin-1; j <= kbmax-1; j++)  {
          printf(" %lf  ", width[i*nk + j]);
        }
        printf("\n\n");    
       }
       for (i=0; i < nrdtot; i++) {
        printf("amp   n = %d\n", i);
        for (j = kbmin-1; j <= kbmax-1; j++)  {
          printf(" %lf  ", amp[i*nk + j]);
        }
        printf("\n\n");    
       }
       for (i=0; i < nrdtot; i++) {
        printf("bkg   n = %d\n", i);
        for (j = kbmin-1; j <= kbmax-1; j++)  {
          printf(" %lf  ", bkg[i*nk + j]);
        }
        printf("\n\n");    
       }
        printf("still in make_table kbmin = %d   kbmax = %d \n", kbmin, kbmax);
      }  */
      return 0;
 }     
