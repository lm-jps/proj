static double median(int n, float arr[]);

int fstats(int n, float arr[], double *min, double *max, double *medn,
	   double *mean, double *sig, double *skew, double *kurt, int *ngood);

int fstats2(int n, float arr[], double *min, double *max, double *medn,
	   double *mean, double *sig, double *skew, double *kurt, int *ngood);

int dstats(int n, double arr[], double *min, double *max, double *medn,
	   double *mean, double *sig, double *skew, double *kurt, int *ngood);

int dstats2(int n, double arr[], double *min, double *max, double *medn,
	   double *mean, double *sig, double *skew, double *kurt, int *ngood);

int set_statistics(DRMS_Segment_t *seg, DRMS_Array_t *data, int mode);

