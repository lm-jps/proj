#ifndef KAHAN_SUM
#define KAHAN_SUM
#endif

void kahan_sum (double *arr, int len, double *sum)
{
    double c = 0.0;
    double sum_t = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:sum_t,c)
#endif
    for (int i = 0; i < len; i++) {
        double y = arr[i] - c;
        double t = sum_t + y;
        c = (t - sum_t) - y;
        sum_t = t;
    }
    *sum = sum_t;
}

void kahan_sum_sgl (double *arr, int len, double *sum)
{
    double c = 0.0;
    double sum_t = 0.0;
    for (int i = 0; i < len; i++) {
        double y = arr[i] - c;
        double t = sum_t + y;
        c = (t - sum_t) - y;
        sum_t = t;
    }
    *sum = sum_t;
}

