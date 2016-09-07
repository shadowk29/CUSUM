#ifndef BESSEL_H_INCLUDED
#define BESSEL_H_INCLUDED
#include"utils.h"
#include<stdio.h>
#include<stdlib.h>
#include<complex.h>
#include<math.h>
#include<inttypes.h>

struct Bessel
{
    double *dcof;
    double *ccof;
    double cutoff;
    int64_t order;
    int64_t padding;
    double *paddedsignal;
    double *temp;
};
typedef struct Bessel bessel;

void besselap(int64_t N, double complex *poles, double complex *zeros);
void polymult(double complex *roots, int64_t N, double *rcoefs);
double scale_filter(double complex *poles, int64_t N, double warped, double scale);
double bilinear(double complex *poles, int64_t N, double scale, double fs);
void transform_filter(double complex *poles, double complex *zeros, int64_t N, double scale, double *b, double *a);
bessel *initialize_filter(bessel *lpfilter, int64_t order, double cutoff, int64_t length, int64_t samplingfreq);
void free_filter(bessel *lpfilter);
void filter_signal(double *signal, bessel *lpfilter, int64_t length);

#endif //BESSEL_H_INCLUDED
