#ifndef BESSEL_H_INCLUDED
#define BESSEL_H_INCLUDED
#include"utils.h"
#include<stdio.h>
#include<stdlib.h>
#include<complex.h>
#include<math.h>
#include<inttypes.h>



void besselap(uint64_t N, double complex *poles, double complex *zeros);
void polymult(double complex *roots, uint64_t N, double *rcoefs);
double scale_filter(double complex *poles, uint64_t N, double warped, double scale);
double bilinear(double complex *poles, uint64_t N, double scale, double fs);
void transform_filter(double complex *poles, double complex *zeros, uint64_t N, double scale, double *b, double *a);
bessel *initialize_filter(bessel *lpfilter, uint64_t order, double cutoff, uint64_t length);
void free_filter(bessel *lpfilter);

#endif //BESSEL_H_INCLUDED
