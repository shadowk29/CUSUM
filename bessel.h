/*
bessel.c is a C port of select functions from the
python library bessel filter implementation, scipy.signal.bessel.

Copyright (c) 2001, 2002 Enthought, Inc.
All rights reserved.

Copyright (c) 2003-2016 SciPy Developers.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

  a. Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.
  b. Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
  c. Neither the name of Enthought nor the names of the SciPy Developers
     may be used to endorse or promote products derived from this software
     without specific prior written permission.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
THE POSSIBILITY OF SUCH DAMAGE.
*/

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
    double *temp;
};
typedef struct Bessel bessel;

void besselap(int64_t N, double complex *poles, double complex *zeros);
void polymult(double complex *roots, int64_t N, double *rcoefs);
double scale_filter(double complex *poles, int64_t N, double warped, double scale);
double bilinear(double complex *poles, int64_t N, double scale, double fs);
void transform_filter(double complex *poles, double complex *zeros, int64_t N, double scale, double *b, double *a);
bessel *initialize_filter(int usefilter, int eventfilter, int64_t order, double cutoff, int64_t length);
void free_filter(bessel *lpfilter);
void filter_signal(double *signal, double *paddedsignal, bessel *lpfilter, int64_t length);

#endif //BESSEL_H_INCLUDED
