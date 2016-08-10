#ifndef STEPFIT_H_INCLUDED
#define STEPFIT_H_INCLUDED
#include <stdlib.h>
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include<math.h>
#include"utils.h"
#include"detector.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
//test
struct data {
    int64_t n;
    double *y;
    double *weight;
};

int expb_f (const gsl_vector * x, void *data, gsl_vector * f);
int expb_df (const gsl_vector * x, void *data, gsl_matrix * J);
int expb_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J);
void step_response(event *current, double risetime, int64_t maxiters, double minstep);

#endif //STEPFIT_H_INCLUDED
