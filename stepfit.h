#ifndef STEPFIT_H_INCLUDED
#define STEPFIT_H_INCLUDED
#include <stdlib.h>
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

struct data {
    uint64_t n;
    double *y;
};

int expb_f (const gsl_vector * x, void *data, gsl_vector * f);
int expb_df (const gsl_vector * x, void *data, gsl_matrix * J);
int expb_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J);


#endif //STEPFIT_H_INCLUDED
