#include "stepfit.h"

int expb_f (const gsl_vector * x, void *data, gsl_vector * f)
{
    uint64_t n = ((struct data *)data)->n;
    double *y = ((struct data *)data)->y;

    double i0 = gsl_vector_get (x, 0);
    double a = gsl_vector_get (x, 1);
    double u1 = gsl_vector_get (x, 2);
    double tau1 = gsl_vector_get (x, 3);
    double b = gsl_vector_get (x, 4);
    double u2 = gsl_vector_get (x, 5);
    double tau2 = gsl_vector_get (x, 6);
    uint64_t i;
    for (i = 0; i < n; i++)
    {
        double t = i;
        double Yi = i0;
        if (t > u1)
        {
            Yi += a*(exp(-(t-u1)/tau1)-1);
        }
        if (t > u2)
        {
            Yi += b*(1-exp(-(t-u2)/tau2));
        }
        gsl_vector_set (f, i, (Yi - y[i]));
    }
    return GSL_SUCCESS;
}

int expb_df (const gsl_vector * x, void *data, gsl_matrix * J)
{
    uint64_t n = ((struct data *)data)->n;

    double a = gsl_vector_get (x, 1);
    double u1 = gsl_vector_get (x, 2);
    double tau1 = gsl_vector_get (x, 3);
    double b = gsl_vector_get (x, 4);
    double u2 = gsl_vector_get (x, 5);
    double tau2 = gsl_vector_get (x, 6);

    uint64_t i;

    for (i = 0; i < n; i++)
    {
        double t = i;
        double s = 1;
        gsl_matrix_set (J, i, 0, 1);
        gsl_matrix_set (J, i, 1, t < u1 ? 0 : 1/s*(exp(-(t-u1)/tau1)-1));
        gsl_matrix_set (J, i, 2, t < u1 ? 0 : 1/s*(a/tau1 * exp(-(t-u1)/tau1)));
        gsl_matrix_set (J, i, 3, t < u1 ? 0 : 1/s*(a * (t-u1)/(tau1*tau1) * exp(-(t-u1)/tau1)));
        gsl_matrix_set (J, i, 4, t < u2 ? 0 : 1/s*(1-exp(-(t-u2)/tau2)));
        gsl_matrix_set (J, i, 5, t < u2 ? 0 : 1/s*(b/tau2 * (-exp(-(t-u2)/tau2))));
        gsl_matrix_set (J, i, 6, t < u2 ? 0 :  1/s*(-b * (t-u2)/(tau2*tau2) * exp(-(t-u2)/tau2)));
    }
    return GSL_SUCCESS;
}

int expb_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J)
{
    expb_f (x, data, f);
    expb_df (x, data, J);

    return GSL_SUCCESS;
}
