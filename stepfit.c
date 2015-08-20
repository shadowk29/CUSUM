/*

                                COPYRIGHT
    Copyright (C) 2015 Kyle Briggs (kbrig035<at>uottawa.ca)

    This file is part of CUSUM.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "stepfit.h"

int expb_f (const gsl_vector * x, void *data, gsl_vector * f)
{
    uint64_t n = ((struct data *)data)->n;
    double *y = ((struct data *)data)->y;
    double *weight = ((struct data *)data)->weight;

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
        gsl_vector_set (f, i, (Yi - y[i])/weight[i]);
    }
    return GSL_SUCCESS;
}

int expb_df (const gsl_vector * x, void *data, gsl_matrix * J)
{
    uint64_t n = ((struct data *)data)->n;
    double *weight = ((struct data *)data)->weight;
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
        double s = weight[i];
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


int step_response(event *current, double risetime, uint64_t maxiters, double minstep)
{
    const gsl_multifit_fdfsolver_type *T;
    gsl_multifit_fdfsolver *s;
    int status = GSL_CONTINUE;
    uint64_t i,iter = 0;
    uint64_t p = 7;
    uint64_t n = current->length + current->padding_before + current->padding_after;
    double *weight;
    if ((weight = malloc(n*sizeof(double)))==NULL)
    {
        printf("Cannot allocate error array\n");
        abort();
    }

    gsl_matrix *covar = gsl_matrix_alloc (p, p);

    gsl_multifit_function_fdf f;



    double maxsignal = signal_max(current->signal, current->length + current->padding_before + current->padding_after);
    double minsignal = signal_min(current->signal, current->length + current->padding_before + current->padding_after);
    double baseline = signal_average(current->signal,current->padding_before);
    int sign = signum(baseline);
    uint64_t start = current->padding_before - intmin(current->length/2, current->padding_before/4);
    uint64_t end = current->padding_before+current->length;



    for (i=0; i<n; i++)
    {
        if (i < start)
        {
            weight[i] = 1;
        }
        else if (i < end)
        {
            weight[i] = 0.3;
        }
        else
        {
            weight[i] = 1;
        }
    }
    struct data d = {n, current->signal,weight};

    gsl_vector *x = gsl_vector_alloc(p);
    gsl_vector_set(x,0,baseline);
    gsl_vector_set(x,1,sign < 0 ? maxsignal - minsignal: minsignal - maxsignal);
    gsl_vector_set(x,2,start);
    gsl_vector_set(x,3,risetime);
    gsl_vector_set(x,4,sign < 0 ? maxsignal - minsignal: minsignal - maxsignal);
    gsl_vector_set(x,5,end - current->length/2);
    gsl_vector_set(x,6,risetime);


    f.f = &expb_f;
    f.df = &expb_df;
    f.fdf = &expb_fdf;
    f.n = n;
    f.p = p;
    f.params = &d;

    T = gsl_multifit_fdfsolver_lmsder;
    s = gsl_multifit_fdfsolver_alloc (T, n, p);
    gsl_multifit_fdfsolver_set (s, &f, x);

    while (iter < maxiters && status == GSL_CONTINUE)
    {
        iter++;
        status = gsl_multifit_fdfsolver_iterate (s);
        if (status)
        {
            break;
        }
        status = gsl_multifit_test_delta (s->dx, s->x,1e-3, 1e-3);
    }


    double i0 = FIT(0);
    double a = FIT(1);
    uint64_t u1 = (uint64_t) FIT(2);
    double rc1 = FIT(3);
    double b = FIT(4);
    uint64_t u2 = (uint64_t) FIT(5);
    double rc2 = FIT(6);

    //printf("i0 = %g\na=%g\nu1=%"PRIu64"\nb=%g\nu2=%"PRIu64"\n",i0,a,u1,b,u2);

    if (u1 > u2) //if the fit got stuck in a variable swapped equivalent minimum, we can just reverse the paramters
    {
        uint64_t temp;
        double dtemp;
        temp = u1;
        u1 = u2;
        u2 = temp;
        dtemp = a;
        a = -b;
        b = -dtemp;
        dtemp = rc1;
        rc1 = rc2;
        rc2 = dtemp;
    }
    if (signum(a) != signum(b)) //if it is not a step-and-return event
    {
        current->type = BADFIT;
    }
    if (d_abs(b) < minstep || d_abs(b) < minstep) //if it is not a step-and-return event
    {
        current->type = BADFIT;
    }
    if (u2 > n || u1 <= 0) //if for some reason we are out of range
    {
        current->type = BADFIT;
        return BADFIT;
    }
    for (i=0; i<u1; i++) //if all went well, populate the filtered trace
    {
        current->filtered_signal[i] = i0;
    }
    for (i=u1; i<u2; i++)
    {
        current->filtered_signal[i] = i0-a;
    }
    for (i=u2; i<n; i++)
    {
        current->filtered_signal[i] = i0-a+b;
    }
    current->rc1 = rc1;
    current->rc2 = rc2;



    gsl_multifit_fdfsolver_free (s);
    gsl_matrix_free (covar);
    gsl_vector_free(x);
    free(weight);
    return status;
}

void step_response_events(event *current, double risetime, uint64_t maxiters, double minstep)
{
    event *head = current;
    uint64_t numevents = 0;
    int status;
    while (current)
    {
        numevents++;
        current = current->next;
    }
    current = head;
    while (current)
    {
        progressbar(current->index, numevents);
        if (current->type == STEPRESPONSE)
        {
            status = step_response(current, risetime, maxiters, minstep);
            if (status != GSL_SUCCESS)
            {
                current->type = BADFIT;
            }
        }
        current = current->next;
    }
}
