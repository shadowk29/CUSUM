/*

                                COPYRIGHT
    Copyright (C) 2015-2016 Kyle Briggs (kbrig035<at>uottawa.ca)

    This section of the code implements an algorithm used in MOSAIC,
    which can be found here: http://usnistgov.github.io/mosaic/html/
    and takes advantage of a modified version of lmfit-6.1, which
    can be found here: apps.jcns.fz-juelich.de/lmfit

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

double stepfunc(double time, const double *p, double maxlength, double maxstep, double maxbaseline, double risetime, int sign)
{
    double sigma1, sigma2, a, b;
    double t1 = maxlength/2.0 * (1.0 + tanh(p[2])); //constrain t1 to (0, maxlength)
    double t2 = maxlength/2.0 * (1.0 + tanh(p[5])); //constrint t2 to (0, maxlength)
    double fitval = sign*maxbaseline/2.0 * (1.0 + tanh(p[0])); //constrain baseline to sign * (0, maxbaseline)
    if (time > t1)
    {
        sigma1 = risetime * exp(p[3]); //constrain sigma1 to (0, inf)
        a = sign*maxstep/2.0 * (1.0 + tanh(p[1])); //constrain a to sign * (0, maxstep)
        fitval -= a*(1.0-exp(-(time-t1)/sigma1));
    }
    if (time > t2)
    {
        sigma2 = risetime * exp(p[6]); //constrain sigma2 to (0, inf)
        b = sign*maxstep/2.0 * (1.0 + tanh(p[4])); //constrain b to sign * (0, maxstep)
        fitval += b*(1.0-exp(-(time-t2)/sigma2));
    }
    return fitval;
}

void time_array(double *time, int64_t m)
{
    int64_t i;
    for (i=0; i<m; i++)
    {
        time[i] = i;
    }
}

void evaluate(const double *p, int64_t length, const void *data, double *fvec, int64_t *userbreak)
{
    data_struct *D;
    D = (data_struct*)data;
    *userbreak=0;
    int64_t i;
    for (i=0; i<length; i++)
    {
        fvec[i] = D->signal[i] - D->stepfunc(D->time[i], p, D->maxlength, D->maxstep, D->maxbaseline, D->risetime, D->sign);
    }
}

void step_response(event *current, double risetime, int64_t maxiters, double minstep)
{
#ifdef DEBUG
    printf("StepResponse\n");
    fflush(stdout);
#endif // DEBUG
    if (current->type == STEPRESPONSE)
    {
        minstep *= current->local_stdev;
        int64_t length = current->length + current->padding_before + current->padding_after; //number of data points
        double *time;
        time = calloc_and_check(length, sizeof(double), "cannot allocate stepfit time array"); //time array
        time_array(time, length);
        int64_t n = 7; // number of parameters in model function f
        double par[n];  // parameter array


        double maxsignal = signal_max(current->signal, current->length + current->padding_before + current->padding_after);
        double minsignal = signal_min(current->signal, current->length + current->padding_before + current->padding_after);
        double baseline = signal_average(current->signal,current->padding_before);
        int sign = signum(baseline);
        int64_t start = current->padding_before;
        int64_t end = current->length + current->padding_before;
        double stepguess = 0.66 * sign * (minsignal - maxsignal);
        double maxlength = (double) length;
        double maxstep = d_abs(maxsignal - minsignal);
        double maxbaseline = my_max(d_abs(maxsignal),d_abs(minsignal));


        data_struct data = {time, current->signal, maxlength, maxstep, maxbaseline, risetime, sign, stepfunc};

        par[0] = atanh(2.0*sign*baseline/maxbaseline-1.0);
        par[1] = atanh(2.0*stepguess/maxstep-1.0);
        par[2] = atanh(2.0*(start - (int64_t) risetime)/maxlength - 1.0);
        par[3] = 0;
        par[4] = atanh(2.0*stepguess/maxstep-1.0);
        par[5] = atanh(2.0*(end - (int64_t) risetime)/maxlength - 1.0);
        par[6] = 0;

        int64_t i;

        lm_control_struct control = lm_control_float;
        control.patience = maxiters;
        lm_status_struct status = {0,0,0,0};

        lmmin_int64(n, par, length, (const void*) &data, evaluate, &control, &status );


        if (status.outcome < 1 || status.outcome > 3)
        {
            current->type = status.outcome + 10;
            free(time);
            return;
        }

        double i0 = sign*maxbaseline/2.0 * (1.0 + tanh(par[0]));
        double a = sign*maxstep/2.0 * (1.0 + tanh(par[1]));
        int64_t u1 = maxlength/2.0 * (1.0 + tanh(par[2]));
        double rc1 = risetime * exp(par[3]);
        double b = sign*maxstep/2.0 * (1.0 + tanh(par[4]));
        int64_t u2 = maxlength/2.0 * (1.0 + tanh(par[5]));
        double rc2 = risetime * exp(par[6]);
        double residual = 0;

        if (d_abs(a) < minstep || d_abs(b) < minstep)
        {
            current->type = FITSTEP;
            printf("\n%"PRId64":\t%"PRId64"\t%g\t%g\t%g\n",current->index, current->start, a, b, minstep);
            free(time);
            return;
        }

        double t;
        for (i=0; i<u1; i++) //if all went well, populate the filtered trace
        {
            current->filtered_signal[i] = i0;
            residual += (current->signal[i]-i0)*(current->signal[i]-i0);
        }
        for (i=(int64_t) u1; i<u2; i++)
        {
            t = i;
            current->filtered_signal[i] = i0-a;
            residual += (current->signal[i]-(i0+a*(exp(-(t-u1)/rc1)-1)))*(current->signal[i]-(i0+a*(exp(-(t-u1)/rc1)-1)));
        }
        for (i=(int64_t) u2; i<length; i++)
        {
            t = i;
            current->filtered_signal[i] = i0-a+b;
            residual += (current->signal[i]-(i0+a*(exp(-(t-u1)/rc1)-1)+b*(1-exp(-(t-u2)/rc2))))*(current->signal[i]-(i0+a*(exp(-(t-u1)/rc1)-1)+b*(1-exp(-(t-u2)/rc2))));
        }
        current->rc1 = rc1;
        current->rc2 = rc2;
        current->residual = sqrt(residual/(current->length + current->padding_before + current->padding_after));
        free(time);
    }
}

