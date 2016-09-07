/*

                                COPYRIGHT
    Copyright (C) 2015 Kyle Briggs (kbrig035<at>uottawa.ca)

    This section of the code implements an algorithm used in MOSAIC,
    which can be found here: http://usnistgov.github.io/mosaic/html/

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

double heaviside(double x)
{
    return x>0 ? 1 : 0;
}

double stepfunc(double time, const double *p, double maxlength, double maxstep, double risetime, int sign)
{
    double sigma1 = risetime * exp(p[3]);
    double sigma2 = risetime * exp(p[6]);
    double t1 = maxlength/2.0 * (1.0 + tanh(p[2]));
    double t2 = maxlength/2.0 * (1.0 + tanh(p[5]));
    double a = sign*maxstep/2.0 * (1.0 + tanh(p[1]));
    double b = sign*maxstep/2.0 * (1.0 + tanh(p[4]));
    return p[0] - a*heaviside(time-t1)*(1.0-exp(-(time-t1)/sigma1)) + b*heaviside(time-t2)*(1.0-exp(-(time-t2)/sigma2));
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
    fvec[i] = D->signal[i] - D->stepfunc(D->time[i], p, D->maxlength, D->maxstep, D->risetime, D->sign);
}

void step_response(event *current, double risetime, int64_t maxiters, double minstep)
{
#ifdef DEBUG
    printf("StepResponse\n");
    fflush(stdout);
#endif // DEBUG
    if (current->type == STEPRESPONSE)
    {
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


        data_struct data = {time, current->signal, maxlength, maxstep, risetime, sign, stepfunc};

        par[0] = baseline;
        par[1] = atanh(2.0*stepguess/maxstep-1.0);
        par[2] = atanh(2.0*(start - (int64_t) risetime)/maxlength - 1.0);
        par[3] = 0;
        par[4] = atanh(2.0*stepguess/maxstep-1.0);
        par[5] = atanh(2.0*(end - (int64_t) risetime)/maxlength - 1.0);
        par[6] = 0;

        int64_t i;

        lm_control_struct control = lm_control_double;
        control.patience = maxiters;
        lm_status_struct status;

        lmmin_int64(n, par, length, (const void*) &data, evaluate, &control, &status );


        if (status.outcome < 1 || status.outcome > 3)
        {
            if (status.outcome != 12)
            {
                current->type = FITERR;
            }
            else
            {
                current->type = FITNAN;
            }
            return;
        }

        double i0 = par[0];
        double a = sign*maxstep/2.0 * (1.0 + tanh(par[1]));
        int64_t u1 = maxlength/2.0 * (1.0 + tanh(par[2]));
        double rc1 = risetime * exp(par[3]);
        double b = sign*maxstep/2.0 * (1.0 + tanh(par[4]));
        int64_t u2 = maxlength/2.0 * (1.0 + tanh(par[5]));
        double rc2 = risetime * exp(par[6]);
        double residual = 0;

        if (d_abs(b) < minstep || d_abs(b) < minstep)
        {
            current->type = FITSTEP;
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

