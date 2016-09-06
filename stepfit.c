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

double stepfunc(double time, const double *p)
{
    return p[0] - p[1]*heaviside(time-p[2])*(1.0-exp(-(time-p[2])/p[3])) + p[4]*heaviside(time-p[5])*(1.0-exp(-(time-p[5])/p[6]));
}

void time_array(double *time, int64_t m)
{
    int64_t i;
    for (i=0; i<m; i++)
    {
        time[i] = i;
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

        par[0] = baseline;
        par[1] = (sign > 0 ? maxsignal - minsignal: minsignal - maxsignal);
        par[2] = start - (int64_t) ( risetime);
        par[3] = risetime;
        par[4] = (sign > 0 ? maxsignal - minsignal: minsignal - maxsignal);
        par[5] = end - (int64_t) (risetime);
        par[6] = risetime;

        int64_t i;

        lm_control_struct control = lm_control_double;
        control.patience = maxiters;
        lm_status_struct status;

        lmcurve_int64(n, par, length, time, current->signal, stepfunc, &control, &status );


        if (status.outcome < 1 || status.outcome > 3)
        {
            current->type = -10-status.outcome;
            return;
        }

        double i0 = par[0];
        double a = par[1];
        int64_t u1 = par[2];
        double rc1 = par[3];
        double b = par[4];
        int64_t u2 = par[5];
        double rc2 = par[6];
        double residual = 0;

        //printf("i0 = %g\na=%g\nu1=%"PRId64"\nb=%g\nu2=%"PRId64"\n",i0,a,u1,b,u2);

        if (u1 > u2) //if the fit got stuck in a variable swapped equivalent minimum, we can just reverse the paramters
        {
            int64_t temp;
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

        if (u2 >= length || u1 <= 0) //if for some reason we are out of range
        {
            current->type = FITRANGE;
            return;
        }
        else if (signum(a) != signum(b)) //if it is not a step-and-return event
        {
            current->type = FITSIGN;
            return;
        }
        else if (d_abs(b) < minstep || d_abs(b) < minstep)
        {
            current->type = FITSTEP;
            return;
        }
        else if (signum(a) != sign || signum(b) != sign)

        {
            current->type = FITDIR;
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

