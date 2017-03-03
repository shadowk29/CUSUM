/*

                                COPYRIGHT
    Copyright (C) 2015-2016 Kyle Briggs (kbrig035<at>uottawa.ca)

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

#ifndef STEPFIT_H_INCLUDED
#define STEPFIT_H_INCLUDED
#include <stdlib.h>
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include<math.h>
#include"utils.h"
#include"detector.h"
#include"lmstruct_int64.h"
#include"lmmin_int64.h"

struct Data_Struct
{
    long double *time;
    double *signal;
    long double maxlength;
    long double maxstep;
    long double maxbaseline;
    long double risetime;
    int sign;
    long double (*stepfunc)(long double time, const long double *p, long double maxlength, long double maxstep, long double maxbaseline, long double risetime, int sign);
};
typedef struct Data_Struct data_struct;

long double stepfunc(long double time, const long double *p, long double maxlength, long double maxstep,long  double maxbaseline, long double risetime, int sign);
void time_array(long double *time, int64_t m);
void step_response(event *current, long double risetime, int64_t maxiters, long double minstep);
void evaluate(const long double *p, int64_t length, const void *data, long double *fvec, int64_t *userbreak);

#endif //STEPFIT_H_INCLUDED
