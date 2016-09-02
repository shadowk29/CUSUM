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
#include"lmcurve_int64.h"
#include"lmmin_int64.h"

double heaviside(double x);
double stepfunc(double time, const double *p);
void time_array(double *time, double timestep, int64_t m);


void step_response(event *current, double risetime, int64_t maxiters, double minstep, double timestep);

#endif //STEPFIT_H_INCLUDED
