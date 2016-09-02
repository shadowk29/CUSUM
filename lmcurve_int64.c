/*
 * Library:   lmfit (Levenberg-Marquardt least squares fitting)
 *
 * File:      lmcurve_int64.c
 *
 * Contents:  Implements lmcurve_int64, a simplified API for curve fitting
 *            using the generic Levenberg-Marquardt routine lmmin.
 *
 * Copyright: Joachim Wuttke, Forschungszentrum Juelich GmbH (2004-2013)
 *
 * License:   FreeBSD
 *
 * Homepage:  apps.jcns.fz-juelich.de/lmfit
 *
 * Note to programmers: Don't patch and fork, but copy and variate!
 *   If you need to compute residues differently, then please do not patch
 * lmcurve.c, but copy it to a differently named file, and rename lmcurve(),
 * lmcurve_evaluate() and lmcurve_data_struct before adapting them to your
 * needs, like we have done in lmcurve_tyd.c.
 */

/*
 * Copyright (c) Joachim Wuttke, Forschungszentrum Juelich GmbH (2004-2013)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation and/or
 * other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 * OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
 * WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * The views and conclusions contained in the software and documentation are those of the
 * authors and should not be interpreted as representing official policies, either expressed
 * or implied, of the FreeBSD Project.
 */

 /*
 * Changelog:   KB      Aug 5, 2016     Updated all instances of int to int64_t to handle long data sets
 *
 *
 *
 */

#include "lmmin_int64.h"

typedef struct {
    const double* t;
    const double* y;
    double (*f)(const double t, const double* par);
} lmcurve_int64_data_struct;

void lmcurve_int64_evaluate(
    const double* par, const int64_t m_dat, const void* data, double* fvec,
    int64_t* info)
{
    info = info;
    lmcurve_int64_data_struct* D = (lmcurve_int64_data_struct*)data;
    int64_t i;
    for (i = 0; i < m_dat; i++)
        fvec[i] = D->y[i] - D->f(D->t[i], par);
}

void lmcurve_int64(
    const int64_t n_par, double* par, const int64_t m_dat,
    const double* t, const double* y,
    double (*f)(const double t, const double* par),
    const lm_control_struct* control, lm_status_struct* status)
{
    lmcurve_int64_data_struct data = { t, y, f };
    lmmin_int64(n_par, par, m_dat, (const void*)&data, lmcurve_int64_evaluate,
          control, status);
}
