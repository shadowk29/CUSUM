/*
 * Library:   lmfit (Levenberg-Marquardt least squares fitting)
 *
 * File:      lmstruct_int64.h
 *
 * Contents:  Declarations of parameter records, used in lmmin.h and lmcurve_int64.h
 *
 * Copyright: Joachim Wuttke, Forschungszentrum Juelich GmbH (2004-2013)
 *
 * License:   FreeBSD
 *
 * Homepage:  apps.jcns.fz-juelich.de/lmfit
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



#ifndef LMSTRUCT_INT64_H
#define LMSTRUCT_INT64_H
#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS /* empty */
#define __END_DECLS   /* empty */
#endif
__BEGIN_DECLS

#include <stdio.h>
#include<inttypes.h>

/* Collection of input parameters for fit control. */
typedef struct {
    double ftol;      /* Relative error desired in the sum of squares.
                         Termination occurs when both the actual and
                         predicted relative reductions in the sum of squares
                         are at most ftol. */
    double xtol;      /* Relative error between last two approximations.
                         Termination occurs when the relative error between
                         two consecutive iterates is at most xtol. */
    double gtol;      /* Orthogonality desired between fvec and its derivs.
                         Termination occurs when the cosine of the angle
                         between fvec and any column of the Jacobian is at
                         most gtol in absolute value. */
    double epsilon;   /* Step used to calculate the Jacobian, should be
                         slightly larger than the relative error in the
                         user-supplied functions. */
    double stepbound; /* Used in determining the initial step bound. This
                         bound is set to the product of stepbound and the
                         Euclidean norm of diag*x if nonzero, or else to
                         stepbound itself. In most cases stepbound should lie
                         in the interval (0.1,100.0). Generally, the value
                         100.0 is recommended. */
    int64_t patience;     /* Used to set the maximum number of function evaluations
                         to patience*(number_of_parameters+1). */
    int64_t scale_diag;   /* If 1, the variables will be rescaled internally.
                         Recommended value is 1. */
    FILE* msgfile;    /* Progress messages will be written to this file. */
    int64_t verbosity;    /* OR'ed: 1: print64_t some messages; 2: print64_t Jacobian. */
    int64_t n_maxpri;     /* -1, or max number of parameters to print. */
    int64_t m_maxpri;     /* -1, or max number of residuals to print. */
} lm_control_struct;

/* Collection of output parameters for status info. */
typedef struct {
    double fnorm;  /* norm of the residue vector fvec. */
    int64_t nfev;      /* actual number of iterations. */
    int64_t outcome;   /* Status indicator. Nonnegative values are used as index
                      for the message text lm_infmsg, set in lmmin.c. */
    int64_t userbreak; /* Set when function evaluation requests termination. */
} lm_status_struct;

/* Preset (and recommended) control parameter settings. */
extern const lm_control_struct lm_control_double;
extern const lm_control_struct lm_control_float;

/* Preset message texts. */

extern const char* lm_infmsg[];
extern const char* lm_shortmsg[];

__END_DECLS
#endif /* LMSTRUCT_H */
