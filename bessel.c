/*

                                COPYRIGHT
    Copyright (C) 2015 Kyle Briggs (kbrig035<at>uottawa.ca), and is based on
    the scipy.signal implementation of bessel

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

#include"bessel.h"

void filter_signal(double *signal, bessel *lpfilter, uint64_t length)
{
    uint64_t i;
    uint64_t p;
    uint64_t end;
    uint64_t order = lpfilter->order;
    uint64_t padding = lpfilter->padding;
    double *paddedsignal = lpfilter->paddedsignal;
    double *temp = lpfilter->temp;
    double *tempback = lpfilter->tempback;
    double *dcof = lpfilter->dcof;
    double *ccof = lpfilter->ccof;
    end = length+2*(order+padding);


    memcpy(&paddedsignal[order+padding],signal,length*sizeof(double));
    for (i=0; i<order+padding; i++)
    {
        temp[i] = signal_average(signal,padding);
        paddedsignal[i] = signal_average(signal,padding);
        paddedsignal[end-1-i]=signal_average(signal,padding);
        temp[end-1-i]=signal_average(signal,padding);
    }




    for (i=order; i<end; i++)
    {
        temp[i] = ccof[0]*paddedsignal[i];
        for (p=1; p<=order; p++)
        {
            temp[i] += ccof[p]*paddedsignal[i-p] - dcof[p]*temp[i-p];
        }
    }




    for (i=0; i<order+padding; i++)
    {
        tempback[end-1-i] = signal_average(&temp[order],padding);
        tempback[i] = signal_average(&temp[order],padding);
    }
    for (i=order; i<end; i++)
    {
        tempback[end-1-i] = ccof[0]*temp[end-1-i];
        for (p=1; p<=order; p++)
        {
            tempback[end-1-i] += ccof[p]*temp[end-1-i+p] - dcof[p]*tempback[end-1-i+p];
        }
    }




    memcpy(signal,&tempback[order+padding],length*sizeof(double));
    memset(paddedsignal,'0',(length+2*(order+padding))*sizeof(double));
    memset(temp,'0',(length+2*(order+padding))*sizeof(double));
    memset(tempback,'0',(length+2*(order+padding))*sizeof(double));
}


void besselap(uint64_t N, double complex *poles, double complex *zeros)
{
    switch(N)
    {
        case 2:
            poles[0] = -.8660254037844386467637229 + .4999999999999999999999996 * I;
            poles[1] = -.8660254037844386467637229 - .4999999999999999999999996 * I;
            zeros[0] = -1.0 + 0*I;
            zeros[1] = -1.0 + 0*I;
            break;
        case 4:
            poles[0] = -.6572111716718829545787781 - .8301614350048733772399715 * I;
            poles[1] = -.6572111716718829545787788 + .8301614350048733772399715 * I;
            poles[2] = -.9047587967882449459642637 - .2709187330038746636700923 * I;
            poles[3] = -.9047587967882449459642624 + .2709187330038746636700926 * I;
            zeros[0] = -1.0 + 0*I;
            zeros[1] = -1.0 + 0*I;
            zeros[2] = -1.0 + 0*I;
            zeros[3] = -1.0 + 0*I;
            break;
        case 6:
            poles[0] = -.9093906830472271808050953 - .1856964396793046769246397 * I;
            poles[1] = -.9093906830472271808050953 + .1856964396793046769246397 * I;
            poles[2] = -.7996541858328288520243325 - .5621717346937317988594118 * I;
            poles[3] = -.7996541858328288520243325 + .5621717346937317988594118 * I;
            poles[4] = -.5385526816693109683073792 - .9616876881954277199245657 * I;
            poles[5] = -.5385526816693109683073792 + .9616876881954277199245657 * I;
            zeros[0] = -1.0 + 0*I;
            zeros[1] = -1.0 + 0*I;
            zeros[2] = -1.0 + 0*I;
            zeros[3] = -1.0 + 0*I;
            zeros[4] = -1.0 + 0*I;
            zeros[5] = -1.0 + 0*I;
            break;
        case 8:
            poles[0] = -.9096831546652910216327629 - .1412437976671422927888150 * I;
            poles[1] = -.9096831546652910216327629 + .1412437976671422927888150 * I;
            poles[2] = -.8473250802359334320103023 - .4259017538272934994996429 * I;
            poles[3] = -.8473250802359334320103023 + .4259017538272934994996429 * I;
            poles[4] = -.7111381808485399250796172 - .7186517314108401705762571 * I;
            poles[5] = -.7111381808485399250796172 + .7186517314108401705762571 * I;
            poles[6] = -.4621740412532122027072175 - 1.034388681126901058116589 * I;
            poles[7] = -.4621740412532122027072175 + 1.034388681126901058116589 * I;
            zeros[0] = -1.0 + 0*I;
            zeros[1] = -1.0 + 0*I;
            zeros[2] = -1.0 + 0*I;
            zeros[3] = -1.0 + 0*I;
            zeros[4] = -1.0 + 0*I;
            zeros[5] = -1.0 + 0*I;
            zeros[6] = -1.0 + 0*I;
            zeros[7] = -1.0 + 0*I;
            break;
        case 10:
            poles[0] = -.9091347320900502436826431 - .1139583137335511169927714 * I;
            poles[1] = -.9091347320900502436826431 + .1139583137335511169927714 * I;
            poles[2] = -.8688459641284764527921864 - .3430008233766309973110589 * I;
            poles[3] = -.8688459641284764527921864 + .3430008233766309973110589 * I;
            poles[4] = -.7837694413101441082655890 - .5759147538499947070009852 * I;
            poles[5] = -.7837694413101441082655890 + .5759147538499947070009852 * I;
            poles[6] = -.6417513866988316136190854 - .8175836167191017226233947 * I;
            poles[7] = -.6417513866988316136190854 + .8175836167191017226233947 * I;
            poles[8] = -.4083220732868861566219785 - 1.081274842819124562037210 * I;
            poles[9] = -.4083220732868861566219785 + 1.081274842819124562037210 * I;
            zeros[0] = -1.0 + 0*I;
            zeros[1] = -1.0 + 0*I;
            zeros[2] = -1.0 + 0*I;
            zeros[3] = -1.0 + 0*I;
            zeros[4] = -1.0 + 0*I;
            zeros[5] = -1.0 + 0*I;
            zeros[6] = -1.0 + 0*I;
            zeros[7] = -1.0 + 0*I;
            zeros[8] = -1.0 + 0*I;
            zeros[9] = -1.0 + 0*I;
            break;
    }
}

void polymult(double complex *roots, uint64_t N, double *rcoefs)
{
    uint64_t i, j;
    double complex *polycoefs;
    polycoefs = calloc_and_check(N+1, sizeof(double complex),"Cannot allocate polynomial coefficients");
    polycoefs[N]=1.0;
    for (i=0; i<N; i++)
    {
        for (j=0; j<N+1; j++)
        {
            polycoefs[j] = -roots[i]*polycoefs[j] + (j == N ? 0 : polycoefs[j+1]);
        }
    }

    for (i=0; i<N+1; i++)
    {
        rcoefs[i] = creal(polycoefs[i]);
    }
    free(polycoefs);
}

double scale_filter(double complex *poles, uint64_t N, double warped, double scale)
{
    uint64_t i;
    for (i=0; i<N; i++)
    {
        poles[i] *= warped;
    }
    return scale * pow(warped, (double) N);
}

double bilinear(double complex *poles, uint64_t N, double scale, double fs)
{
    double fs2 = 2.0 * fs;
    double complex poleprod = 1;
    uint64_t i;
    for (i=0; i<N; i++)
    {
        poleprod *= fs2 - poles[i];
        poles[i] = (fs2 + poles[i])/(fs2-poles[i]);
    }
    return scale * creal(1.0 / poleprod);
}

void transform_filter(double complex *poles, double complex *zeros, uint64_t N, double scale, double *b, double *a)
{
    uint64_t i;
    polymult(zeros, N, b);
    polymult(poles, N, a);
    for (i=0; i<N+1; i++)
    {
        b[i] *= scale;
    }
}

bessel *initialize_filter(bessel *lpfilter, uint64_t order, double cutoff, uint64_t length, uint64_t samplingfreq)
{



    if ((lpfilter = malloc(sizeof(bessel)))==NULL)
    {
        printf("Cannot allocate filter memory\n");
        exit(36);
    }



    lpfilter->order = order;
    lpfilter->cutoff = cutoff;
    lpfilter->padding = (uint64_t) (100e-6*samplingfreq);


    double complex *poles;
    double complex *zeros;

    poles = calloc_and_check(order,sizeof(double complex),"Cannot allocate Bessel poles");

    zeros = calloc_and_check(order,sizeof(double complex),"Cannot allocate Bessel zeros");

    lpfilter->ccof = calloc_and_check(order+1,sizeof(double complex),"Cannot allocate ccof");

    lpfilter->dcof = calloc_and_check(order+1,sizeof(double complex),"Cannot allocate dcof");


    double fs = 2.0;
    double warped = 2.0 * fs * tan(M_PI * cutoff / fs);
    double scale = 1;

    besselap(order, poles, zeros);
    scale = scale_filter(poles,order, warped, scale);
    scale = bilinear(poles, order, scale, fs);
    transform_filter(poles, zeros, order, scale, lpfilter->ccof, lpfilter->dcof);

    lpfilter->temp = NULL;
    lpfilter->tempback = NULL;
    lpfilter->paddedsignal = NULL;

    lpfilter->temp = calloc_and_check(length+2*(order+lpfilter->padding), sizeof(double),"Cannot allocate filter temp");

    lpfilter->tempback = calloc_and_check(length+2*(order+lpfilter->padding), sizeof(double),"Cannot allocate filter tempback");

    lpfilter->paddedsignal = calloc_and_check(length+2*(order+lpfilter->padding), sizeof(double),"Cannot allocate paddedsignal");

    free(poles);
    free(zeros);
    return lpfilter;
}

void free_filter(bessel *lpfilter)
{
    free(lpfilter->dcof);
    free(lpfilter->ccof);
    free(lpfilter->temp);
    free(lpfilter->tempback);
    free(lpfilter->paddedsignal);
    free(lpfilter);
}


