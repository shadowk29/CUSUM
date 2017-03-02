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

#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED
#include<stdio.h>
#include<stdlib.h>
#include<inttypes.h>
#include<stdint.h>
#include<math.h>
#include<string.h>
#include<time.h>
#include<limits.h>
#define EPS 1.2e-16
#define STRLENGTH 1024
#define HEAD -1000
#define NUMTYPES 10

#define SECONDS_TO_MICROSECONDS 1e6
#define AMPS_TO_PICOAMPS 1e12
#define FRACTION_TO_PERCENTAGE 100

#define CUSUM 0
#define STEPRESPONSE 1
#define BADBASELINE 2
#define TOOLONG 3
#define TOOSHORT 4
#define BADLEVELS 5
#define BADTRACE 6
#define BADPADDING 7
#define FITERR 8
#define FITSTEP 9

//#define DEBUG

struct Chimera
{
    double samplerate;
    double TIAgain;
    double preADCgain;
    double currentoffset;
    double ADCvref;
    int ADCbits;
};
typedef struct Chimera chimera;


struct Cusumlevel
{
    double current;
    double stdev;
    int64_t length;
    struct Cusumlevel *next;
};
typedef struct Cusumlevel cusumlevel;


struct Baseline_struct
{
    double *histogram;
    double *current;
    int64_t numbins;
    double baseline_min;
    double baseline_max;
    double range;
    double delta;
    double mean;
    double stdev;
    double amplitude;
};
typedef struct Baseline_struct baseline_struct;

struct Event
{
    int64_t index;
    int64_t start;
    int64_t finish;
    int64_t length;
    int type;
    double area;
    double baseline_before;
    double baseline_after;
    double average_blockage;
    double max_blockage;
    int64_t max_length;
    double *signal;
    double *paddedsignal;
    double *filtered_signal;
    int64_t padding_before;
    int64_t padding_after;
    int numlevels;
    double threshold;
    double rc1;
    double rc2;
    double residual;
    double maxdeviation;
    double local_stdev;
    struct Edge *first_edge;
    struct Cusumlevel *first_level;
};
typedef struct Event event;

struct Edge
{
    int64_t location;
    int64_t type;
    double local_stdev;
    struct Edge *next;
};
typedef struct Edge edge;


struct Configuration
{
    char filepath[STRLENGTH]; //input file
    char outputfolder[STRLENGTH];
    char eventsfolder[STRLENGTH];
    char eventsfile[STRLENGTH];
    char ratefile[STRLENGTH];
    char logfile[STRLENGTH];
    char baselinefile[STRLENGTH];

    //file reading parameters
    int64_t start;
    int64_t finish;
    int64_t readlength;

    //filter parameters
    int usefilter;
    int eventfilter;
    double cutoff;
    int64_t samplingfreq;
    int64_t order; //must be even

    //detection parameters
    double threshold;
    double hysteresis;
    int64_t event_minpoints;
    int64_t event_maxpoints;

    double baseline_min;
    double baseline_max;

    int event_direction;

    double cusum_min_threshold;
    double cusum_max_threshold;
    double cusum_delta;
    double cusum_minstep;
    int64_t subevent_minpoints;

    int64_t stepfit_samples;
    int64_t maxiters;
    int attempt_recovery;
    int datatype;
    chimera *daqsetup;
};
typedef struct Configuration configuration;

FILE *fopen64_and_check(const char *fname, const char *mode, int error);
void *calloc_and_check(int64_t num, int64_t size, char *msg);
int signum(double num);
double my_min(double a, double b);
double my_max(double a, double b);
int64_t intmin(int64_t a, int64_t b);
int64_t intmax(int64_t a, int64_t b);
double d_abs(double num); //absolute value of a number

double ARL(int64_t length, double sigma, double mun, double h);

int64_t count_edges(edge *head_edge);
edge *initialize_edges(void);
edge *add_edge(edge *current, int64_t location, int type, double stdev);
void free_edges(edge *current);

event *initialize_events(void);
event *add_event(event *current, int64_t start, int64_t finish, int64_t index, double local_stdev);
void free_single_event(event *current);

cusumlevel *add_cusum_level(cusumlevel *lastlevel, double current, int64_t length);
void free_levels(cusumlevel *current);
cusumlevel *initialize_levels(void);

double signal_max(double *signal, int64_t length);
double signal_min(double *signal, int64_t length);
double signal_average(double *signal, int64_t length);
double signal_extreme(double *signal, int64_t length, double sign);
double signal_variance(double *signal, int64_t length);


int64_t get_filesize(FILE *input, int datatype);
void progressbar(int64_t pos, int64_t finish, const char *msg, double elapsed);



void invert_matrix(double m[3][3], double inverse[3][3]);
baseline_struct *initialize_baseline(baseline_struct *baseline, configuration *config);
void free_baseline(baseline_struct *baseline);

int64_t locate_min(double *signal, int64_t length);
int64_t locate_max(double *signal, int64_t length);
#endif // UTILS_H_INCLUDED
