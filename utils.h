#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED
#include<stdio.h>
#include<stdlib.h>
#include<inttypes.h>
#include<stdint.h>
#include<math.h>
#include"butterworth.h"
#define EPS 1e-50
#define STRLENGTH 1024
#define HEAD -1000
#define STEPRESPONSE 1
#define CUSUM 0
#define BADBASELINE -1
#define TOOLONG -2
#define TOOSHORT -3
#define BADLEVELS -4
#define BADTRACE -5
#define BADFIT -6
#define FITRANGE -7
#define FITSIGN -8
#define FITSTEP -9
#define FITDIR -10

struct Butterworth
{
    double *dcof;
    int *ccof;
    double scale;
    uint64_t order;
    double *paddedsignal;
    double *temp;
    double *tempback;
};
typedef struct Butterworth butterworth;

struct Cusumlevel
{
    double current;
    double stdev;
    uint64_t length;
    struct Cusumlevel *next;
};
typedef struct Cusumlevel cusumlevel;

struct Histostruct
{
    double **histogram;
    uint64_t numbins;
    double offset;
    double delta;
};
typedef struct Histostruct histostruct;

struct Event
{
    int64_t index;
    uint64_t start;
    uint64_t finish;
    uint64_t length;
    int type;
    double area;
    double baseline_before;
    double baseline_after;
    double average_blockage;
    double max_blockage;
    uint64_t max_length;
    double *signal;
    double *filtered_signal;
    double binsize;
    uint64_t padding_before;
    uint64_t padding_after;
    int numlevels;
    double threshold;
    double rc1;
    double rc2;
    double residual;
    struct Edge *first_edge;
    struct Cusumlevel *first_level;
    struct Event *next;
    struct Event *prev;
};
typedef struct Event event;

struct Edge
{
    uint64_t location;
    int64_t type;
    struct Edge *next;
};
typedef struct Edge edge;


struct Configuration
{
    char filepath[STRLENGTH]; //input file
    char tracefile[STRLENGTH];
    //file reading parameters
    uint64_t start;
    uint64_t finish;
    uint64_t readlength;

    //filter parameters
    int usefilter;
    double cutoff;
    uint64_t samplingfreq;
    uint64_t order; //must be even

    //detection parameters
    double threshold;
    double hysteresis;

    double binsize;

    uint64_t event_minpoints;
    uint64_t event_maxpoints;
    uint64_t subevent_minpoints;

    uint64_t export_trace_start;
    uint64_t export_trace_end;

    double baseline_min;
    double baseline_max;

    int event_direction;

    double cusum_min_threshold;
    double cusum_max_threshold;
    double cusum_delta;
    double cusum_minstep;
    int refine_estimates;
    uint64_t stepfit_samples;
    uint64_t maxiters;
    int attempt_recovery;

    int datatype;

};
typedef struct Configuration configuration;


inline int signum(double num);
inline double my_min(double a, double b);
inline double my_max(double a, double b);
inline int64_t intmin(int64_t a, int64_t b);
inline int64_t intmax(int64_t a, int64_t b);
inline double d_abs(double num); //absolute value of a number

double ARL(uint64_t length, double sigma, double mun, double h);


edge *initialize_edges(void);
edge *add_edge(edge *current, uint64_t location, int type);
void free_edges(edge *current);

event *initialize_events(void);
event *add_event(event *current, uint64_t start, uint64_t finish);
void free_events(event *current);
void free_single_event(event *current);
event *delete_bad_events(event *head);

cusumlevel *add_cusum_level(cusumlevel *lastlevel, double current, uint64_t length);
void free_levels(cusumlevel *current);
cusumlevel *initialize_levels(void);


uint64_t get_filesize(FILE *input, int datatype);
inline void progressbar(uint64_t pos, uint64_t finish);

butterworth *initialize_filter(butterworth *lpfilter, uint64_t order, double cutoff, uint64_t length);
void free_filter(butterworth *lpfilter);
#endif // UTILS_H_INCLUDED
