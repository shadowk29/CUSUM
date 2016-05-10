#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED
#include<stdio.h>
#include<stdlib.h>
#include<inttypes.h>
#include<stdint.h>
#include<math.h>
#define EPS 1e-50
#define STRLENGTH 1024
#define HEAD -1000
#define NUMTYPES 12
#define CUSUM 0
#define STEPRESPONSE 1
#define BADBASELINE 2
#define TOOLONG 3
#define TOOSHORT 4
#define BADLEVELS 5
#define BADTRACE 6
#define BADFIT 7
#define FITRANGE 8
#define FITSIGN 9
#define FITSTEP 10
#define FITDIR 11

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

struct Bessel
{
    double *dcof;
    double *ccof;
    double cutoff;
    uint64_t order;
    uint64_t padding;
    double *paddedsignal;
    double *temp;
    double *tempback;
};
typedef struct Bessel bessel;

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
    //file reading parameters
    uint64_t start;
    uint64_t finish;
    uint64_t readlength;

    //filter parameters
    int usefilter;
    int eventfilter;
    double cutoff;
    uint64_t samplingfreq;
    uint64_t order; //must be even

    //detection parameters
    double threshold;
    double hysteresis;
    uint64_t event_minpoints;
    uint64_t event_maxpoints;

    double binsize;

    double baseline_min;
    double baseline_max;

    int event_direction;

    double cusum_min_threshold;
    double cusum_max_threshold;
    double cusum_delta;
    double cusum_minstep;
    uint64_t subevent_minpoints;

    uint64_t stepfit_samples;
    uint64_t maxiters;
    int attempt_recovery;
    int datatype;
    chimera *daqsetup;
};
typedef struct Configuration configuration;

FILE *fopen64_and_check(const char *fname, const char *mode, int error);
void *calloc_and_check(size_t num, size_t size, int error);
inline int signum(double num);
inline double my_min(double a, double b);
inline double my_max(double a, double b);
inline int64_t intmin(int64_t a, int64_t b);
inline int64_t intmax(int64_t a, int64_t b);
inline double d_abs(double num); //absolute value of a number

double ARL(uint64_t length, double sigma, double mun, double h);

uint64_t count_edges(edge *head_edge);
edge *initialize_edges(void);
edge *add_edge(edge *current, uint64_t location, int type);
void free_edges(edge *current);

event *initialize_events(void);
event *add_event(event *current, uint64_t start, uint64_t finish, uint64_t index);
void free_single_event(event *current);

cusumlevel *add_cusum_level(cusumlevel *lastlevel, double current, uint64_t length);
void free_levels(cusumlevel *current);
cusumlevel *initialize_levels(void);


uint64_t get_filesize(FILE *input, int datatype);
inline void progressbar(uint64_t pos, uint64_t finish, const char *msg);
#endif // UTILS_H_INCLUDED
