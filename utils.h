#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED
#include<stdio.h>
#include<stdlib.h>
#include<inttypes.h>
#include<stdint.h>
#define EPS 1e-50
#define STRLENGTH 1024
#define HEAD -1000

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
    double *signal;
    double *filtered_signal;
    double binsize;
    uint64_t padding;
    struct Edge *first_edge;
    struct Event *next;
    struct Event *prev;
};
typedef struct Event event;

struct Edge
{
    uint64_t location;
    int type;
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
    uint64_t poles; //must be even

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

    double cusum_threshold;
    double cusum_delta;

};
typedef struct Configuration configuration;


int signum(double num);
double my_min(double a, double b);
double my_max(double a, double b);
int64_t intmin(int64_t a, int64_t b);
int64_t intmax(int64_t a, int64_t b);
double d_abs(double num); //absolute value of a number

edge *initialize_edges(void);
edge *add_edge(edge *current, uint64_t location, int type);
void free_edges(edge *current);

event *initialize_events(void);
event *add_event(event *current, uint64_t start, uint64_t finish);
void free_events(event *current);
#endif // UTILS_H_INCLUDED
