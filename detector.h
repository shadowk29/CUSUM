#ifndef DETECTOR_H_INCLUDED
#define DETECTOR_H_INCLUDED
#include"utils.h"
#include"io.h"
#define HISTOGRAM 0
#define FIRST_DERIV 1
#define SCND_DERIV 2
#define BADBASELINE -1
#define TOOLONG -2
#define TOOSHORT -3
#define BADLEVELS -4
#define BADTRACE -5


void populate_all_levels(event *current);
void populate_event_levels(event *current);


void print_error_summary(event *current, FILE *logfile);

void filter_event_length(event *current, uint64_t maxpoints, uint64_t minpoints, FILE *logfile);
void detect_subevents(event *current_event, double delta, double minthreshold, double maxthreshold, uint64_t subevent_minpoints);
uint64_t locate_min(double *signal, uint64_t length);

void assign_cusum_levels(event *current, uint64_t subevent_minpoints, double cusum_minstep);
void average_cusum_levels(event *current, uint64_t subevent_minpoints, double cusum_minstep);
void cusum(event *current_event, double delta, double minthreshold, double maxthreshold, uint64_t subevent_minpoints);
double get_cusum_threshold(uint64_t length, double minthreshold, double maxthreshold, double sigma, double mun);

void refine_all_estimates(event *current);
void refine_event_estimates(event *current);

void find_max_blockages(event *current);
void event_max_blockage(event *current);

void populate_event_traces(FILE *input, event *current_event, int datatype, FILE *logfile);


event *process_edges(edge *current_edge, event *current_event);
edge *detect_edges(double *signal, double baseline, uint64_t length, edge *current, double threshold, double hysteresis, uint64_t position, int event_direction);
double build_histogram(double *signal, histostruct *histogram, uint64_t length, double delta, double baseline_max, double baseline_min);

void generate_trace(FILE *input, event *current, int datatype, FILE *logfile);



double signal_max(double *signal, uint64_t length);
double signal_min(double *signal, uint64_t length);
double signal_average(double *signal, uint64_t length);
double signal_extreme(double *signal, uint64_t length, double sign);
double signal_variance(double *signal, uint64_t length);

void event_baseline(event *current_event, double baseline_min, double baseline_max);
void assign_event_baselines(event *current_event, FILE *logfile, double baseline_min, double baseline_max);
void event_area(event *current_event, double timestep);
void assign_event_areas(event *current_event, double timestep);
#endif // DETECTOR_H_INCLUDED
