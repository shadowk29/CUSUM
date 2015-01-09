#ifndef DETECTOR_H_INCLUDED
#define DETECTOR_H_INCLUDED
#include"utils.h"
#include"io.h"
#define HISTOGRAM 0
#define FIRST_DERIV 1
#define SCND_DERIV 2

void count_all_levels(event *current);
int count_levels(event *current);
void detect_subevents(event *current_event, double delta, double threshold);
uint64_t locate_min(double *signal, uint64_t length);

void assign_cusum_levels(event *current);
void average_cusum_levels(event *current);
void cusum(event *current_event, double delta, double threshold);



void populate_event_traces(FILE *input, event *current_event, uint64_t order);


event *process_edges(edge *current_edge, event *current_event);
edge *detect_edges(double *signal, double baseline, uint64_t length, edge *current, double threshold, double hysteresis, uint64_t position, int event_direction);
double build_histogram(double *signal, histostruct *histogram, uint64_t length, double delta, uint64_t pos, double baseline_max, double baseline_min);

void generate_trace(FILE *input, event *current, uint64_t order);



double signal_max(double *signal, uint64_t length);
double signal_min(double *signal, uint64_t length);
double signal_average(double *signal, uint64_t length);
double signal_extreme(double *signal, uint64_t length, double sign);
double signal_variance(double *signal, uint64_t length);

void event_baseline(event *current_event);
void assign_event_baselines(event *current_event);
void event_area(event *current_event, double timestep);
void assign_event_areas(event *current_event, double timestep);
#endif // DETECTOR_H_INCLUDED
