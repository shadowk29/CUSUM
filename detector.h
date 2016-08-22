#ifndef DETECTOR_H_INCLUDED
#define DETECTOR_H_INCLUDED
#include"utils.h"
#include"io.h"
#include"bessel.h"
#include<inttypes.h>
#include<stdint.h>
#include<string.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<float.h>
#define HISTOGRAM 0
#define FIRST_DERIV 1
#define SCND_DERIV 2
int64_t get_next_event(event *current_event, edge *current_edge, int64_t index);
int64_t get_next_event_start(edge *current_edge);

void calculate_level_noise(event *current, int64_t minpoints);
void populate_event_levels(event *current);
void identify_step_events(event *current, int64_t stepfit_samples, int64_t subevent_minpoints, int attempt_recovery);



void filter_short_events(event *current, int64_t minpoints);
void filter_long_events(event *current, int64_t event_maxpoints);
int64_t locate_min(double *signal, int64_t length);

int64_t average_cusum_levels(event *current, int64_t subevent_minpoints, double cusum_minstep, int attempt_recovery);
void cusum(event *current_event, double delta, double minthreshold, double maxthreshold, int64_t subevent_minpoints);
double get_cusum_threshold(int64_t length, double minthreshold, double maxthreshold, double sigma, double mun);

void refine_event_estimates(event *current);

void event_max_blockage(event *current);



edge *detect_edges(double *signal, double baseline, int64_t length, edge *current, double threshold, double hysteresis, int64_t position, int event_direction);
double build_histogram(double *signal, histostruct *histogram, int64_t length, double delta, double baseline_max, double baseline_min);
double baseline_averaging(double *signal, int64_t length, double baseline_min, double baseline_max);
void generate_trace(FILE *input, event *current, int datatype, void *rawsignal, FILE *logfile, bessel *lpfilter, int eventfilter, chimera *daqsetup, int64_t samplingfreq, edge *current_edge, int64_t last_end, int64_t start, int64_t subevent_minpoints);




void event_baseline(event *current_event, double baseline_min, double baseline_max);
void event_area(event *current_event, double timestep);
#endif // DETECTOR_H_INCLUDED
