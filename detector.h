#ifndef DETECTOR_H_INCLUDED
#define DETECTOR_H_INCLUDED
#include"utils.h"
#include"io.h"
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
uint64_t get_next_event(event *current_event, edge *current_edge, uint64_t index);

void calculate_level_noise(event *current, uint64_t minpoints);
void populate_event_levels(event *current);


void filter_signal(double *signal, double *filtered, bessel *lpfilter, uint64_t length);

void filter_event_length(event *current, uint64_t maxpoints, uint64_t minpoints, uint64_t stepfit_samples);
uint64_t locate_min(double *signal, uint64_t length);

uint64_t average_cusum_levels(event *current, uint64_t subevent_minpoints, double cusum_minstep, int attempt_recovery);
void cusum(event *current_event, double delta, double minthreshold, double maxthreshold, uint64_t subevent_minpoints);
double get_cusum_threshold(uint64_t length, double minthreshold, double maxthreshold, double sigma, double mun);

void refine_event_estimates(event *current);

void event_max_blockage(event *current);



edge *detect_edges(double *signal, double baseline, uint64_t length, edge *current, double threshold, double hysteresis, uint64_t position, int event_direction);
double build_histogram(double *signal, histostruct *histogram, uint64_t length, double delta, double baseline_max, double baseline_min);
double baseline_averaging(double *signal, uint64_t length, double baseline_min, double baseline_max);
void generate_trace(FILE *input, event *current, int datatype, FILE *logfile, bessel *lpfilter, int eventfilter, chimera *daqsetup, uint64_t samplingfreq);


double signal_max(double *signal, uint64_t length);
double signal_min(double *signal, uint64_t length);
double signal_average(double *signal, uint64_t length);
double signal_extreme(double *signal, uint64_t length, double sign);
double signal_variance(double *signal, uint64_t length);

void event_baseline(event *current_event, double baseline_min, double baseline_max);
void event_area(event *current_event, double timestep);
#endif // DETECTOR_H_INCLUDED
