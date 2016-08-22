#ifndef IO_H_INCLUDED
#define IO_H_INCLUDED
#include"utils.h"
#include"detector.h"
#include<string.h>
#include<inttypes.h>
#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>

union doublebits
{
    double currentval;
    uint64_t bits;
};

void configure_defaults(configuration *config);
void config_sanity_check(configuration *config, FILE *logfile);
void print_error_summary(FILE *logfile, int64_t *error_summary, int64_t numevents);
void initialize_events_file(FILE *events, FILE *rate);
void print_event_line(FILE *events, FILE *rate, event *current, double timestep, int64_t );


//void print_events(event *current, double timestep);
void print_histogram(char *filename, histostruct *histogram);
FILE *read_config(configuration *config);
void swapByteOrder(double *current, uint64_t *rawsignal, int64_t length);
void swapByteOrder_int16(double *current, uint16_t *rawsignal, int64_t length);
void chimera_gain(double *current, uint16_t *rawsignal, int64_t length, chimera *daqsetup);
int64_t read_current_chimera(FILE *input, double *current, uint16_t *rawsignal, int64_t position, int64_t length, chimera *daqsetup);
int64_t read_current_double(FILE *input, double *current, uint64_t *rawsignal, int64_t position, int64_t length);
int64_t read_current_int16(FILE *input, double *current, uint16_t *rawsignal, int64_t position, int64_t length);
int64_t read_current(FILE *input, double *signal, void *rawsignal, int64_t position, int64_t length, int datatype, chimera *daqsetup);
void print_event_signal(int64_t index, event *current, double timestep, char *eventsfolder);
void print_signal(event *current, int64_t length, char *filename, double timestep);
//void print_error_summary(event *current, FILE *logfile);
#endif // IO_H_INCLUDED
