#ifndef IO_H_INCLUDED
#define IO_H_INCLUDED
#include"utils.h"
#include"detector.h"
#include<string.h>
#include<inttypes.h>
#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>
void configure_defaults(configuration *config);
void config_sanity_check(configuration *config, FILE *logfile);
void print_error_summary(FILE *logfile, uint64_t *error_summary, uint64_t numevents);
void initialize_events_file(FILE *events, FILE *rate);
void print_event_line(FILE *events, FILE *rate, event *current, double timestep, uint64_t );


//void print_events(event *current, double timestep);
void print_histogram(char *filename, histostruct *histogram);
void read_config(configuration *config, FILE *logfile);
inline void swapByteOrder(uint64_t *ull);
inline void swapByteOrder_int16(uint16_t *ull);
double chimera_gain(uint64_t sample, chimera *daqsetup);
uint64_t read_current_chimera(FILE *input, double *current, uint64_t position, uint64_t length, chimera *daqsetup);
uint64_t read_current_double(FILE *input, double *current, uint64_t position, uint64_t length);
uint64_t read_current_int16(FILE *input, double *current, uint64_t position, uint64_t length);
uint64_t read_current(FILE *input, double *signal, uint64_t position, uint64_t length, int datatype, chimera *daqsetup);
void print_event_signal(int index, event *current, double timestep);
void print_signal(event *current, int length, char *filename, double timestep);
//void print_error_summary(event *current, FILE *logfile);
#endif // IO_H_INCLUDED
