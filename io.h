#ifndef IO_H_INCLUDED
#define IO_H_INCLUDED
#include"utils.h"
#include<stdio.h>
#include<stdlib.h>

void print_events(event *current, double timestep);
void print_all_signals(event *current_event, double timestep);
void print_histogram(char *filename, histostruct *histogram);
void read_config(configuration *config);
void export_trace(double *signal, uint64_t length, char *file, double timestep, double start_time);
inline void swapByteOrder(uint64_t *ull);
inline void swapByteOrder_int16(uint16_t *ull);
uint64_t read_current(FILE *input, double *current, uint64_t position, uint64_t length);
uint64_t read_current_int16(FILE *input, double *current, uint64_t position, uint64_t length);
void print_event_signal(int index, event *current, double timestep);
void print_signal(event *current, int length, char *filename, double timestep);

#endif // IO_H_INCLUDED
