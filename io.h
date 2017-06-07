/*

                                COPYRIGHT
    Copyright (C) 2015-2017 Kyle Briggs (kbrig035<at>uottawa.ca)

    This file is part of CUSUM.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

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

union int16bits
{
    uint16_t bits;
    int16_t currentval;
};
void print_license(FILE *logfile);
void check_filesize(configuration *config, FILE *input);
void initialize_files(io_struct *io, configuration *config);
void free_io(io_struct *io);
void configure_defaults(configuration *config);
void config_sanity_check(configuration *config, FILE *logfile);
void print_error_summary(FILE *logfile, int64_t *error_summary, int64_t numevents);
void initialize_events_file(FILE *events, FILE *rate, FILE *baselinefile);
void print_event_line(FILE *events, FILE *rate, event *current, double timestep, int64_t );


void output_baseline_stats(FILE *baselinefile, baseline_struct *baseline_stats, int64_t pos, double samplingfreq);
FILE *read_config(configuration *config, const char *version);
void swapByteOrder(double *current, uint64_t *rawsignal, int64_t length);
void swapByteOrder_int16(double *current, uint16_t *rawsignal, int64_t length);
void chimera_gain(double *current, uint16_t *rawsignal, int64_t length, chimera *daqsetup);
int64_t read_current_chimera(FILE *input, double *current, uint16_t *rawsignal, int64_t position, int64_t length, chimera *daqsetup);
int64_t read_current_double(FILE *input, double *current, uint64_t *rawsignal, int64_t position, int64_t length);
int64_t read_current_int16(FILE *input, double *current, uint16_t *rawsignal, int64_t position, int64_t length);
int64_t read_current(FILE *input, double *signal, void *rawsignal, int64_t position, int64_t length, int datatype, chimera *daqsetup);
void print_event_signal(int64_t index, event *current, double timestep, char *eventsfolder);
void print_signal(event *current, int64_t length, char *filename, double timestep);
#endif // IO_H_INCLUDED
