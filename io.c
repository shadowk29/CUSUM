/*

                                COPYRIGHT
    Copyright (C) 2015 Kyle Briggs (kbrig035<at>uottawa.ca)

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
#include"io.h"
void print_license(FILE *logfile)
{
    char license[] =
"Copyright (C) 2015 Kyle Briggs (kbrig035<at>uottawa.ca)\n\nThis file is part of CUSUM.\n\n\
This program is free software: you can redistribute it and/or modify \
it under the terms of the GNU General Public License as published by \
the Free Software Foundation, either version 3 of the License, or \
(at your option) any later version.\n\n\
This program is distributed in the hope that it will be useful, \
but WITHOUT ANY WARRANTY; without even the implied warranty of \
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the \
GNU General Public License for more details.\n\n\
You should have received a copy of the GNU General Public License \
along with this program. If not, see <http://www.gnu.org/licenses/>.";

    fprintf(logfile,"%s\n\n",license);
}
void check_filesize(configuration *config, FILE *input, chimera_file *chimera_input)
{
    int64_t filesize = get_filesize(input, config->datatype, chimera_input);
    if (config->finish == 0)
    {
        config->finish = filesize;
    }
    else
    {
        config->finish = filesize < config->finish ? filesize : config->finish;
    }
}
void initialize_files(io_struct *io, configuration *config)
{
    io->input = fopen64_and_check(config->filepath,"rb", 4);
    io->events = fopen64_and_check(config->eventsfile,"w",21);
    io->rate = fopen64_and_check(config->ratefile,"w",21);
    io->baselinefile = fopen64_and_check(config->baselinefile,"w",21);
    initialize_events_file(io->events, io->rate, io->baselinefile);
}

void free_io(io_struct *io)
{
    fclose(io->logfile);
    fclose(io->input);
    fclose(io->events);
    fclose(io->rate);
    fclose(io->baselinefile);
    free(io);
}

void output_baseline_stats(FILE *baselinefile, baseline_struct *baseline_stats, int64_t pos, double samplingfreq)
{
    fprintf(baselinefile, "%g,%g,%g\n",pos/(double) samplingfreq, baseline_stats->mean, baseline_stats->stdev);
}

void print_error_summary(FILE *logfile, int64_t *error_summary, int64_t numevents)
{
    const char* error_type[] = {
    "cusum success",
    "stepfit success",
    "baseline differs",
    "too long",
    "too short",
    "too few levels",
    "cannot read data",
    "cannot pad event",
    "fitted step too small",
    "stepfit found zero",
    "stepfit degenerate",
    "maxiters reached",
    "stepfit failed (f)",
    "stepfit failed (o)",
    "stepfit failed (p)",
    "stepfit out of memory",
    "stepfit invalid input",
    "stepfit user interrupt",
    "stepfit found NaN"};

    int i;


    printf("\n\n--------------------------------------------------------\nEvent Summary: %"PRId64" events detected\n\n",numevents);
    fprintf(logfile,"\n\n--------------------------------------------------------\nEvent Summary: %"PRId64" events detected\n\n",numevents);


    printf("Success: %.3g %%\nFailed: %.3g %%\t\t%s\n--------------------------------------------------------\n",FRACTION_TO_PERCENTAGE*(error_summary[CUSUM]+error_summary[STEPRESPONSE])/(double) numevents, FRACTION_TO_PERCENTAGE*(numevents - (error_summary[CUSUM]+error_summary[STEPRESPONSE]))/(double) numevents,error_type[CUSUM]);
    fprintf(logfile,"Success: %.3g\nFailed: %.3g\t\t%s\n--------------------------------------------------------\n",FRACTION_TO_PERCENTAGE*(error_summary[CUSUM]+error_summary[STEPRESPONSE])/(double) numevents, FRACTION_TO_PERCENTAGE*(numevents - (error_summary[CUSUM]+error_summary[STEPRESPONSE]))/(double) numevents,error_type[STEPRESPONSE]);

    printf("Event Type\tCount\tPercentage\t\tReason\n\n");
    fprintf(logfile,"Event Type\tCount\t\tPercentage\tReason\n\n");

    for (i=0; i<NUMTYPES; i++)
    {
        printf("%d\t\t%"PRId64"\t%.3g %%\t\t%s\n",i,error_summary[i],FRACTION_TO_PERCENTAGE*error_summary[i]/(double)numevents,error_type[i]);
        fprintf(logfile,"%d\t\t%"PRId64"\t%.3g %%\t\t%s\n",i,error_summary[i],FRACTION_TO_PERCENTAGE*error_summary[i]/(double)numevents,error_type[i]);
        if (i==1)
        {
            printf("--------------------------------------------------------\n");
            fprintf(logfile,"--------------------------------------------------------\n");
        }
    }
    printf("--------------------------------------------------------\n");
    fprintf(logfile,"--------------------------------------------------------\n");
}

int64_t read_current(FILE *input, double *signal, void *rawsignal, int64_t position, int64_t length, int datatype, chimera *daqsetup, double savegain, chimera_file *start)
{
    int64_t read = 0;
    if (datatype == 64)
    {
        read = read_current_double(input, signal, (uint64_t *) rawsignal, position, length);
    }
    else if (datatype == 16)
    {
        read = read_current_int16(input, signal, (uint16_t *) rawsignal, position, length, savegain);
    }
    else if (datatype == 0)
    {
        read = read_current_chimera(input, signal, (uint16_t *) rawsignal, position, length, daqsetup);
    }
    else if (datatype == -1)
    {
        read = read_current_chimera_native(start, signal, (uint16_t *) rawsignal, position, length);
    }
    return read;
}

void chimera_gain(double *current, uint16_t *rawsignal, int64_t length, chimera *daqsetup)
{
    int64_t i;
    double closed_loop_gain = daqsetup->TIAgain*daqsetup->preADCgain;
    uint16_t bitmask = (uint16_t) ((1 << 16)-1) - (uint16_t) ((1 << (16-daqsetup->ADCbits)) - 1);
    uint16_t sample;
    for (i=0; i<length; i++)
    {
        current[i] = 0;
        sample = (uint16_t) (rawsignal[i] & bitmask);
        current[i] = daqsetup->ADCvref - (2.0 * daqsetup->ADCvref) * (double) sample / (double) (1<<16);
        current[i] = -current[i] / closed_loop_gain + daqsetup->currentoffset;
        current[i] *= AMPS_TO_PICOAMPS;
    }
}



int64_t read_current_chimera(FILE *input, double *current, uint16_t *rawsignal, int64_t position, int64_t length, chimera *daqsetup)
{
    int64_t test;
    int64_t read = 0;

    if (fseeko64(input,(off64_t) position*sizeof(uint16_t),SEEK_SET))
    {
        return 0;
    }
    test = fread(rawsignal, sizeof(uint16_t), length, input);
    read = test;
    if (test != length)
    {
        printf("Tried to read %"PRId64" samples at position %"PRId64": \n", length, (off64_t) position*sizeof(uint16_t));
        perror("End of file reached\n");
    }

    chimera_gain(current, rawsignal, read, daqsetup);
    return read;
}


chimera_file *chimera_file_by_index(chimera_file *current, int64_t position)
{
    int64_t i = 0;
    while (current && i + current->length < position)
    {
        i += current->length;
        current = current->next;
    }
    if (current == NULL)
    {
        printf("Attempting to read past end up indexed data\n");
        pause_and_exit(-1);
    }
    current->offset = position - i;
    return current;
}

int64_t read_current_chimera_native(chimera_file *start, double *current, uint16_t *rawsignal, int64_t position, int64_t length)
{
    int64_t test;
    int64_t read = 0;
    chimera_file *current_file = start;
    current_file = chimera_file_by_index(start, position);

    if (fseeko64(current_file->data_file,(off64_t) current_file->offset*sizeof(uint16_t),SEEK_SET))
    {
        return 0;
    }
    test = fread(rawsignal, sizeof(uint16_t), intmin(length - read, current_file->length - current_file->offset), current_file->data_file);
    chimera_gain(current, rawsignal, test, current_file->daqsetup);
    read += test;


    while (read < length)
    {
        current_file = current_file->next;
        current_file->offset = 0;
        test = fread(&rawsignal[read], sizeof(uint16_t), intmin(length - read, current_file->length), current_file->data_file);
        chimera_gain(&current[read], &rawsignal[read], test, current_file->daqsetup);
        read += test;
    }
    return read;
}



void swapByteOrder(double *current, uint64_t *rawsignal, int64_t length)
{
    union doublebits bitval;
    int64_t i;
    for (i=0; i<length; i++)
    {
        bitval.bits = rawsignal[2*i];
        bitval.bits = (bitval.bits >> 56) |
          ((bitval.bits<<40) & 0x00FF000000000000) |
          ((bitval.bits<<24) & 0x0000FF0000000000) |
          ((bitval.bits<<8) & 0x000000FF00000000) |
          ((bitval.bits>>8) & 0x00000000FF000000) |
          ((bitval.bits>>24) & 0x0000000000FF0000) |
          ((bitval.bits>>40) & 0x000000000000FF00) |
          (bitval.bits << 56);
        current[i] = bitval.currentval;
    }
}


int64_t read_current_double(FILE *input, double *current, uint64_t *rawsignal, int64_t position, int64_t length)
{
    int64_t test;

    int64_t read = 0;

    if (fseeko64(input,(off64_t) position*2*sizeof(double),SEEK_SET))
    {
        return 0;
    }
    test = fread(rawsignal, sizeof(uint64_t), 2*length, input);
    read = test/2;
    if (test != 2*length)
    {
        printf("Tried to read %"PRId64" samples at position %"PRId64": \n", length, (off64_t) position*2*sizeof(double));
        perror("End of file reached\n");
    }
    swapByteOrder(current, rawsignal, length);
    return read;
}


void swapByteOrder_int16(double *current, uint16_t *rawsignal, int64_t length, double savegain)
{
    union int16bits bitval;

    int64_t i;
    for (i=0; i<length; i++)
    {
        bitval.bits = rawsignal[2*i];
        bitval.bits = (bitval.bits << 8) | (bitval.bits >> 8);
        current[i] = (double) bitval.currentval;
        current[i] *= savegain;
    }
}



int64_t read_current_int16(FILE *input, double *current, uint16_t *rawsignal, int64_t position, int64_t length, double savegain)
{
    int64_t test;

    int64_t read = 0;

    if (fseeko64(input,(off64_t) position*2*sizeof(uint16_t),SEEK_SET))
    {
        return 0;
    }
    test = fread(rawsignal, sizeof(uint16_t), 2*length, input);
    read = test/2;
    if (test != 2*length)
    {
        printf("Tried to read %"PRId64" samples at position %"PRId64": \n", length, (off64_t) position*2*sizeof(uint16_t));
        perror("End of file reached\n");
    }
    swapByteOrder_int16(current, rawsignal, length, savegain);
    return read;
}

void initialize_events_file(FILE *events, FILE *rate, FILE *baselinefile)
{
    fprintf(events,"id,\
type,\
start_time_s,\
event_delay_s,\
duration_us,\
threshold,\
delta,\
baseline_before_pA,\
baseline_after_pA,\
effective_baseline_pA,\
area_pC,\
average_blockage_pA,\
relative_average_blockage,\
max_blockage_pA,\
relative_max_blockage,\
max_blockage_duration_us,\
n_levels,\
rc_const1_us,\
rc_const2_us,\
residual_pA,\
max_deviation_pA,\
min_blockage_pA,\
relative_min_blockage,\
min_blockage_duration_us,\
intra_crossings,\
level_current_pA,\
level_duration_us,\
blockages_pA,\
stdev_pA,\
intra_crossing_times_us\n");

fprintf(rate,"id,\
type,\
start_time_s,\
end_time_s,\
local_stdev,\
local_baseline,\
intra_crossings,\
intra_crossing_times_us\n");

fprintf(baselinefile, "time_s,\
baseline_pA,\
stdev_pA\n");
}

void print_event_line(FILE *events, FILE *rate, event *current, double timestep, int64_t lasttime)
{
#ifdef DEBUG
    printf("Print Line\n");
    fflush(stdout);
#endif // DEBUG
    edge *intraedge = current->intra_edges;
    fprintf(rate,"%"PRId64",\
            %d,\
            %g,\
            %g,\
            %g,\
            %g,\
            %"PRId64",",
            current->index, \
            current->type, \
            current->start * timestep, \
            current->finish * timestep, \
            current->local_stdev, \
            current->local_baseline, \
            current->intracrossings);
    while (intraedge)
    {
        fprintf(rate,"%g",intraedge->location * timestep * SECONDS_TO_MICROSECONDS);
        if (intraedge->next)
        {
            fprintf(rate,";");
        }
        intraedge = intraedge->next;
    }
    fprintf(rate,"\n");


    if (current->type == CUSUM || current->type == STEPRESPONSE)
    {
        intraedge = current->intra_edges;
        cusumlevel *level = current->first_level;
        fprintf(events,"%"PRId64",\
            %d,\
            %g,\
            %g,\
            %g,\
            %g,\
            %g,\
            %g,\
            %g,\
            %g,\
            %g,\
            %g,\
            %g,\
            %g,\
            %g,\
            %g,\
            %d,\
            %g,\
            %g,\
            %g,\
            %g,\
            %g,\
            %g,\
            %g,\
            %"PRId64",",
            current->index, \
            current->type, \
            current->start * timestep, \
            (current->start - lasttime) * timestep, \
            current->length * timestep * SECONDS_TO_MICROSECONDS, \
            current->threshold, \
            current->delta, \
            current->baseline_before, \
            current->baseline_after, \
            0.5 * (current->baseline_after + current->baseline_before),\
            current->area, \
            current->average_blockage, \
            d_abs(current->average_blockage / (0.5 * (current->baseline_before + current->baseline_after))), \
            current->max_blockage, \
            d_abs(current->max_blockage / (0.5 * (current->baseline_before + current->baseline_after))), \
            current->max_length * timestep * SECONDS_TO_MICROSECONDS, \
            current->numlevels, \
            current->rc1 * timestep * SECONDS_TO_MICROSECONDS, \
            current->rc2 * timestep * SECONDS_TO_MICROSECONDS, \
            current->residual, \
            current->maxdeviation, \
            current->min_blockage, \
            d_abs(current->min_blockage / (0.5 * (current->baseline_before + current->baseline_after))), \
            current->min_length * timestep * SECONDS_TO_MICROSECONDS,\
            current->intracrossings);
        while (level)
        {
            fprintf(events,"%g",level->current);
            if (level->next)
            {
                fprintf(events,";");
            }
            level = level->next;
        }
        fprintf(events,",");
        level = current->first_level;
        while (level)
        {
            fprintf(events,"%g",level->length * timestep * SECONDS_TO_MICROSECONDS);
            if (level->next)
            {
                fprintf(events,";");
            }
            level = level->next;
        }
        fprintf(events,",");
        level = current->first_level;
        while (level)
        {
            fprintf(events,"%g",-1.0*signum(current->baseline_before)*(level->current-0.5*(current->baseline_after+current->baseline_before)));
            if (level->next)
            {
                fprintf(events,";");
            }
            level = level->next;
        }
        fprintf(events,",");
        level = current->first_level;
        while (level)
        {
            fprintf(events,"%g",level->stdev);
            if (level->next)
            {
                fprintf(events,";");
            }
            level = level->next;
        }
        fprintf(events,",");
        while (intraedge)
        {
            fprintf(events,"%g",intraedge->location * timestep * SECONDS_TO_MICROSECONDS);
            if (intraedge->next)
            {
                fprintf(events,";");
            }
            intraedge = intraedge->next;
        }
        fprintf(events,"\n");
        fflush(events);
    }
}


void print_signal(event *current, int64_t length, char *filename, double timestep)
{
    FILE *output; //a test file for plotting output
    if ((output = fopen64(filename,"w"))==NULL)
    {
        printf("Cannot open output file\n");
        pause_and_exit(24);
    }
    int64_t i;
    if (current->type == STEPRESPONSE)
    {
        double stepfit;
        double u1 = current->first_level->length;
        double u2 = u1 + current->first_level->next->length;
        for (i=0; i<length; i++)
        {
            stepfit = current->first_level->current;
            if (i > u1)
            {
                stepfit -= (current->first_level->current-current->first_level->next->current)*(1-exp(-(i-u1)/current->rc1));
            }
            if (i > u2)
            {
                stepfit += (current->first_level->next->next->current-current->first_level->next->current)*(1-exp(-(i-u2)/current->rc2));
            }
            fprintf(output,"%g,%g,%g,%g\n",i*timestep,current->signal[i], current->filtered_signal[i],stepfit);
        }
    }
    else
    {
        for (i=0; i<length; i++)
        {
            fprintf(output,"%g,%g,%g\n",i*timestep,current->signal[i], current->filtered_signal[i]);
        }
    }
    fflush(output);
    fclose(output);
}

void print_event_signal(int64_t index, event *current, double timestep, char *eventsfolder, int print_bad)
{
#ifdef DEBUG
    printf("Print Signal\n");
    fflush(stdout);
#endif // DEBUG
    if(!print_bad)
    {
        if (current->type == CUSUM || current->type == STEPRESPONSE)
        {
            char eventname[1024];
            sprintf(eventname,"%s/event_%08"PRId64".csv",eventsfolder,index);
            print_signal(current, current->length + current->padding_before + current->padding_after, eventname, timestep);
        }
    }
    else
    {
        if (current->signal && current->type != TOOLONG)
        {
            char eventname[1024];
            sprintf(eventname,"%s/event_%08"PRId64".csv",eventsfolder,index);
            print_signal(current, current->length + current->padding_before + current->padding_after, eventname, timestep);
        }

    }
}


void configure_defaults(configuration *config)
{
    //deaults for config file, can be overwritten if provided
    config->start = 0;
    config->finish = 0;
    config->usefilter = 0;
    config->eventfilter = 0;
    config->event_direction = 0;
    config->cusum_min_threshold = 0.1;
    config->cusum_max_threshold = 100.0;
    config->maxiters = 2000;
    config->stepfit_samples = 0;
    config->attempt_recovery = 0;
    config->padding_wait = 0;
}

void config_sanity_check(configuration *config, FILE *logfile)
{
    if (config->datatype != 16 && config->datatype != 64 && config->datatype !=0 && config->datatype != -1)
    {
        printf("datatype currently can only be 0, 16, or 64\n");
        pause_and_exit(43);
    }

    int correctionflag = 0;
    printf("\nVerifying config parameters\n\n");
    fprintf(logfile,"\nVerifying config parameters\n\n");
    if (config->threshold < 1 || config->hysteresis < 1)
    {
        printf("Fractional threshold and hysteresis are no longer supported! Use a multiple of the baseline standard deviation.\n");
        correctionflag = 1;
        pause_and_exit(1);
    }

    if (config->order > 10 && (config->usefilter || config->eventfilter))
    {
        printf("Bessel filters of order >10 are not supported. Correction:\npoles=10\n");
        fprintf(logfile,"Bessel filters of order >10 are not supported. Correction:\npoles=10\n");
        config->order = 10;
        correctionflag = 1;
    }
    else if (config->order < 2 && (config->usefilter || config->eventfilter))
    {
        printf("Bessel filters of order <2 are not supported. Correction:\npoles=2\n");
        fprintf(logfile,"Bessel filters of order <2 are not supported. Correction:\npoles=2\n");
        config->order = 2;
        correctionflag = 1;
    }
    else if (config->order % 2 == 1 && (config->usefilter || config->eventfilter))
    {
        printf("Bessel filters of order >10 are not supported. Correction:\npoles=%"PRId64"\n",config->order + 1);
        fprintf(logfile,"Bessel filters of order >10 are not supported. Correction:\npoles=%"PRId64"\n",config->order + 1);
        config->order += 1;
        correctionflag = 1;
    }
    if (config->datatype==0 && config->samplingfreq != (int64_t) config->daqsetup->samplerate)
    {
        printf("Sampling rate does not match Chimera setup. Correction:\nsamplingfreq=%"PRId64"\n",config->samplingfreq);
        fprintf(logfile,"Sampling rate does not match Chimera setup. Correction:\nsamplingfreq=%"PRId64"\n",config->samplingfreq);
        config->samplingfreq = (int64_t) config->daqsetup->samplerate;
        correctionflag = 1;
    }
    if (config->stepfit_samples > 0 && config->stepfit_samples < config->subevent_minpoints)
    {
        printf("Stepfit samples should be at least as large as subevent_minpoints. Correction:\nstepfit_samples=%"PRId64"\n",config->subevent_minpoints);
        fprintf(logfile,"Stepfit samples should be at least as large as subevent_minpoints. Correction:\nstepfit_samples=%"PRId64"\n",config->subevent_minpoints);
        config->stepfit_samples = config->subevent_minpoints;
        correctionflag = 1;
    }
    if (config->stepfit_samples && !config->attempt_recovery)
    {
        printf("Stepfit_samples is on, but attempt_recovery is off. Correction:\nattempt_recovery=1\n");
        fprintf(logfile,"Stepfit_samples is on, but attempt_recovery is off. Correction:\nattempt_recovery=1\n");
        config->attempt_recovery = 1;
        correctionflag = 1;
    }
    if (config->stepfit_samples <= 0 && config->attempt_recovery)
    {
        printf("Attempt_recovery is on, but step_fit samples is off. Correction:\nstepfit_samples=%"PRId64"\n",config->subevent_minpoints);
        fprintf(logfile,"Attempt_recovery is on, but step_fit samples is off. Correction:\nstepfit_samples=%"PRId64"\n",config->subevent_minpoints);
        config->stepfit_samples = config->subevent_minpoints;
        correctionflag = 1;
    }
    if (config->stepfit_samples <= 0 && !config->attempt_recovery && config->event_minpoints < config->subevent_minpoints)
    {
        printf("Stepfit is off, and event_minpoints is less than subevent_minpoints. Correction:\nevent_minpoints=%"PRId64"\n",config->subevent_minpoints);
        fprintf(logfile,"Stepfit is off, and event_minpoints is less than subevent_minpoints. Correction:\nevent_minpoints=%"PRId64"\n",config->subevent_minpoints);
        config->event_minpoints = config->subevent_minpoints;
        correctionflag = 1;
    }
    if (!correctionflag)
    {
        printf("No corrections\n");
    }
    printf("\nDone config check\n\n");
    fprintf(logfile,"\nDone config check\n\n");
}


FILE * read_config(configuration *config, const char *version)
{
    configure_defaults(config);
    char configline[STRLENGTH];
    char *name;
    char *value;
    long int cutoff = 0;
    FILE *configfile;
    configfile = fopen64_and_check("config.txt","r",1);
    config->outputfolder[0]='\0';



    //initialize some defaults
    config->start = 0;
    config->usefilter = 0;
    config->eventfilter = 0;
    config->savegain = 1;
    config->cusum_elasticity = 0;

    while ((fgets(configline, STRLENGTH, configfile)) != NULL)
    {
        name = strtok(configline,"=");
        value = strtok(NULL,"=\n");
        if (strcmp(name,"readlength") == 0)
        {
            config->readlength = strtoull(value,NULL,10);
        }
        if (strcmp(name,"print_bad_events") == 0)
        {
            config->print_bad = strtoull(value,NULL,10);
        }
        else if (strcmp(name,"start") == 0)
        {
            config->start = strtoull(value,NULL,10);
        }
        else if (strcmp(name,"finish") == 0)
        {
            config->finish = strtoull(value,NULL,10);
        }
        else if (strcmp(name,"cutoff") == 0)
        {
            cutoff = strtol(value,NULL,10);
        }
        else if (strcmp(name,"samplingfreq") == 0)
        {
            config->samplingfreq = strtoull(value,NULL,10);
        }
        else if (strcmp(name,"poles") == 0)
        {
            config->order = strtoull(value,NULL,10);
        }
        else if (strcmp(name,"threshold") == 0)
        {
            config->threshold = strtod(value,NULL);
        }
        else if (strcmp(name,"hysteresis") == 0)
        {
            config->hysteresis = strtod(value,NULL);
        }
        else if (strcmp(name,"padding_wait") == 0)
        {
            config->padding_wait = strtod(value,NULL);
        }
        else if (strcmp(name,"baseline_min") == 0)
        {
            config->baseline_min = strtod(value,NULL);
        }
        else if (strcmp(name,"baseline_max") == 0)
        {
            config->baseline_max = strtod(value,NULL);
        }
        else if (strcmp(name,"cusum_min_threshold") == 0)
        {
            config->cusum_min_threshold = strtod(value,NULL);
        }
        else if (strcmp(name,"cusum_max_threshold") == 0)
        {
            config->cusum_max_threshold = strtod(value,NULL);
        }
        else if (strcmp(name,"cusum_delta") == 0)
        {
            config->cusum_delta = strtod(value,NULL);
        }
        else if (strcmp(name,"cusum_elasticity") == 0)
        {
            config->cusum_elasticity = strtod(value,NULL);
        }
        else if (strcmp(name,"cusum_minstep") == 0)
        {
            config->cusum_minstep = strtod(value,NULL);
        }
        else if (strcmp(name,"intra_threshold") == 0)
        {
            config->intra_threshold = strtod(value,NULL);
        }
        else if (strcmp(name,"intra_hysteresis") == 0)
        {
            config->intra_hysteresis = strtod(value,NULL);
        }
        else if (strcmp(name,"event_minpoints") == 0)
        {
            config->event_minpoints = strtoull(value,NULL,10);
        }
        else if (strcmp(name,"event_maxpoints") == 0)
        {
            config->event_maxpoints = strtoull(value,NULL,10);
        }
        else if (strcmp(name,"subevent_minpoints") == 0)
        {
            config->subevent_minpoints = strtoull(value,NULL,10);
        }
        else if (strcmp(name,"event_direction") == 0)
        {
            config->event_direction = strtol(value,NULL,10);
        }
        else if (strcmp(name,"use_filter") == 0)
        {
            config->usefilter = strtol(value,NULL,10);
        }
        else if (strcmp(name,"event_filter") == 0)
        {
            config->eventfilter = strtol(value,NULL,10);
        }
        else if (strcmp(name,"input_file") == 0)
        {
            strncpy(config->filepath,value,STRLENGTH-1);
            config->filepath[STRLENGTH-1]='\0';
        }
        else if (strcmp(name,"output_folder") == 0)
        {
            strncpy(config->outputfolder,value,STRLENGTH-1);
            config->outputfolder[STRLENGTH-1]='\0';
        }
        else if (strcmp(name,"datatype") == 0)
        {
            config->datatype = strtol(value,NULL,10);
        }
        else if (strcmp(name,"savegain") == 0)
        {
            config->savegain = strtod(value,NULL);
        }
        else if (strcmp(name,"stepfit_samples") == 0)
        {
            config->stepfit_samples = strtol(value,NULL,10);
        }
        else if (strcmp(name,"maxiters") == 0)
        {
            config->maxiters = strtol(value,NULL,10);
        }
        else if (strcmp(name,"attempt_recovery") == 0)
        {
            config->attempt_recovery = strtol(value,NULL,10);
        }
        else if (strcmp(name,"SETUP_ADCSAMPLERATE") == 0)
        {
            config->daqsetup->samplerate = strtod(value,NULL);
        }
        else if (strcmp(name,"SETUP_pAoffset") == 0)
        {
            config->daqsetup->currentoffset = strtod(value,NULL);
        }
        else if (strcmp(name,"SETUP_TIAgain") == 0)
        {
            config->daqsetup->TIAgain = strtod(value,NULL);
        }
        else if (strcmp(name,"SETUP_ADCVREF") == 0)
        {
            config->daqsetup->ADCvref = strtod(value,NULL);
        }
        else if (strcmp(name,"SETUP_ADCBITS") == 0)
        {
            config->daqsetup->ADCbits = strtol(value,NULL,10);
        }
        else if (strcmp(name,"SETUP_preADCgain") == 0)
        {
            config->daqsetup->preADCgain = strtod(value,NULL);
        }
    }

    if (config->usefilter == 0 && config->eventfilter == 0)
    {
        config->order = 0;
        config->cutoff = 0;
    }
    else
    {
        config->cutoff = 2.0 *(double) cutoff/(double) config->samplingfreq;
    }


    int test;
    if (config->outputfolder[0] == '\0' || config->outputfolder[0] == '\n')
    {
        if ((test=snprintf(config->outputfolder,STRLENGTH-1,"output")) < 0 || test >= STRLENGTH)
        {
            printf("Cannot write default output folder string\n");
            pause_and_exit(1);
        }
    }
    if ((test=snprintf(config->eventsfolder,STRLENGTH-1,"%s/events",config->outputfolder)) < 0 || test >= STRLENGTH)
    {
        printf("Cannot write eventsfolder string\n");
        pause_and_exit(1);
    }
    if ((test=snprintf(config->eventsfile,STRLENGTH-1,"%s/events.csv",config->outputfolder)) < 0 || test >= STRLENGTH)
    {
        printf("Cannot write eventsfile string\n");
        pause_and_exit(1);
    }
    if ((test=snprintf(config->ratefile,STRLENGTH-1,"%s/rate.csv",config->outputfolder)) < 0 || test >= STRLENGTH)
    {
        printf("Cannot write ratefile string\n");
        pause_and_exit(1);
    }
    if ((test=snprintf(config->logfile,STRLENGTH-1,"%s/summary.txt",config->outputfolder)) < 0 || test >= STRLENGTH)
    {
        printf("Cannot write logfile string\n");
        pause_and_exit(1);
    }
    if ((test=snprintf(config->baselinefile,STRLENGTH-1,"%s/baseline.csv",config->outputfolder)) < 0 || test >= STRLENGTH)
    {
        printf("Cannot write baselinefile string\n");
        pause_and_exit(1);
    }

    config->chimera_input = NULL;
    if (config->datatype == -1)
    {
        config->chimera_input = index_chimera_files(config->filepath);
    }

    FILE *logfile;
    logfile = fopen64_and_check(config->logfile,"w", 4);
    printf("Using CUSUM version %s\n",version);
    fprintf(logfile,"Using CUSUM version %s\n",version);
    fseeko64(configfile,0,SEEK_SET);
    fprintf(logfile, "<----CONFIGURATION BEGINS---->\n\n");
    while ((fgets(configline, STRLENGTH, configfile)) != NULL)
    {
        fprintf(logfile, "%s", configline);
    }
    fprintf(logfile, "<----CONFIGURATION ENDS---->\n\n");
    config_sanity_check(config, logfile);
    fclose(configfile);
    return logfile;
}

chimera_file *index_chimera_files(char *filepath)
{
    char basedir[STRLENGTH];
    char basename[STRLENGTH];
    char settingsfilename[STRLENGTH];
    char datafilename[STRLENGTH];
    char expname[STRLENGTH];
    char templatename[STRLENGTH];
    char *pch;
    char *prev = NULL;
    int forwardflag = 0;
    pch = strchr(filepath, '\\');
    if (pch == NULL)
    {
        pch = strchr(filepath, '/');
        forwardflag = 1;
    }
    while (pch != NULL)
    {
        prev = pch + 1;
        if (forwardflag == 1)
        {
            pch = strchr(pch+1, '/');
        }
        else
        {
            pch = strchr(pch+1, '\\');
        }
    }
    strncpy(basedir, filepath, (size_t) (prev - filepath));
    basedir[prev-filepath] = '\0';


    if ((pch = strstr(filepath, ".log")) != NULL)
    {
        snprintf(templatename, strlen(filepath) - strlen(basedir) - 19, "%s", &filepath[strlen(basedir)]);
    }
    DIR *directory;
    struct dirent *dir;
    directory = opendir(basedir);
    if (directory == NULL)
    {
        printf("Could not open data directory %s\n", basedir);
        pause_and_exit(-1);
    }


    chimera_file *head = initialize_chimera_files();
    while ((dir = readdir(directory)) != NULL)
    {
        if ((pch = strstr(dir->d_name, ".log")) != NULL)
        {
            strncpy(basename, dir->d_name, (size_t) (pch - dir->d_name));
            basename[pch - dir->d_name] = '\0';
            snprintf(settingsfilename, STRLENGTH, "%s%s.settings", basedir, basename);
            snprintf(datafilename, STRLENGTH, "%s%s.log", basedir, basename);
            snprintf(expname, strlen(basename)-15, "%s", basename);
            if (strcmp(expname, templatename) == 0)
            {
                head = add_chimera_file(head, datafilename, settingsfilename);
            }
        }
    }
    closedir(directory);
    fflush(stdout);
    chimera_file *current;
    current = head;
    double offset = current->daqsetup->currentoffset;
    while (current != NULL)
    {
        if (current->daqsetup->currentoffset != offset)
        {
            printf("Current offset changes during recording, limit dataset to a single offset value\n");
            pause_and_exit(-1);
        }
        current = current->next;
    }
    return head;
}


