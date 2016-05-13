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

void print_error_summary(FILE *logfile, uint64_t *error_summary, uint64_t numevents)
{
    int i;
    printf("\n\nEvent Summary: %"PRIu64" events detected\n\nEvent Type\tCount\tPercentage\n\n",numevents);
    fprintf(logfile,"\n\nEvent Summary: %"PRIu64" events detected\n\nEvent Type\tCount\tPercentage\n\n",numevents);
    for (i=0; i<NUMTYPES; i++)
    {
        printf("%d\t\t%"PRIu64"\t%.3g %%\n",i,error_summary[i],100.0*error_summary[i]/(double)numevents);
        fprintf(logfile,"%d\t\t%"PRIu64"\t%.3g %%\n",i,error_summary[i],100.0*error_summary[i]/(double)numevents);
    }
}

uint64_t read_current(FILE *input, double *signal, uint64_t position, uint64_t length, int datatype, chimera *daqsetup)
{
    uint64_t read = 0;
    if (datatype == 64)
    {
        read = read_current_double(input, signal, position, length);
    }
    else if (datatype == 16)
    {
        read = read_current_int16(input, signal, position, length);
    }
    else if (datatype == 0)
    {
        read = read_current_chimera(input, signal, position, length, daqsetup);
    }
    return read;
}

double chimera_gain(uint64_t sample, chimera *daqsetup)
{
    double current = 0;
    double closed_loop_gain = daqsetup->TIAgain*daqsetup->preADCgain;
    uint16_t bitmask = (uint16_t) (1 << 16) - (uint16_t) ((1 << (16-daqsetup->ADCbits)) - 1);
    sample = (uint16_t) (sample & bitmask);
    current = daqsetup->ADCvref - (2.0 * daqsetup->ADCvref) * (double) sample / (double) (1<<16);
    current = -current / closed_loop_gain + daqsetup->currentoffset;
    return current * 1e12;
}

uint64_t read_current_chimera(FILE *input, double *current, uint64_t position, uint64_t length, chimera *daqsetup)
{
    uint64_t test;
    uint16_t iv;

    uint64_t i;
    uint64_t read = 0;

    if (fseeko64(input,(off64_t) position*sizeof(uint16_t),SEEK_SET))
    {
        return 0;
    }


    for (i = 0; i < length; i++)
    {
        test = fread(&iv, sizeof(uint16_t), 1, input);
        if (test == 1)
        {
            read++;
            current[i] = chimera_gain(iv, daqsetup);
        }
        else
        {
            perror("End of file reached");
            break;
        }
    }
    return read;
}




inline void swapByteOrder_int16(uint16_t *ull)
{
    *ull = (*ull>>8)|(*ull<<8);
}


uint64_t read_current_int16(FILE *input, double *current, uint64_t position, uint64_t length)
{
    uint64_t test;
    int16_t iv[2];

    uint64_t i;
    uint64_t read = 0;

    if (fseeko64(input,(off64_t) position*2*sizeof(int16_t),SEEK_SET))
    {
        return 0;
    }


    for (i = 0; i < length; i++)
    {
        test = fread(iv, sizeof(int16_t), 2, input);
        if (test == 2)
        {
            read++;
            swapByteOrder_int16((uint16_t *) &iv[0]);
            current[i] = (double) iv[0];
        }
        else
        {
            perror("End of file reached");
            break;
        }

    }
    return read;
}

void initialize_events_file(FILE *events, FILE *rate)
{
    fprintf(events,"id,\
type,\
start_time_s,\
event_delay_s,\
duration_us,\
threshold,\
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
level_current_pA,\
level_duration_us,\
blockages_pA,\
stdev_pA\n");

fprintf(rate,"id,\
type,\
start_time_s\n");
}

void print_event_line(FILE *events, FILE *rate, event *current, double timestep, uint64_t lasttime)
{
#ifdef DEBUG
    printf("Print Line\n");
    fflush(stdout);
#endif // DEBUG
    fprintf(rate,"%"PRId64",\
            %d,\
            %.6f\n",\
            current->index, \
            current->type, \
            current->start * timestep);

    if (current->type == CUSUM || current->type == STEPRESPONSE)
    {
        cusumlevel *level = current->first_level;
        fprintf(events,"%"PRId64",\
            %d,\
            %.6f,\
            %.6f,\
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
            %g,",\
            current->index, \
            current->type, \
            current->start * timestep, \
            (current->start - lasttime) * timestep, \
            current->length * timestep * 1e6, \
            current->threshold, \
            current->baseline_before, \
            current->baseline_after, \
            0.5 * (current->baseline_after + current->baseline_before),\
            current->area, \
            current->average_blockage, \
            d_abs(current->average_blockage / (0.5 * (current->baseline_before + current->baseline_after))), \
            current->max_blockage, \
            d_abs(current->max_blockage / (0.5 * (current->baseline_before + current->baseline_after))), \
            current->max_length * timestep * 1e6, \
            current->numlevels, \
            current->rc1 * timestep * 1e6, \
            current->rc2 * timestep * 1e6, \
            current->residual);
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
            fprintf(events,"%g",level->length * timestep * 1e6);
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
            fprintf(events,"%g",level->current-0.5*(current->baseline_after+current->baseline_before));
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
        fprintf(events,"\n");
        fflush(events);
    }
}


void print_signal(event *current, int length, char *filename, double timestep)
{
    FILE *output; //a test file for plotting output
    if ((output = fopen(filename,"w"))==NULL)
    {
        printf("Cannot open output file\n");
        exit(24);
    }
    int i;
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

void print_event_signal(int index, event *current, double timestep)
{
#ifdef DEBUG
    printf("Print Signal\n");
    fflush(stdout);
#endif // DEBUG
    if (current->type == CUSUM || current->type == STEPRESPONSE)
    {
        char eventname[1024];
        sprintf(eventname,"output/events/event_%05d.csv",index);
        print_signal(current, current->length + current->padding_before + current->padding_after, eventname, timestep);
    }
}



void print_histogram(char *filename, histostruct *histogram)
{
    uint64_t i;
    FILE *output;
    if ((output = fopen(filename, "w"))==NULL)
    {
        printf("Cannot open histogram file\n");
        exit(25);
    }
    for (i=0; i<histogram->numbins; i++)
    {
        fprintf(output,"%g,%g,%g,%g\n",(i+0.5)*histogram->delta+histogram->offset, histogram->histogram[i][0], histogram->histogram[i][1], histogram->histogram[i][2]);
    }
    fclose(output);
}


inline void swapByteOrder(uint64_t *ull)
{
    *ull = (*ull >> 56) |
          ((*ull<<40) & 0x00FF000000000000) |
          ((*ull<<24) & 0x0000FF0000000000) |
          ((*ull<<8) & 0x000000FF00000000) |
          ((*ull>>8) & 0x00000000FF000000) |
          ((*ull>>24) & 0x0000000000FF0000) |
          ((*ull>>40) & 0x000000000000FF00) |
          (*ull << 56);
}


uint64_t read_current_double(FILE *input, double *current, uint64_t position, uint64_t length)
{
    uint64_t test;
    double iv[2];

    uint64_t i;
    uint64_t read = 0;

    if (fseeko64(input,(off64_t) position*2*sizeof(double),SEEK_SET))
    {
        return 0;
    }


    for (i = 0; i < length; i++)
    {
        test = fread(iv, sizeof(double), 2, input);
        if (test == 2)
        {
            read++;
            swapByteOrder((uint64_t *) &iv[0]);
            current[i] = iv[0];
        }
        else
        {
            perror("End of file reached: ");
            break;
        }

    }
    return read;
}


void configure_defaults(configuration *config)
{
    //deaults for config file, can be overwritten if provided
    config->start = 0;
    config->finish = 0;
    config->usefilter = 0;
    config->eventfilter = 0;
    config->binsize = 10;
    config->event_direction = 0;
    config->cusum_min_threshold = 0.1;
    config->cusum_max_threshold = 100.0;
    config->maxiters = 5000;
    config->stepfit_samples = 0;
    config->attempt_recovery = 0;
}

void config_sanity_check(configuration *config, FILE *logfile)
{
    printf("Verifying config parameters\nAny notes below will modify the config file above for the actual run\n\n");
    fprintf(logfile,"Verifying config parameters\nAny notes below will modify the config file above for the actual run\n\n");
    if (config->readlength < 2 * config->event_maxpoints)
    {
        printf("Readlength should be at least 2 times event_maxpoints. Correction:\nreadlength=%"PRIu64"\n",2 * config->event_maxpoints);
        fprintf(logfile,"Readlength should be at least 2 times event_maxpoints. Correction:\nreadlength=%"PRIu64"\n",2 * config->event_maxpoints);
        config->readlength = 2 * config->event_maxpoints;
    }

    if (config->order > 10)
    {
        printf("Bessel filters of order >10 are not supported. Correction:\npoles=10\n");
        fprintf(logfile,"Bessel filters of order >10 are not supported. Correction:\npoles=10\n");
        config->order = 10;
    }
    else if (config->order < 2)
    {
        printf("Bessel filters of order >10 are not supported. Correction:\npoles=2\n");
        fprintf(logfile,"Bessel filters of order >10 are not supported. Correction:\npoles=2\n");
        config->order = 2;
    }
    else if (config->order % 2 == 1)
    {
        printf("Bessel filters of order >10 are not supported. Correction:\npoles=%"PRIu64"\n",config->order + 1);
        fprintf(logfile,"Bessel filters of order >10 are not supported. Correction:\npoles=%"PRIu64"\n",config->order + 1);
        config->order += 1;
    }
    if (config->datatype==0 && config->samplingfreq != (uint64_t) config->daqsetup->samplerate)
    {
        printf("Sampling rate does not match Chimera setup. Correction:\nsamplingfreq=%"PRIu64"\n",config->samplingfreq);
        fprintf(logfile,"Sampling rate does not match Chimera setup. Correction:\nsamplingfreq=%"PRIu64"\n",config->samplingfreq);
        config->samplingfreq = (uint64_t) config->daqsetup->samplerate;
    }
    if (config->stepfit_samples > 0 && config->stepfit_samples < config->subevent_minpoints)
    {
        printf("Stepfit samples should be at least as large as subevent_minpoints. Correction:\nstepfit_samples=%"PRIu64"\n",config->subevent_minpoints);
        fprintf(logfile,"Stepfit samples should be at least as large as subevent_minpoints. Correction:\nstepfit_samples=%"PRIu64"\n",config->subevent_minpoints);
        config->stepfit_samples = config->subevent_minpoints;
    }
    if (config->stepfit_samples && !config->attempt_recovery)
    {
        printf("Stepfit_samples is on, but attempt_recovery is off. Correction:\nattempt_recovery=1\n");
        fprintf(logfile,"Stepfit_samples is on, but attempt_recovery is off. Correction:\nattempt_recovery=1\n");
        config->attempt_recovery = 1;
    }
    if (!config->stepfit_samples && config->attempt_recovery)
    {
        printf("Attempt_recovery is on, but step_fit samples is off. Correction:\nstepfit_samples=%"PRIu64"\n",config->subevent_minpoints);
        fprintf(logfile,"Attempt_recovery is on, but step_fit samples is off. Correction:\nstepfit_samples=%"PRIu64"\n",config->subevent_minpoints);
        config->stepfit_samples = config->subevent_minpoints;
    }
    if (!config->stepfit_samples && !config->attempt_recovery && config->event_minpoints < config->subevent_minpoints)
    {
        printf("Stepfit is off, and event_minpoints is less than subevent_minpoints. Correction:\nevent_minpoints=%"PRIu64"\n",config->subevent_minpoints);
        fprintf(logfile,"Stepfit is off, and event_minpoints is less than subevent_minpoints. Correction:\nevent_minpoints=%"PRIu64"\n",config->subevent_minpoints);
        config->event_minpoints = config->subevent_minpoints;
    }
    if (config->usefilter || config->eventfilter)
    {
        if (config->subevent_minpoints < 8.0/config->cutoff)
        {
            printf("Warning: subevent_minpoints is less than 4RC, levels might be underestimated. Suggest increaseing subevent_minpoints to %"PRIu64"\nNo Correction\n",(uint64_t) (8.0/config->cutoff));
            fprintf(logfile,"Warning: subevent_minpoints is less than 4RC, levels might be underestimated. Suggest increaseing subevent_minpoints to %"PRIu64"\nNo Correction\n",(uint64_t) (8.0/config->cutoff));
        }
        if (config->event_minpoints < 8.0/config->cutoff && config->stepfit_samples == 0)
        {
            printf("Warning: event_minpoints is less than 4RC, short events will not be fit accurately. Suggest increaseing event_minpoints to %"PRIu64"\nNo Correction\n",(uint64_t) (8.0/config->cutoff));
            fprintf(logfile,"Warning: event_minpoints is less than 4RC, short events will not be fit accurately. Suggest increaseing event_minpoints to %"PRIu64"\nNo Correction\n",(uint64_t) (8.0/config->cutoff));
        }
    }
    printf("\nDone config check\n\n");
    fprintf(logfile,"\nDone config check\n\n");
}


void read_config(configuration *config, FILE *logfile)
{
    configure_defaults(config);
    char configline[STRLENGTH];
    char *name;
    char *value;
    long int cutoff = 0;
    FILE *configfile;
    if ((configfile = fopen("config.txt","r"))==NULL)
    {
        printf("Cannot find config file: \"config.txt\"!");
        exit(27);
    }
    fprintf(logfile, "<----CONFIGURATION BEGINS---->\n\n");


    //initialize some defaults
    config->start = 0;
    config->usefilter = 0;
    config->eventfilter = 0;

    while ((fgets(configline, STRLENGTH, configfile)) != NULL)
    {
        fprintf(logfile, "%s", configline);
        name = strtok(configline,"=");
        value = strtok(NULL,"=\n");
        if (strcmp(name,"readlength") == 0)
        {
            config->readlength = strtoull(value,NULL,10);
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
        else if (strcmp(name,"baseline_min") == 0)
        {
            config->baseline_min = strtod(value,NULL);
        }else if (strcmp(name,"baseline_max") == 0)
        {
            config->baseline_max = strtod(value,NULL);
        }
        else if (strcmp(name,"binsize") == 0)
        {
            config->binsize = strtod(value,NULL);
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
        else if (strcmp(name,"cusum_minstep") == 0)
        {
            config->cusum_minstep = strtod(value,NULL);
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
        else if (strcmp(name,"datatype") == 0)
        {
            config->datatype = strtol(value,NULL,10);
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
    fprintf(logfile, "<----CONFIGURATION ENDS---->\n\n");
    if (config->usefilter == 0 && config->eventfilter == 0)
    {
        config->order = 0;
        config->cutoff = 0;
    }
    else
    {
        config->cutoff = 2.0 *(double) cutoff/(double) config->samplingfreq;
    }
    config_sanity_check(config, logfile);
    fclose(configfile);
}

