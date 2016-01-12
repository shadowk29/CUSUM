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


void print_events(event *current, double timestep)
{
    FILE *events;
    if ((events = fopen("output/events.csv","w"))==NULL)
    {
        printf("Cannot open event summary file\n");
        abort();
    }
    FILE *rejected;
    if ((rejected = fopen("output/rejected.csv","w"))==NULL)
    {
        printf("Cannot open event rejected events file\n");
        abort();
    }
    FILE *rate;
    if ((rate = fopen("output/rate.csv","w"))==NULL)
    {
        printf("Cannot open event rate events file\n");
        abort();
    }
    uint64_t lasttime = 0;
    uint64_t lasttime_rate = 0;
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

fprintf(rejected,"id,\
type,\
start_time_s\n");
    while (current)
    {
        fprintf(rate,"%"PRId64",\
                %d,\
                %.6f,\
                %.6f,\
                %.6f\n", \
                current->index, \
                current->type, \
                current->start * timestep,\
                current->finish * timestep,\
                (current->start - lasttime_rate) * timestep);
        lasttime_rate = current->start;
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
            lasttime = current->start;
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
        }
        else
        {
            fprintf(rejected,"%"PRId64",\
                    %d,\
                    %.6f\n", \
                    current->index, \
                    current->type, \
                    current->start * timestep);
        }
        current = current->next;
    }
    fclose(events);
    fclose(rejected);
    fclose(rate);
}



void print_all_signals(event *current_event, double timestep)
{
    event *head_event = current_event;
    uint64_t numevents = 0;
    while (current_event)
    {
        numevents++;
        current_event = current_event->next;
    }
    current_event = head_event;
    while (current_event != NULL)
    {
        progressbar(current_event->index, numevents);
        if (current_event->type == CUSUM || current_event->type == STEPRESPONSE)
        {
            print_event_signal(current_event->index, current_event, timestep);
        }
        current_event = current_event->next;
    }
}

void print_signal(event *current, int length, char *filename, double timestep)
{
    FILE *output; //a test file for plotting output
    if ((output = fopen(filename,"w"))==NULL)
    {
        printf("Cannot open output file\n");
        abort();
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
    fclose(output);
}

void print_event_signal(int index, event *current, double timestep)
{
    while (current && current->index != index)
    {
        current = current->next;
    }
    if (!current || current->index != index)
    {
        printf("Cannot locate event with index %d\n",index);
        return;
    }
    if (current->index != -1 && current->type != -1)
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
        abort();
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


uint64_t read_current(FILE *input, double *current, uint64_t position, uint64_t length)
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

void export_trace(double *signal, uint64_t length, char *file, double timestep, double start_time)
{
    FILE *output;
    if ((output = fopen64(file,"w"))==NULL)
    {
        printf("Cannot open input file\n");
        abort();
    }
    uint64_t i;

    for (i=0; i<length; i++)
    {
        fprintf(output,"%lf,%lf\n",start_time + i * timestep, signal[i]);
    }
    fclose(output);
}


void read_config(configuration *config, FILE *logfile)
{
    char configline[STRLENGTH];
    char *name;
    char *value;
    long int cutoff = 0;
    FILE *configfile;
    if ((configfile = fopen("config.txt","r"))==NULL)
    {
        printf("Cannot find config file: \"config.txt\"!");
        abort();
    }
    fprintf(logfile, "<----CONFIGURATION BEGINS---->\n\n");
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
        else if (strcmp(name,"input_file") == 0)
        {
            strncpy(config->filepath,value,STRLENGTH-1);
            config->filepath[STRLENGTH-1]='\0';
        }
        else if (strcmp(name,"tracefile") == 0)
        {
            strncpy(config->tracefile,value,STRLENGTH-1);
            config->tracefile[STRLENGTH-1]='\0';
        }
        else if (strcmp(name,"export_trace_start") == 0)
        {
            config->export_trace_start = strtoull(value,NULL,10);
        }
        else if (strcmp(name,"export_trace_end") == 0)
        {
            config->export_trace_end = strtoull(value,NULL,10);
        }
        else if (strcmp(name,"datatype") == 0)
        {
            config->datatype = strtol(value,NULL,10);
        }
        else if (strcmp(name,"refine_estimates") == 0)
        {
            config->refine_estimates = strtol(value,NULL,10);
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

    }
    fprintf(logfile, "<----CONFIGURATION ENDS---->\n\n");
    if (config->usefilter == 0)
    {
        config->order = 0;
        cutoff = 0;
    }
    config->cutoff = 2.0 *(double) cutoff/(double) config->samplingfreq;
    fclose(configfile);
}


void print_error_summary(event *current, FILE *logfile)
{
    uint64_t good = 0;
    uint64_t total = 0;
    uint64_t bad = 0;
    uint64_t badbaseline = 0;
    uint64_t tooshort = 0;
    uint64_t toolong = 0;
    uint64_t badlevels = 0;
    uint64_t badtrace = 0;
    uint64_t badfit = 0;
    uint64_t fitrange = 0;
    uint64_t fitsign = 0;
    uint64_t fitstep = 0;
    uint64_t fitdir = 0;
    while (current)
    {
        total++;
        switch (current->type)
        {
            case CUSUM:
                good++;
                break;
            case STEPRESPONSE:
                good++;
                break;
            case BADBASELINE:
                badbaseline++;
                bad++;
                break;
            case TOOSHORT:
                tooshort++;
                bad++;
                break;
            case TOOLONG:
                toolong++;
                bad++;
                break;
            case BADTRACE:
                badtrace++;
                bad++;
                break;
            case BADLEVELS:
                badlevels++;
                bad++;
                break;
            case BADFIT:
                badfit++;
                bad++;
                break;
            case FITRANGE:
                fitrange++;
                bad++;
                break;
            case FITSIGN:
                fitsign++;
                bad++;
                break;
            case FITSTEP:
                fitstep++;
                bad++;
                break;
            case FITDIR:
                fitdir++;
                bad++;
                break;
        }
        current = current->next;
    }
    printf("\n\nError summary:\n%"PRIu64" (%.2f%%) good events and %"PRIu64" (%.2f%%) bad events\n",good, (100.0 * (double) good)/total, bad, (100.0 * (double) bad)/total);
    printf("\n%"PRIu64" (%.2f%%) too short (type %d)\n",tooshort, (100.0 * (double) tooshort)/total,TOOSHORT);
    if ((100.0 * (double) tooshort)/total > 5)
    {
        printf("\tYou seem to have a lot of short events\n");
        printf("\tIf you are picking up events that look like noise, try increasing threshold\n");
        printf("\tIf you have lots of short events, try reducing event_minpoints\n\n");
    }
    printf("\n%"PRIu64" (%.2f%%) too long (type %d)\n",toolong, (100.0 * (double) toolong)/total,TOOLONG);
    if ((100.0 * (double) toolong)/total > 5)
    {
        printf("\tYou seem to have a lot of long events\n");
        printf("\tIf your experiment is relatively free of clogs, try increasing event_maxpoints\n\n");
    }
    printf("\n%"PRIu64" (%.2f%%) too few levels (type %d)\n",badlevels, (100.0 * (double) badlevels)/total,BADLEVELS);
    if ((100.0 * (double) badlevels)/total > 5)
    {
        printf("\tYou seem to have a lot of events that are too short for CUSUM to process properly with the current settings\n");
        printf("\tTo simply get rid of them, try increasing event_minpoints\n");
        printf("\tTo try to process them properly, try reducing cusum_delta, reducing cusum_minstep, or reducing subevent_minpoints to increase sensitivity\n\n");
    }
    printf("\n%"PRIu64" (%.2f%%) current trace could not be populated (type %d)\n",badtrace, (100.0 * (double) badtrace)/total,BADTRACE);
    if ((100.0 * (double) badtrace)/total > 5)
    {
        printf("\tI have no idea why this number is nonzero\n");
        printf("\tCongratulations, you may have broken logic\n");
    }
    printf("\n%"PRIu64" (%.2f%%) event did return to baseline (type %d)\n",badbaseline, (100.0 * (double) badbaseline)/total,BADBASELINE);
    if ((100.0 * (double) badbaseline)/total > 5)
    {
        printf("\tYou seem to have a lot of events that do not have similar baseline before and after\n");
        printf("\tTry increasing hysteresis first, then threshold as well if that doesn't help\n");
    }
    printf("\n%"PRIu64" (%.2f%%) stepresponse failed to fit a single-level event:\n",badfit+fitrange+fitstep+fitdir+fitsign, (100.0 * (double) (badfit+fitrange+fitstep+fitdir+fitsign))/total);
    printf("\t%"PRIu64" (%.2f%%) fitted jump locations were out of range (type %d)\n",fitrange, (100.0 * (double) fitrange)/total,FITRANGE);
    printf("\t%"PRIu64" (%.2f%%) sign of the fitted step did not match (type %d)\n",fitsign, (100.0 * (double) fitsign)/total,FITSIGN);
    printf("\t%"PRIu64" (%.2f%%) fitted step was smaller than cusum_minstep (type %d)\n",fitstep, (100.0 * (double) fitstep)/total,FITSTEP);
    printf("\t%"PRIu64" (%.2f%%) fitted step was in the wrong direction (type %d)\n",fitdir, (100.0 * (double) fitdir)/total,FITDIR);
    printf("\t%"PRIu64" (%.2f%%) fit failed to converge (type %d)\n",badfit, (100.0 * (double) badfit)/total,BADFIT);
    if ((100.0 * (double) badfit)/total > 5)
    {
        printf("\tYou seem to have a lot of events that cannot be fit to a single-level event\n");
        printf("\tTry reducing stepfit_samples so that you are not fitting multi-level events with it\n");
    }


    fprintf(logfile,"\n\nError summary:\n%"PRIu64" (%.2f%%) good events and %"PRIu64" (%.2f%%) bad events\n",good, (100.0 * (double) good)/total, bad, (100.0 * (double) bad)/total);
    fprintf(logfile,"\n%"PRIu64" (%.2f%%) too short (type %d)\n",tooshort, (100.0 * (double) tooshort)/total,TOOSHORT);
    if ((100.0 * (double) tooshort)/total > 5)
    {
        fprintf(logfile,"\tYou seem to have a lot of short events\n");
        fprintf(logfile,"\tIf you are picking up events that look like noise, try increasing threshold\n");
        fprintf(logfile,"\tIf you have lots of short events, try reducing event_minpoints\n\n");
    }
    fprintf(logfile,"\n%"PRIu64" (%.2f%%) too long (type %d)\n",toolong, (100.0 * (double) toolong)/total,TOOLONG);
    if ((100.0 * (double) toolong)/total > 5)
    {
        fprintf(logfile,"\tYou seem to have a lot of long events\n");
        fprintf(logfile,"\tIf your experiment is relatively free of clogs, try increasing event_maxpoints\n\n");
    }
    fprintf(logfile,"\n%"PRIu64" (%.2f%%) too few levels (type %d)\n",badlevels, (100.0 * (double) badlevels)/total,BADLEVELS);
    if ((100.0 * (double) badlevels)/total > 5)
    {
        fprintf(logfile,"\tYou seem to have a lot of events that are too short for CUSUM to process properly with the current settings\n");
        fprintf(logfile,"\tTo simply get rid of them, try increasing event_minpoints\n");
        fprintf(logfile,"\tTo try to process them properly, try reducing cusum_delta, reducing cusum_minstep, or reducing subevent_minpoints to increase sensitivity\n\n");
    }
    fprintf(logfile,"\n%"PRIu64" (%.2f%%) wcurrent trace could not be populated (type %d)\n",badtrace, (100.0 * (double) badtrace)/total,BADTRACE);
    if ((100.0 * (double) badtrace)/total > 5)
    {
        fprintf(logfile,"\tI have no idea why this number is nonzero\n");
        fprintf(logfile,"\tCongratulations, you may have broken logic\n");
    }
    fprintf(logfile,"\n%"PRIu64" (%.2f%%) event did return to baseline (type %d)\n",badbaseline, (100.0 * (double) badbaseline)/total,BADBASELINE);
    if ((100.0 * (double) badbaseline)/total > 5)
    {
        fprintf(logfile,"\tYou seem to have a lot of events that do not have similar baseline before and after\n");
        fprintf(logfile,"\tTry increasing hysteresis first, then threshold as well if that doesn't help\n");
    }
    fprintf(logfile,"\n%"PRIu64" (%.2f%%) stepresponse failed to fit a single-level event:\n",badfit+fitrange+fitstep+fitdir+fitsign, (100.0 * (double) (badfit+fitrange+fitstep+fitdir+fitsign))/total);
    fprintf(logfile,"\t%"PRIu64" (%.2f%%) fitted jump locations were out of range (type %d)\n",fitrange, (100.0 * (double) fitrange)/total,FITRANGE);
    fprintf(logfile,"\t%"PRIu64" (%.2f%%) sign of the fitted step did not match (type %d)\n",fitsign, (100.0 * (double) fitsign)/total,FITSIGN);
    fprintf(logfile,"\t%"PRIu64" (%.2f%%) fitted step was smaller than cusum_minstep (type %d)\n",fitstep, (100.0 * (double) fitstep)/total,FITSTEP);
    fprintf(logfile,"\t%"PRIu64" (%.2f%%) fitted step was in the wrong direction (type %d)\n",fitdir, (100.0 * (double) fitdir)/total,FITDIR);
    fprintf(logfile,"\t%"PRIu64" (%.2f%%) fit failed to converge (type %d)\n",badfit, (100.0 * (double) badfit)/total,BADFIT);
    if ((100.0 * (double) badfit)/total > 5)
    {
        fprintf(logfile,"\tYou seem to have a lot of events that cannot be fit to a single-level event\n");
        fprintf(logfile,"\tTry reducing stepfit_samples so that you are not fitting multi-level events with it\n");
    }
}

