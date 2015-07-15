#include"io.h"
#include<string.h>
#include<inttypes.h>
#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>


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

    FILE *cusumlevels;
    if ((cusumlevels = fopen("output/levels.csv","w"))==NULL)
    {
        printf("Cannot open levels file\n");
        abort();
    }

    FILE *blockages;
    if ((blockages = fopen("output/blockages.csv","w"))==NULL)
    {
        printf("Cannot open blockages file\n");
        abort();
    }
    uint64_t lasttime = 0;

    fprintf(events,"Index,\
Type,\
Start Time (s),\
Time Since Last (s),\
Length (us),\
Threshold ,\
Baseline Before (pA),\
Baseline After,\
Effective Baseline,\
Area (pC),\
Average Blockage (pA),\
Relative Average Blockage (pA) ,\
Max Blockage (pA),\
Relative Max Blockage,\
Max Blockage Length (us),\
Num Levels,\
Level Current (pA),\
Level Length (us),\
Blockages (pA) \n");
    while (current)
    {
        if (current->type == 0)
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
                    %d,",\
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
                    current->numlevels);
            lasttime = current->start;
            while (level)
            {
                fprintf(events,"%g",level->current);
                if (level->next)
                {
                    fprintf(events,";");
                }
                fprintf(cusumlevels,"%g\n",level->current);
                fprintf(blockages,"%g\n",level->current-0.5*(current->baseline_after+current->baseline_before));
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
            fprintf(events,"\n");
        }
        current = current->next;
    }
    fclose(events);
    fclose(cusumlevels);
    fclose(blockages);
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
        if (current_event->type == 0)
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
    for (i=0; i<length; i++)
    {
        fprintf(output,"%g,%g,%g\n",i*timestep,current->signal[i], current->filtered_signal[i]);
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
        else if (strcmp(name,"order") == 0)
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
    while (current)
    {
        total++;
        switch (current->type)
        {
            case 0:
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
        }
        current = current->next;
    }
    printf("\n\nError summary:\n%"PRIu64" (%.2f%%) good events and %"PRIu64" (%.2f%%) bad events\n",good, (100.0 * (double) good)/total, bad, (100.0 * (double) bad)/total);
    printf("\n%"PRIu64" (%.2f%%) were discarded for being too short\n",tooshort, (100.0 * (double) tooshort)/total);
    if ((100.0 * (double) tooshort)/total > 5)
    {
        printf("\tYou seem to have a lot of short events\n");
        printf("\tIf you are picking up events that look like noise, try increasing threshold\n");
        printf("\tIf you have lots of short events, try reducing event_minpoints\n\n");
    }
    printf("\n%"PRIu64" (%.2f%%) were discarded for being too long\n",toolong, (100.0 * (double) toolong)/total);
    if ((100.0 * (double) toolong)/total > 5)
    {
        printf("\tYou seem to have a lot of long events\n");
        printf("\tIf your experiment is relatively free of clogs, try increasing event_maxpoints\n\n");
    }
    printf("\n%"PRIu64" (%.2f%%) were discarded for having too few levels\n",badlevels, (100.0 * (double) badlevels)/total);
    if ((100.0 * (double) badlevels)/total > 5)
    {
        printf("\tYou seem to have a lot of events that are too short for CUSUM to process properly with the current settings\n");
        printf("\tTo simply get rid of them, try increasing event_minpoints\n");
        printf("\tTo try to process them properly, try reducing cusum_delta, reducing cusum_minstep, or reducing subevent_minpoints to increase sensitivity\n\n");
    }
    printf("\n%"PRIu64" (%.2f%%) were discarded because the current trace could not be populated\n",badtrace, (100.0 * (double) badtrace)/total);
    if ((100.0 * (double) badtrace)/total > 5)
    {
        printf("\tI have no idea why this number is nonzero\n");
        printf("\tCongratulations, you may have broken logic\n");
    }
    printf("\n%"PRIu64" (%.2f%%) were discarded because the event did return to baseline\n",badbaseline, (100.0 * (double) badbaseline)/total);
    if ((100.0 * (double) badbaseline)/total > 5)
    {
        printf("\tYou seem to have a lot of events that do not have similar baseline before and after\n");
        printf("\tTry increasing hysteresis first, then threshold as well if that doesn't help\n");
    }


    fprintf(logfile,"\n\nError summary:\n%"PRIu64" (%.2f%%) good events and %"PRIu64" (%.2f%%) bad events\n",good, (100.0 * (double) good)/total, bad, (100.0 * (double) bad)/total);
    fprintf(logfile,"\n%"PRIu64" (%.2f%%) were discarded for being too short\n",tooshort, (100.0 * (double) tooshort)/total);
    if ((100.0 * (double) tooshort)/total > 5)
    {
        fprintf(logfile,"\tYou seem to have a lot of short events\n");
        fprintf(logfile,"\tIf you are picking up events that look like noise, try increasing threshold\n");
        fprintf(logfile,"\tIf you have lots of short events, try reducing event_minpoints\n\n");
    }
    fprintf(logfile,"\n%"PRIu64" (%.2f%%) were discarded for being too long\n",toolong, (100.0 * (double) toolong)/total);
    if ((100.0 * (double) toolong)/total > 5)
    {
        fprintf(logfile,"\tYou seem to have a lot of long events\n");
        fprintf(logfile,"\tIf your experiment is relatively free of clogs, try increasing event_maxpoints\n\n");
    }
    fprintf(logfile,"\n%"PRIu64" (%.2f%%) were discarded for having too few levels\n",badlevels, (100.0 * (double) badlevels)/total);
    if ((100.0 * (double) badlevels)/total > 5)
    {
        fprintf(logfile,"\tYou seem to have a lot of events that are too short for CUSUM to process properly with the current settings\n");
        fprintf(logfile,"\tTo simply get rid of them, try increasing event_minpoints\n");
        fprintf(logfile,"\tTo try to process them properly, try reducing cusum_delta, reducing cusum_minstep, or reducing subevent_minpoints to increase sensitivity\n\n");
    }
    fprintf(logfile,"\n%"PRIu64" (%.2f%%) were discarded because the current trace could not be populated\n",badtrace, (100.0 * (double) badtrace)/total);
    if ((100.0 * (double) badtrace)/total > 5)
    {
        fprintf(logfile,"\tI have no idea why this number is nonzero\n");
        fprintf(logfile,"\tCongratulations, you may have broken logic\n");
    }
    fprintf(logfile,"\n%"PRIu64" (%.2f%%) were discarded because the event did return to baseline\n",badbaseline, (100.0 * (double) badbaseline)/total);
    if ((100.0 * (double) badbaseline)/total > 5)
    {
        fprintf(logfile,"\tYou seem to have a lot of events that do not have similar baseline before and after\n");
        fprintf(logfile,"\tTry increasing hysteresis first, then threshold as well if that doesn't help\n");
    }
}

