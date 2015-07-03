#include"io.h"
#include"detector.h"
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
    double currentlevel;
    double currenttime;
    uint64_t i;
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

    fprintf(events,"Index,\
Type,\
Padding (us),\
Start Time (s),\
Length (us),\
Threshold ,\
Baseline Before (pA),\
Baseline After,\
Effective Baseline,\
Area (pC),\
Average Blockage (pA),\
Relative Average Blockage (pA) ,\
Max Blockage (pA),\
Relative Max Blockage ,\
Num Levels,\
Level Current (pA),\
Level Length (us),\
Level Current (pA),\
Level Length (us),\
Level Current (pA),\
Level Length (us),\
Level Current (pA),\
Level Length (us),\
Level Current (pA),\
Level Length (us),\
Level Current (pA),\
Level Length (us),\
Level Current (pA),\
Level Length (us)\n");
    while (current)
    {
        if (current->type == 0)
        {
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
                    %d,",\
                    current->index, \
                    current->type, \
                    current->padding * timestep * 1e6, \
                    current->start * timestep, \
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
                    current->numlevels);

            currenttime = 0;
            currentlevel = current->filtered_signal[0];
            fprintf(events,"%g,",currentlevel);
            for (i=0; i<current->length + 2*current->padding; i++)
            {
                if (signum(current->filtered_signal[i] - currentlevel) != 0)
                {
                    currentlevel = current->filtered_signal[i];
                    fprintf(events,"%g,%g,",(i-currenttime) * timestep * 1e6, currentlevel);
                    fprintf(cusumlevels,"%g\n",currentlevel);
                    currenttime = i;
                }
            }
            fprintf(events,"%g\n",(i-currenttime) * timestep * 1e6);
        }
        current = current->next;
    }
    fclose(events);
    fclose(cusumlevels);
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
        sprintf(eventname,"output/event_%05d.csv",index);

        print_signal(current, current->length+2*current->padding, eventname, timestep);
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
            config->poles = strtoull(value,NULL,10);
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
        config->poles = 0;
        cutoff = 0;
    }
    config->cutoff = 2.0 *(double) cutoff/(double) config->samplingfreq;
    fclose(configfile);
}
