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
#include"utils.h"
#include"detector.h"
#include"bessel.h"
#include"stepfit.h"


int main()
{
    //read the configuration file
    configuration *config;
    if ((config = malloc(sizeof(configuration)))==NULL)
    {
        printf("Cannot allocate config structure\n");
        exit(1);
    }
    if ((config->daqsetup = malloc(sizeof(chimera)))==NULL)
    {
        printf("Cannot allocate chimera setup structure\n");
        exit(2);
    }

    FILE *logfile;
    if ((logfile = fopen64("output/summary.txt","w"))==NULL)
    {
        printf("Cannot open summary file\n");
        exit(3);
    }
    read_config(config, logfile);


    //open the input file
    FILE *input;
    if ((input = fopen64(config->filepath,"rb"))==NULL)
    {
        printf("Cannot open input file\n");
        exit(4);
    }

    FILE *events;
    if ((events = fopen("output/events.csv","w"))==NULL)
    {
        printf("Cannot open event summary file\n");
        exit(21);
    }


    FILE *rate;
    if ((rate = fopen("output/rate.csv","w"))==NULL)
    {
        printf("Cannot open event rate file\n");
        exit(21);
    }

    initialize_events_file(events, rate);

    //file reading variables
    uint64_t i;
    uint64_t read;
    uint64_t pos;
    int endflag;

    endflag = 0;
    read = 0;
    pos = 0;


    //event location identification parameters
    double threshold;
    double hysteresis;
    double baseline_min;
    double baseline_max;
    int event_direction;
    double binsize;

    threshold = config->threshold;
    hysteresis = config->hysteresis;
    event_direction = config->event_direction;
    baseline_min = config->baseline_min;
    baseline_max = config->baseline_max;
    binsize = config->binsize;

    //low-pass filter parameters
    uint64_t order;
    int usefilter;
    int eventfilter;
    double cutoff;
    uint64_t samplingfreq;

    order = config->order;
    cutoff = config->cutoff;
    usefilter = config->usefilter;
    eventfilter = config->eventfilter;
    samplingfreq = config->samplingfreq;

    bessel *lpfilter = NULL;
    double *filtered = NULL;

    //CUSUM parameters
    int attempt_recovery;

    attempt_recovery = config->attempt_recovery;
    double risetime = 5; //FIXME
    //event requirement paramaters
    uint64_t typeswitch = 0;

    //initialize the low-pass filter and allocate necessary memory
    if (usefilter || eventfilter)
    {
        lpfilter = initialize_filter(lpfilter, order, cutoff, config->readlength, config->samplingfreq);

        if ((filtered = (double *) calloc(config->readlength,sizeof(double)))==NULL)
        {
            printf("Cannot allocate filtered signal array\n");
            exit(5);
        }
    }
    //allocate memory for file reading
    double *signal = calloc_and_check(config->readlength,sizeof(double));
    /*if ((signal = (double *) calloc(config->readlength,sizeof(double)))==NULL)
    {
        printf("Cannot allocate signal array\n");
        exit(6);
    }*/

    if (config->datatype != 16 && config->datatype != 64 && config->datatype !=0)
    {
        printf("datatype currently can only be 0, 16, or 64\n");
        exit(43);
    }



    histostruct *histogram;
    if ((histogram = malloc(sizeof(histostruct)))==NULL)
    {
        printf("cannot allocate histogram structure\n");
        exit(7);
    }
    histogram->histogram = NULL;
    histogram->numbins = 0;

    //find out how big the file is for use in a progressbar
    uint64_t filesize = get_filesize(input, config->datatype);
    config->finish = filesize < config->finish ? filesize : config->finish;


    //dummy variables
    double baseline;
    double badbaseline = 0;
    double goodbaseline = 0;

    //initialize linked list to store the locations of edges in the input file
    edge *head_edge;
    edge *current_edge;
    head_edge = initialize_edges();
    current_edge = head_edge;

    //initialize linked list to store the information about events found using the edge list
    event *current_event;
    current_event = initialize_events();




    fprintf(logfile, "<----RUN LOG BEGINS---->\n\n");
    printf("Locating events... \n");
    fprintf(logfile, "Locating events...\n ");
    fflush(stdout);
    for (pos = config->start; pos < config->finish; pos += read)
    {
        progressbar(pos-config->start,config->finish-config->start);
        if (config->datatype == 64)
        {
            read = read_current(input, signal, pos, intmin(config->readlength,config->finish - pos));
        }
        else if (config->datatype == 16)
        {
            read = read_current_int16(input, signal, pos, intmin(config->readlength,config->finish - pos));
        }
        else if (config->datatype == 0)
        {
            read = read_current_chimera(input, signal, pos, intmin(config->readlength,config->finish - pos), config->daqsetup);
        }

        if (read < config->readlength || feof(input))
        {
            endflag = 1;
        }




        if (usefilter)
        {
            filter_signal(signal, filtered, lpfilter, read);

            /*FILE *filtertest;
            if ((filtertest=fopen64("filtertest.csv","w"))==NULL)
            {
                printf("Cannot open output file for filter data\n");
                exit(1);
            }
            for (i=0; i<read; i++)
            {
                fprintf(filtertest,"%g,%g,%g\n",i/(double)samplingfreq,signal[i],filtered[i]);
            }
            exit(1);*/
            //baseline = baseline_averaging(filtered, read, baseline_min, baseline_max);
            baseline = build_histogram(filtered, histogram, read, binsize, baseline_max, baseline_min);
            if (baseline < baseline_min || baseline > baseline_max)
            {
                badbaseline += read;
            }
            else
            {
                goodbaseline += read;
                current_edge = detect_edges(filtered, baseline, read, current_edge, threshold, hysteresis, pos, event_direction);
            }

            if (endflag)
            {
                pos += read;
                break;
            }
            memset(filtered,'0',(config->readlength)*sizeof(double));
        }
        else
        {
            //baseline = baseline_averaging(signal, read, baseline_min, baseline_max);
            baseline = build_histogram(signal, histogram, read, binsize, baseline_max, baseline_min);
            if (baseline < baseline_min || baseline > baseline_max)
            {
                badbaseline += read;
            }
            else
            {
                goodbaseline += read;
                current_edge = detect_edges(signal, baseline, read, current_edge, threshold, hysteresis, pos, event_direction);
            }

            if (endflag)
            {
                pos += read;
                break;
            }
        }
        memset(signal,'0',(config->readlength)*sizeof(double));
    }
    progressbar(pos-config->start,config->finish-config->start);
    printf("\nRead %g seconds of good baseline\nRead %g seconds of bad baseline\n", goodbaseline/(double) samplingfreq, badbaseline / (double) samplingfreq);
    fprintf(logfile, "\nRead %g seconds of good baseline\nRead %g seconds of bad baseline\n", goodbaseline/(double) samplingfreq, badbaseline / (double) samplingfreq);


    current_edge = head_edge;

    if (!current_edge || current_edge->type == HEAD)
    {
        printf("It is %"PRId64"\n",current_edge->type);
        printf("No edges found in signal, exiting\n");
        fprintf(logfile, "No edges found in signal, exiting\n");
        system("pause");
        exit(8);
    }


    free(signal);
    if (usefilter)
    {
        free(filtered);
    }


    uint64_t index = 0;
    uint64_t edgecount;
    uint64_t edgenum = 0;
    uint64_t edges;

    edgecount = count_edges(current_edge);
    current_edge = head_edge;


    uint64_t lasttime = config->start;
    printf("Processing %"PRIu64" edges\n", edgecount);
    while (current_edge)
    {
        progressbar(edgenum, edgecount);
        edges = get_next_event(current_event, current_edge, index);
        edgenum += edges;
        for (i=0; i<edges; i++)
        {
            current_edge = current_edge->next;
        }
        index++;
        filter_event_length(current_event, config->event_maxpoints, config->event_minpoints, config->stepfit_samples);
        generate_trace(input, current_event, config->datatype, logfile, lpfilter, eventfilter, config->daqsetup, samplingfreq);
        cusum(current_event, config->cusum_delta, config->cusum_min_threshold, config->cusum_max_threshold, config->subevent_minpoints);
        typeswitch += average_cusum_levels(current_event, config->subevent_minpoints, config->cusum_minstep, attempt_recovery);
        step_response(current_event, risetime, config->maxiters, config->cusum_minstep);
        populate_event_levels(current_event);
        calculate_level_noise(current_event, config->subevent_minpoints);
        refine_event_estimates(current_event);
        event_baseline(current_event, baseline_min, baseline_max);
        event_max_blockage(current_event);
        event_area(current_event, 1.0/samplingfreq);
        print_event_signal(current_event->index, current_event, 1.0/samplingfreq*1e6);
        print_event_line(events, rate, current_event, 1.0/samplingfreq, lasttime);
        lasttime = current_event->start;
        current_edge = current_edge->next;
        edgenum++;
        free_single_event(current_event);
    }
    progressbar(edgenum, edgecount);


    printf("\nCleaning up memory usage...\n");
    fprintf(logfile, "Cleaning up memory usage...\n");
    free(current_event);


    free_edges(head_edge);
    for (i=0; i<histogram->numbins; i++)
    {
        free(histogram->histogram[i]);
    }
    free(histogram->histogram);
    free(histogram);
    fclose(input);
    free(config->daqsetup);
    free(config);


    if (usefilter)
    {
        free_filter(lpfilter);
    }


    fprintf(logfile, "<----RUN LOG ENDS---->\n\n");
    fclose(logfile);
    fclose(events);
    fclose(rate);
    system("pause");
    return 0;
}
