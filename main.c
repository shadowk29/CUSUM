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
    initialize_events_file(events);

    //file reading variables
    uint64_t i;
    uint64_t read;
    uint64_t pos;
    uint64_t start;
    uint64_t finish;
    uint64_t readlength;
    int datatype;
    int endflag;

    endflag = 0;
    read = 0;
    pos = 0;
    start = config->start;
    finish = config->finish;
    readlength = config->readlength;
    datatype = config->datatype;


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
    double cusum_minstep;
    uint64_t subevent_minpoints;
    uint64_t stepfit_samples;
    uint64_t maxiters;
    int attempt_recovery;

    cusum_minstep = config->cusum_minstep;
    subevent_minpoints = config->subevent_minpoints;
    stepfit_samples = config->stepfit_samples;
    maxiters = config->maxiters;
    attempt_recovery = config->attempt_recovery;
    double risetime = 5; //FIXME
    //event requirement paramaters
    uint64_t maxpoints;
    uint64_t minpoints;
    maxpoints = config->event_maxpoints;
    minpoints = config->event_minpoints;
    uint64_t typeswitch = 0;

    //initialize the low-pass filter and allocate necessary memory
    if (usefilter || eventfilter)
    {
        lpfilter = initialize_filter(lpfilter, order, cutoff, readlength);

        if ((filtered = (double *) calloc(readlength,sizeof(double)))==NULL)
        {
            printf("Cannot allocate filtered signal array\n");
            exit(5);
        }
    }
    //allocate memory for file reading
    double *signal;
    if ((signal = (double *) calloc(readlength,sizeof(double)))==NULL)
    {
        printf("Cannot allocate signal array\n");
        exit(6);
    }

    if (datatype != 16 && datatype != 64 && datatype !=0)
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
    uint64_t filesize = get_filesize(input, datatype);
    finish = filesize < finish ? filesize : finish;


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
    for (pos = start; pos < finish; pos += read)
    {
        progressbar(pos-start,finish-start);
        if (datatype == 64)
        {
            read = read_current(input, signal, pos, intmin(readlength,finish - pos));
        }
        else if (datatype == 16)
        {
            read = read_current_int16(input, signal, pos, intmin(readlength,finish - pos));
        }
        else if (datatype == 0)
        {
            read = read_current_chimera(input, signal, pos, intmin(readlength,finish - pos), config->daqsetup);
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
            memset(filtered,'0',(readlength)*sizeof(double));
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
        memset(signal,'0',(readlength)*sizeof(double));
    }
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



    start = 0;
    finish = 0;
    uint64_t index = 0;
    uint64_t edgecount;
    uint64_t edgenum = 0;

    edgecount = count_edges(current_edge);
    current_edge = head_edge;


    uint64_t lasttime = 0;
    uint64_t lasttime_rate = 0;
    printf("Processing %"PRIu64" edges\n", edgecount);
    while (current_edge)
    {
        progressbar(edgenum, edgecount);
        while (current_edge->type != 0 && current_edge->next) //if for some reason there are multiple of the same type in a row, skip them.
        {
            current_edge = current_edge->next;
            edgenum++;
        }
        start = current_edge->location;
        while (current_edge->type != 1 && current_edge->next) //if for some reason there are multiple of the same type in a row, skip them.
        {
            current_edge = current_edge->next;
            edgenum++;
        }
        finish = current_edge->location;
        if (finish > start)
        {
            current_event = add_event(current_event, start, finish, index);
            index++;
        }


        //filter events on length, and end processing here if they are too short
        //printf("Filtering\n");
        filter_event_length(current_event, maxpoints, minpoints, stepfit_samples);
        //printf("Finished\n");
        //fflush(stdout);

        if (current_event->type != CUSUM && current_event->type != STEPRESPONSE) //verify that we should continue processing
        {
            //printf("Freeing\n");
            current_edge = current_edge->next;
            free_single_event(current_event);
            //printf("Finished\n");
            lasttime_rate = current_event->start;
            continue;
        }

        //read the data trace from the appropriate file position
        //printf("Reading\n");
        generate_trace(input, current_event, datatype, logfile, lpfilter, eventfilter, config->daqsetup, samplingfreq);
        //printf("Finished\n");
        //fflush(stdout);

        if (current_event->type != CUSUM && current_event->type != STEPRESPONSE) //verify that we should continue processing
        {
            current_edge = current_edge->next;
            free_single_event(current_event);
            lasttime_rate = current_event->start;
            continue;
        }

        if (current_event->type == CUSUM)
        {
            //printf("cusum\n");
            cusum(current_event, config->cusum_delta, config->cusum_min_threshold, config->cusum_max_threshold, subevent_minpoints);
            //printf("Finished\n");
            //fflush(stdout);
            //printf("averaging\n");
            typeswitch += average_cusum_levels(current_event, subevent_minpoints, cusum_minstep, attempt_recovery);
            //printf("Finished\n");
            //fflush(stdout);
        }


        if (current_event->type != CUSUM && current_event->type != STEPRESPONSE) //verify that we should continue processing
        {
            current_edge = current_edge->next;
            free_single_event(current_event);
            lasttime_rate = current_event->start;
            continue;
        }


        if (current_event->type == STEPRESPONSE)
        {
            int status;
            //printf("stepfit\n");
            status = step_response(current_event, risetime, maxiters, cusum_minstep);
            //printf("Finished\n");
            //fflush(stdout);
            if (status != GSL_SUCCESS)
            {
                if (current_event->type == STEPRESPONSE)
                {
                    current_event->type = BADFIT;
                }
            }
        }

        if (current_event->type != CUSUM && current_event->type != STEPRESPONSE) //verify that we should continue processing
        {
            current_edge = current_edge->next;
            free_single_event(current_event);
            lasttime_rate = current_event->start;
            continue;
        }

        //printf("Levels\n");
        populate_event_levels(current_event);
        //printf("Finished\n");
        //fflush(stdout);


        if (current_event->type != CUSUM && current_event->type != STEPRESPONSE) //verify that we should continue processing
        {
            current_edge = current_edge->next;
            free_single_event(current_event);
            lasttime_rate = current_event->start;
            continue;
        }

        //printf("Noise\n");
        calculate_level_noise(current_event, subevent_minpoints);
        //printf("Finished\n");
        //fflush(stdout);


        if (current_event->type != CUSUM && current_event->type != STEPRESPONSE) //verify that we should continue processing
        {
            current_edge = current_edge->next;
            free_single_event(current_event);
            lasttime_rate = current_event->start;
            continue;
        }

        //printf("Refining\n");
        refine_event_estimates(current_event);
        //printf("Finished\n");
        //fflush(stdout);

        if (current_event->type != CUSUM && current_event->type != STEPRESPONSE) //verify that we should continue processing
        {
            current_edge = current_edge->next;
            free_single_event(current_event);
            lasttime_rate = current_event->start;
            continue;
        }

        //printf("Baseline\n");
        event_baseline(current_event, baseline_min, baseline_max);
        //printf("Finished\n");
        //fflush(stdout);

        if (current_event->type != CUSUM && current_event->type != STEPRESPONSE) //verify that we should continue processing
        {
            current_edge = current_edge->next;
            free_single_event(current_event);
            lasttime_rate = current_event->start;
            continue;
        }

        //printf("Maximum\n");
        event_max_blockage(current_event);
        //printf("Finished\n");
        //fflush(stdout);


        if (current_event->type != CUSUM && current_event->type != STEPRESPONSE) //verify that we should continue processing
        {
            current_edge = current_edge->next;
            free_single_event(current_event);
            lasttime_rate = current_event->start;
            continue;
        }

        event_area(current_event, 1.0/samplingfreq);
        //printf("Finished\n");
        //fflush(stdout);


        if (current_event->type != CUSUM && current_event->type != STEPRESPONSE) //verify that we should continue processing
        {
            current_edge = current_edge->next;
            free_single_event(current_event);
            lasttime_rate = current_event->start;
            continue;
        }


        //printf("Printign Trace\n");
        print_event_signal(current_event->index, current_event, 1.0/samplingfreq*1e6);
        //printf("Finished\n");
        //fflush(stdout);
        //printf("Printing summary\n");
        print_event_line(events, current_event, 1.0/samplingfreq, lasttime);
        //printf("Finished\n");
        //fflush(stdout);
        lasttime = current_event->start;
        lasttime_rate = current_event->start;


        current_edge = current_edge->next;
        edgenum++;
        free_single_event(current_event);
    }



    //print_error_summary(current_event, logfile);
    //current_event = head_event;

    printf("Cleaning up memory usage...\n");
    fprintf(logfile, "Cleaning up memory usage...\n");
    free_single_event(current_event);
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
    free(signal);
    if (usefilter)
    {
        free_filter(lpfilter);
        free(filtered);
    }
    fprintf(logfile, "<----RUN LOG ENDS---->\n\n");
    fclose(logfile);
    fclose(events);
    system("pause");
    return 0;
}
