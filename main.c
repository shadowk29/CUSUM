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
#include"butterworth.h"
#include"stepfit.h"


int main()
{
    //read the configuration file
    configuration *config;
    if ((config = malloc(sizeof(configuration)))==NULL)
    {
        printf("Cannot allocate config structure\n");
        abort();
    }
    if ((config->daqsetup = malloc(sizeof(chimera)))==NULL)
    {
        printf("Cannot allocate chimera setup structure\n");
        abort();
    }

    FILE *logfile;
    if ((logfile = fopen64("output/summary.txt","w"))==NULL)
    {
        printf("Cannot open summary file\n");
        abort();
    }
    read_config(config, logfile);


    //open the input file
    FILE *input;
    if ((input = fopen64(config->filepath,"rb"))==NULL)
    {
        printf("Cannot open input file\n");
        abort();
    }

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
    double cutoff;
    uint64_t samplingfreq;

    order = config->order;
    cutoff = config->cutoff;
    usefilter = config->usefilter;
    samplingfreq = config->samplingfreq;

    butterworth *lpfilter = NULL;
    double *filtered = NULL;

    //CUSUM parameters
    double cusum_delta;
    double cusum_minstep;
    double cusum_min_threshold;
    double cusum_max_threshold;
    uint64_t subevent_minpoints;
    int refine_estimates;
    uint64_t stepfit_samples;
    uint64_t maxiters;
    int attempt_recovery;

    cusum_delta = config->cusum_delta;
    cusum_minstep = config->cusum_minstep;
    cusum_min_threshold = config->cusum_min_threshold;
    cusum_max_threshold = config->cusum_max_threshold;
    subevent_minpoints = config->subevent_minpoints;
    refine_estimates = config->refine_estimates;
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
    if (usefilter)
    {
        lpfilter = initialize_filter(lpfilter, order, cutoff, readlength);

        if ((filtered = (double *) calloc(readlength,sizeof(double)))==NULL)
        {
            printf("Cannot allocate filtered signal array\n");
            abort();
        }
    }
    //allocate memory for file reading
    double *signal;
    if ((signal = (double *) calloc(readlength,sizeof(double)))==NULL)
    {
        printf("Cannot allocate signal array\n");
        abort();
    }

    if (datatype != 16 && datatype != 64 && datatype !=0)
    {
        printf("datatype currently can only be 0, 16, or 64\n");
        abort();
    }



    histostruct *histogram;
    if ((histogram = malloc(sizeof(histostruct)))==NULL)
    {
        printf("cannot allocate histogram structure\n");
        abort();
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
    event *head_event;
    event *current_event;
    head_event = initialize_events();
    current_event = head_event;

    fflush(stdout);



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
    printf("Finished\n\n");
    fprintf(logfile, "\nRead %g seconds of good baseline\nRead %g seconds of bad baseline\n", goodbaseline/(double) samplingfreq, badbaseline / (double) samplingfreq);
    fprintf(logfile, "Finished\n\n");

    current_edge = head_edge;
    current_event = head_event;

    if (!current_edge || current_edge->type == HEAD)
    {
        printf("No edges found in signal, exiting\n");
        fprintf(logfile, "No edges found in signal, exiting\n");
        system("pause");
        return -1;
    }

    printf("Processing event locations... ");
    fprintf(logfile, "Processing event locations... ");
    current_event = process_edges(current_edge, current_event);
    current_event = head_event;
    printf("Finished\n\n");
    fprintf(logfile, "Finished\n\n");


    if (!current_event || current_event->index == HEAD)
    {
        printf("No events found in signal, exiting\n");
        fprintf(logfile, "No events found in signal, exiting\n");
        system("pause");
        return -2;
    }

    printf("Filtering on event length... ");
    fprintf(logfile, "Filtering on event length... ");
    filter_event_length(current_event, maxpoints, minpoints, logfile, stepfit_samples); //divide by type CUSUM and type STEPRESPONSE based on length
    current_event = head_event;
    printf("Finished\n\n");
    fprintf(logfile, "Finished\n\n");
    fflush(logfile);

    printf("Populating event traces... ");
    fprintf(logfile, "Populating event traces... ");
    populate_event_traces(input, current_event, datatype, logfile, lpfilter, config->eventfilter, config->daqsetup, samplingfreq); //both types
    current_event = head_event;
    printf("Finished\n\n");
    fprintf(logfile, "Finished\n\n");
    fflush(logfile);


    printf("Detecting subevents...\n");
    fprintf(logfile, "Detecting subevents...");
    detect_subevents(current_event, cusum_delta, cusum_min_threshold, cusum_max_threshold, subevent_minpoints); //only CUSUM types
    current_event = head_event;
    printf("\nFinished\n\n");
    fprintf(logfile, "\nFinished\n\n");
    fflush(logfile);

    printf("Processing subevents...");
    fprintf(logfile, "Processing subevents...");
    typeswitch = assign_cusum_levels(current_event, subevent_minpoints, cusum_minstep, attempt_recovery); //only CUSUM types - if < 3 levels found (modify), assign STEPRESPONSE TYPE
    current_event = head_event;
    printf("Finished\n\n");
    fprintf(logfile, "Finished\n\n");
    printf("CUSUM failed to find a blockage in %"PRIu64" events, which will be processed with StepResponse fitting\n",typeswitch);
    fprintf(logfile, "CUSUM failed to find a blockage in %"PRIu64" events, which will be processed with StepResponse fitting\n",typeswitch);
    fflush(logfile);

    //stepresponse function called - make sure no cusum overlap remains
    printf("Fitting short events...\n");
    fprintf(logfile, "Fitting short events...");
    step_response_events(current_event, risetime, maxiters, cusum_minstep);
    current_event = head_event;
    printf("\nFinished\n\n");
    fprintf(logfile, "\nFinished\n\n");
    fflush(logfile);


    printf("Assigning subevents...");
    fprintf(logfile, "Assigning subevents...");
    populate_all_levels(current_event); //both types
    current_event = head_event;
    printf("Finished\n\n");
    fprintf(logfile, "Finished\n\n");
    fflush(logfile);

    printf("Assigning noise...");
    fprintf(logfile, "Assigning noise...");
    calculate_event_noise(current_event, subevent_minpoints);
    current_event = head_event;
    printf("Finished\n\n");
    fprintf(logfile, "Finished\n\n");
    fflush(logfile);



    if (refine_estimates)
    {
        printf("Refining subevent estimates...");
        fprintf(logfile, "Refining subevent estimates...");
        refine_all_estimates(current_event);
        current_event = head_event;
        printf("Finished\n\n");
        fprintf(logfile, "Finished\n\n");
        fflush(logfile);


        printf("Filtering on refined event length... ");
        fprintf(logfile, "Filtering on refined event length... ");
        filter_event_length(current_event, maxpoints, minpoints, logfile, 0);
        current_event = head_event;
        printf("Finished\n\n");
        fprintf(logfile, "Finished\n\n");
        fflush(logfile);
    }

    printf("Assigning event baselines...");
    fprintf(logfile, "Assigning event baselines...");
    assign_event_baselines(current_event, logfile, baseline_min, baseline_max);
    current_event = head_event;
    printf("Finished\n\n");
    fprintf(logfile, "Finished\n\n");
    fflush(logfile);


    printf("Assigning max blockage...");
    fprintf(logfile, "Assigning max blockage...");
    find_max_blockages(current_event);
    current_event = head_event;
    printf("Finished\n\n");
    fprintf(logfile, "Finished\n\n");
    fflush(logfile);

    printf("Assigning event areas...");
    fprintf(logfile, "Assigning event areas...");
    assign_event_areas(current_event, 1.0/samplingfreq);
    current_event = head_event;
    printf("Finished\n\n");
    fprintf(logfile, "Finished\n\n");
    fflush(logfile);


    printf("Printing all signals...\n");
    fprintf(logfile, "Printing all signals...");
    print_all_signals(current_event, 1.0/samplingfreq*1e6);
    current_event = head_event;
    printf("Finished\n\n");
    fprintf(logfile, "Finished\n\n");
    fflush(logfile);


    printf("Printing event summary...");
    fprintf(logfile, "Printing event summary...");
    print_events(current_event, 1.0/samplingfreq);
    current_event = head_event;
    printf("Finished\n\n");
    fprintf(logfile, "Finished\n\n");
    fflush(logfile);

    print_error_summary(current_event, logfile);
    current_event = head_event;

    printf("Cleaning up memory usage...\n");
    fprintf(logfile, "Cleaning up memory usage...\n");
    free_events(head_event);
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
    printf("Finished\n\n");
    fprintf(logfile, "Finished\n\n");
    fprintf(logfile, "<----RUN LOG ENDS---->\n\n");
    fclose(logfile);
    system("pause");
    return 0;
}
