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
#define _VERSION_ "1.9.0"


int main()
{
    //read the configuration file

    configuration *config;
    config = calloc_and_check(1,sizeof(configuration),1);
    config->daqsetup = calloc_and_check(1,sizeof(chimera),2);

    FILE *logfile;
    logfile = fopen64_and_check("output/summary.txt","w",3);

    printf("#### Using CUSUM Version %s ####\n\n",_VERSION_);
    fprintf(logfile,"#### Using CUSUM Version %s ####\n",_VERSION_);

    read_config(config, logfile);

    if (config->datatype != 16 && config->datatype != 64 && config->datatype !=0)
    {
        printf("datatype currently can only be 0, 16, or 64\n");
        exit(43);
    }

    //open the input file
    FILE *input;
    input = fopen64_and_check(config->filepath,"rb", 4);

    FILE *events;
    events = fopen64_and_check("output/events.csv","w",21);

    FILE *rate;
    rate = fopen64_and_check("output/rate.csv","w",21);

    initialize_events_file(events, rate);

    bessel *lpfilter = NULL;
    double *filtered = NULL;

    uint64_t *error_summary = calloc_and_check(NUMTYPES, sizeof(uint64_t), 5);

    //initialize the low-pass filter and allocate necessary memory
    if (config->usefilter || config->eventfilter)
    {
        lpfilter = initialize_filter(lpfilter, config->order, config->cutoff, config->readlength, config->samplingfreq);
    }
    //allocate memory for file reading
    double *signal = calloc_and_check(config->readlength,sizeof(double), 6);

    histostruct *histogram;
    histogram = calloc_and_check(1,sizeof(histostruct),7);
    histogram->histogram = NULL;
    histogram->numbins = 0;

    //find out how big the file is for use in a progressbar
    uint64_t filesize = get_filesize(input, config->datatype);
    if (config->finish == 0)
    {
        config->finish = filesize;
    }
    else
    {
        config->finish = filesize < config->finish ? filesize : config->finish;
    }


    //initialize linked list to store the locations of edges in the input file
    edge *head_edge;
    edge *current_edge;
    head_edge = initialize_edges();
    current_edge = head_edge;

    //initialize struct to store the information about events found using the edge list
    event *current_event;
    current_event = initialize_events();

    fprintf(logfile, "<----RUN LOG BEGINS---->\n\n");
    printf("Locating events... \n");
    fprintf(logfile, "Locating events...\n ");
    fflush(stdout);

    double risetime = 5; //FIXME
    uint64_t typeswitch = 0;
    double baseline;
    double badbaseline = 0;
    double goodbaseline = 0;
    uint64_t i;
    uint64_t read;
    uint64_t pos;
    int endflag;
    endflag = 0;
    read = 0;
    pos = 0;
    char progressmsg[STRLENGTH];
    for (pos = config->start; pos < config->finish; pos += read)
    {
        snprintf(progressmsg,STRLENGTH*sizeof(char)," %g seconds processed",(pos-config->start)/(double) config->samplingfreq);
        progressbar(pos-config->start,config->finish-config->start,progressmsg);
        read = read_current(input, signal, pos, intmin(config->readlength,config->finish - pos), config->datatype, config->daqsetup);
        if (read < config->readlength || feof(input))
        {
            endflag = 1;
        }
        if (config->usefilter)
        {
            filter_signal(signal, lpfilter, read);
        }
        baseline = build_histogram(signal, histogram, read, config->binsize, config->baseline_max, config->baseline_min);
        if (baseline < config->baseline_min || baseline > config->baseline_max)
        {
            badbaseline += read;
        }
        else
        {
            goodbaseline += read;
            current_edge = detect_edges(signal, baseline, read, current_edge, config->threshold, config->hysteresis, pos, config->event_direction);
        }
        if (endflag)
        {
            pos += read;
                break;
        }
        memset(signal,'0',(config->readlength)*sizeof(double));
    }
    snprintf(progressmsg,STRLENGTH*sizeof(char)," %g seconds processed",(pos-config->start)/(double) config->samplingfreq);
    progressbar(pos-config->start,config->finish-config->start,progressmsg);
    printf("\nRead %g seconds of good baseline\nRead %g seconds of bad baseline\n", goodbaseline/(double) config->samplingfreq, badbaseline / (double) config->samplingfreq);
    fprintf(logfile, "\nRead %g seconds of good baseline\nRead %g seconds of bad baseline\n", goodbaseline/(double) config->samplingfreq, badbaseline / (double) config->samplingfreq);


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
    if (config->usefilter)
    {
        free(filtered);
    }


    uint64_t index = 0;
    uint64_t numevents = 0;
    uint64_t edgecount;
    uint64_t edgenum = 0;
    uint64_t edges;


    edgecount = count_edges(current_edge);
    current_edge = head_edge;


    uint64_t lasttime = config->start;
    uint64_t last_end = config->start;
    printf("Processing %"PRIu64" edges\n", edgecount);
    while (current_edge)
    {
        snprintf(progressmsg,STRLENGTH*sizeof(char)," %"PRIu64" events processed",numevents);
        progressbar(edgenum, edgecount, progressmsg);
        edges = get_next_event(current_event, current_edge, index);
        edgenum += edges;
        for (i=0; i<edges; i++)
        {
            current_edge = current_edge->next;
        }
        index++;
        filter_event_length(current_event, config->event_maxpoints, config->event_minpoints, config->stepfit_samples);
        generate_trace(input, current_event, config->datatype, logfile, lpfilter, config->eventfilter, config->daqsetup, config->samplingfreq, current_edge, last_end, config->start);
        last_end = current_event->finish;
        cusum(current_event, config->cusum_delta, config->cusum_min_threshold, config->cusum_max_threshold, config->subevent_minpoints);
        typeswitch += average_cusum_levels(current_event, config->subevent_minpoints, config->cusum_minstep, config->attempt_recovery);
        step_response(current_event, risetime, config->maxiters, config->cusum_minstep);
        populate_event_levels(current_event);
        calculate_level_noise(current_event, config->subevent_minpoints);
        refine_event_estimates(current_event);
        event_baseline(current_event, config->baseline_min, config->baseline_max);
        event_max_blockage(current_event);
        event_area(current_event, 1.0/config->samplingfreq);
        print_event_signal(current_event->index, current_event, 1.0/config->samplingfreq*1e6);
        print_event_line(events, rate, current_event, 1.0/config->samplingfreq, lasttime);
        lasttime = current_event->start;
        current_edge = current_edge->next;
        edgenum++;
        numevents++;
        error_summary[current_event->type]++;
        free_single_event(current_event);
    }
    snprintf(progressmsg,STRLENGTH*sizeof(char)," %"PRIu64" events processed",numevents);
    progressbar(edgenum, edgecount, progressmsg);

    print_error_summary(logfile, error_summary, numevents);


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
    free(error_summary);


    if (config->usefilter)
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
