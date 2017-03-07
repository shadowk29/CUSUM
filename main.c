/*

                                COPYRIGHT
    Copyright (C) 2015-2016 Kyle Briggs (kbrig035<at>uottawa.ca)

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
#define _VERSION_ "3.1.3"

int main()
{
    check_bits(); //verify that doubles are 64-bit for bit-shifting later
    configuration *config = calloc_and_check(1,sizeof(configuration),"Cannot allocate config struct");
    config->daqsetup = calloc_and_check(1,sizeof(chimera),"Cannot allocate DAQ struct");

    io_struct *io = calloc_and_check(1, sizeof(io_struct),"Cannot allocate io_struct");
    io->logfile = read_config(config, _VERSION_); //read config.txt for the run setup and open the logfile
    initialize_files(io, config); //open all other files used for IO

    bessel *lpfilter = initialize_filter(config->usefilter, config->eventfilter, config->order, config->cutoff, config->readlength); //if turned on, intialize the bessel filter

    int64_t *error_summary = calloc_and_check(NUMTYPES, sizeof(int64_t), "Cannot allocate error array"); //error summary to keep track of why analysis fails

    signal_struct *sig = initialize_signal(config, lpfilter ? lpfilter->padding : 0); //allocate memory for all the main signal arrays used for IO

    baseline_struct *baseline_stats = initialize_baseline(config); //keeps track of local baseline

    check_filesize(config, io->input); //check how big the file is and assign file reading parameters accordingly

    edge *head_edge, *current_edge; //initialize linked list to store the locations of edges in the input file
    head_edge = initialize_edges();
    current_edge = head_edge;

    //main loop over input file, finds and logs rough event locations for a more focused pass later
    current_edge = find_edges(config, io, sig, baseline_stats, lpfilter, current_edge, head_edge);

    free(sig->paddedsignal);

    int64_t index = 0;
    int64_t edgecount;
    int64_t edgenum = 0;
    int64_t edges;


    edgecount = count_edges(current_edge);
    printf("Processing %"PRId64" edges\n", edgecount);
    event *current_event;
    current_event = initialize_events();

    int64_t lasttime = config->start;
    int64_t last_end = config->start;

    time_t start_time;
    time_t curr_time;
    char progressmsg[STRLENGTH];
    int64_t numevents = 0;
    int64_t i;
    int64_t typeswitch = 0;
    time(&start_time);
    while (current_edge)
    {
#ifdef DEBUG
    printf("Main Loop\n");
    fflush(stdout);
#endif // DEBUG
        snprintf(progressmsg,STRLENGTH*sizeof(char)," %"PRId64" events processed",numevents);
        progressbar(edgenum, edgecount, progressmsg,difftime(time(&curr_time),start_time));
        edges = get_next_event(current_event, current_edge, index);
        edgenum += edges;
        for (i=0; i<edges; i++)
        {
            current_edge = current_edge->next;
        }
        index++;
        identify_step_events(current_event, config->stepfit_samples, config->subevent_minpoints, config->attempt_recovery);
        filter_long_events(current_event, config->event_maxpoints);
        filter_short_events(current_event, config->event_minpoints);
        generate_trace(io->input, current_event, config->datatype, sig->rawsignal, io->logfile, lpfilter, config->eventfilter, config->daqsetup, current_edge, last_end, config->start, config->subevent_minpoints);
        last_end = current_event->finish;
        cusum(current_event, config->cusum_delta, config->cusum_min_threshold, config->cusum_max_threshold, config->subevent_minpoints);
        typeswitch += average_cusum_levels(current_event, config->subevent_minpoints, config->cusum_minstep, config->attempt_recovery);
        step_response(current_event, config->usefilter || config->eventfilter ? 2.0/config->cutoff : 5, config->maxiters, config->cusum_minstep);
        populate_event_levels(current_event);
        calculate_level_noise(current_event, config->subevent_minpoints);
        refine_event_estimates(current_event);
        event_baseline(current_event, config->baseline_min, config->baseline_max);
        event_max_blockage(current_event);
        event_area(current_event, 1.0/config->samplingfreq);
        print_event_signal(current_event->index, current_event, 1.0/config->samplingfreq*SECONDS_TO_MICROSECONDS,config->eventsfolder);
        print_event_line(io->events, io->rate, current_event, 1.0/config->samplingfreq, lasttime);
        lasttime = current_event->start;
        current_edge = current_edge->next;
        edgenum++;
        numevents++;
        error_summary[current_event->type]++;
        free_single_event(current_event);
#ifdef DEBUG
    printf("Done\n");
    fflush(stdout);
#endif // DEBUG
    }
    snprintf(progressmsg,STRLENGTH*sizeof(char)," %"PRId64" events processed",numevents);
    progressbar(edgenum, edgecount, progressmsg,difftime(time(&curr_time),start_time));

    print_error_summary(io->logfile, error_summary, numevents);

    printf("\nCleaning up memory usage...\n");
    fprintf(io->logfile, "Cleaning up memory usage...\n");
    free(current_event);
    free_edges(head_edge);
    free(config->daqsetup);

    free(error_summary);
    free(sig->rawsignal);
    free(sig);
    free_baseline(baseline_stats);


    if (config->usefilter || config->eventfilter)
    {
        free_filter(lpfilter);
    }
    free(config);

    fprintf(io->logfile, "<----RUN LOG ENDS---->\n\n");
    free_io(io);
    return 0;
}
