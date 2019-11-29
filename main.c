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
#define _VERSION_ "3.3.0b"


int main()
{
    print_license(stdout);
    check_bits(); //verify that doubles are 64-bit for bit-shifting later
    configuration *config = calloc_and_check(1,sizeof(configuration),"Cannot allocate config struct");
    config->daqsetup = calloc_and_check(1,sizeof(chimera),"Cannot allocate DAQ struct");

    io_struct *io = calloc_and_check(1, sizeof(io_struct),"Cannot allocate io_struct");
    io->logfile = read_config(config, _VERSION_); //read config.txt for the run setup and open the logfile
    print_license(io->logfile);
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

    //no longer needed, and a substantial memory hog
    free(sig->paddedsignal);

    //count the number of edges to be processed
    int64_t edgecount;
    edgecount = count_edges(current_edge);
    printf("Processing %"PRId64" edges\n", edgecount);
    current_edge = head_edge;

    //initialize events
    event *current_event;
    current_event = initialize_events();

    //estimate time statistics for variable cusum parameters
    duration_struct *current_duration = initialize_durations();
    duration_struct *head_duration = current_duration;
    timestruct *timestats = calloc_and_check(1, sizeof(timestruct),"Cannot allocate timestruct");
    estimate_time_statistics(current_duration, timestats, current_edge);
    free_durations(head_duration);

    //main loop over edges, actual fitting and event output happens here
    int64_t numevents = fit_events(config, io, sig->rawsignal, current_event, lpfilter, current_edge, error_summary, edgecount, timestats);

    print_error_summary(io->logfile, error_summary, numevents);

    printf("\nCleaning up memory usage...\n");
    fprintf(io->logfile, "Cleaning up memory usage...\n");
    free(current_event);
    free_edges(head_edge);
    free(config->daqsetup);
    free(timestats);
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
    printf("\nFinished successfully\nPress any key to continue\n");
    getchar();
    return 0;
}
