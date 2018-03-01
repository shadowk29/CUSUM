/*

                                COPYRIGHT
    Copyright (C) 2015-2017 Kyle Briggs (kbrig035<at>uottawa.ca)

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
#include"detector.h"

int64_t fit_events(configuration *config, io_struct *io, double *rawsignal, event *current_event, bessel *lpfilter, edge *current_edge, int64_t *error_summary, int64_t edgecount)
{
    int64_t lasttime = config->start;
    int64_t last_end = config->start;
    int64_t edgenum = 0;
    int64_t edges;
    int64_t index = 0;
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
        generate_trace(io->input, current_event, config->datatype, rawsignal, io->logfile, lpfilter, config->eventfilter, config->daqsetup, current_edge, last_end, config->start, config->subevent_minpoints, config->savegain, config->padding_wait);
        last_end = current_event->finish;
        cusum(current_event, config->cusum_delta, config->cusum_min_threshold, config->cusum_max_threshold, config->subevent_minpoints);
        typeswitch += average_cusum_levels(current_event, config->subevent_minpoints, config->cusum_minstep, config->attempt_recovery, config->padding_wait);
        step_response(current_event, config->usefilter || config->eventfilter ? 2.0/config->cutoff : 5, config->maxiters, config->cusum_minstep);
        populate_event_levels(current_event);
        calculate_level_noise(current_event, config->subevent_minpoints);
        refine_event_estimates(current_event);
        event_baseline(current_event, config->baseline_min, config->baseline_max);
        event_max_blockage(current_event, config->cusum_minstep);
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
    return numevents;
}
edge *find_edges(configuration *config, io_struct *io, signal_struct *sig, baseline_struct *baseline_stats, bessel *lpfilter, edge *current_edge, edge *head_edge)
{
    fprintf(io->logfile, "<----RUN LOG BEGINS---->\n\n");
    printf("Locating events... \n");
    fprintf(io->logfile, "Locating events...\n ");
    fflush(stdout);

    double baseline;
    double badbaseline = 0;
    double goodbaseline = 0;
    int64_t read;
    int64_t pos;
    int endflag;
    endflag = 0;
    read = 0;
    pos = 0;
    char progressmsg[STRLENGTH];
    time_t start_time;
    time_t curr_time;
    time(&start_time);
    for (pos = config->start; pos < config->finish; pos += read)
    {
        snprintf(progressmsg,STRLENGTH*sizeof(char)," %.2f seconds processed (%d%% good)",(pos-config->start)/(double) config->samplingfreq, pos > config->start ? (int) (100.0 * (pos-config->start - badbaseline)/(double) (pos-config->start)) : 0);
        progressbar(pos-config->start,config->finish-config->start,progressmsg,difftime(time(&curr_time),start_time));
        read = read_current(io->input, sig->signal, sig->rawsignal, pos, intmin(config->readlength,config->finish - pos), config->datatype, config->daqsetup, config->savegain);
        if (read < config->readlength || feof(io->input))
        {
            endflag = 1;
        }
        if (config->usefilter)
        {
            filter_signal(sig->signal, sig->paddedsignal, lpfilter, read);
        }
        gauss_histogram(sig->signal, baseline_stats, read);
        if (isnan(baseline_stats->mean) || isnan(baseline_stats->stdev))
        {
            baseline_stats->mean = 0;
            baseline_stats->stdev = 0;
        }
        baseline = baseline_stats->mean;
        output_baseline_stats(io->baselinefile, baseline_stats, pos, config->samplingfreq);
        if (baseline < config->baseline_min || baseline > config->baseline_max || d_abs(baseline_stats->stdev) < EPS)
        {
            badbaseline += read;
        }
        else
        {
            goodbaseline += read;
            current_edge = detect_edges(sig->signal, baseline, read, current_edge, config->threshold, baseline_stats->stdev, config->hysteresis, pos, config->event_direction);
        }
        if (endflag)
        {
            pos += read;
                break;
        }
        memset(sig->signal,'0',(config->readlength)*sizeof(double));
    }
    snprintf(progressmsg,STRLENGTH*sizeof(char)," %.2f seconds processed (%d%% good)",(pos-config->start)/(double) config->samplingfreq, pos > config->start ? (int) (100.0 * (pos-config->start - badbaseline)/(double) (pos-config->start)) : 100);
    progressbar(pos-config->start,config->finish-config->start,progressmsg,difftime(time(&curr_time),start_time));
    printf("\nRead %g seconds of good baseline\nRead %g seconds of bad baseline\n", goodbaseline/(double) config->samplingfreq, badbaseline / (double) config->samplingfreq);
    fprintf(io->logfile, "\nRead %g seconds of good baseline\nRead %g seconds of bad baseline\n", goodbaseline/(double) config->samplingfreq, badbaseline / (double) config->samplingfreq);

    current_edge = head_edge;

    if (!current_edge || current_edge->type == HEAD)
    {
        printf("No signals found, exiting\n");
        fprintf(io->logfile, "No signals found, exiting\n");
        pause_and_exit(8);
    }

    current_edge = head_edge;
    return current_edge;
}
int64_t get_next_event_start(edge *current_edge)
{
    while (current_edge->type != 0 && current_edge->next) //if for some reason there are multiple of the same type in a row, skip them.
    {
        current_edge = current_edge->next;
    }
    return current_edge->location;
}

int64_t get_next_event(event *current_event, edge *current_edge, int64_t index)
{
    int64_t edges = 0;
    int64_t start;
    int64_t finish;
    while (current_edge->type != 0 && current_edge->next) //if for some reason there are multiple of the same type in a row, skip them.
    {
        current_edge = current_edge->next;
        edges++;
    }
    start = current_edge->location;
    while (current_edge->type != 1 && current_edge->next) //if for some reason there are multiple of the same type in a row, skip them.
    {
        current_edge = current_edge->next;
        edges++;
    }
    finish = current_edge->location;
    if (finish > start)
    {
        current_event = add_event(current_event, start, finish, index, current_edge->local_stdev, current_edge->local_baseline);
    }
    return edges;
}

void calculate_level_noise(event *current, int64_t minpoints)
{
#ifdef DEBUG
    printf("Noise\n");
    fflush(stdout);
#endif // DEBUG
    if (current->type == CUSUM || current->type == STEPRESPONSE)
    {
        cusumlevel *level = current->first_level;
        int64_t start = minpoints;
        while (level)
        {
            if (minpoints + 2 > level->length)
            {
                level->stdev = 0;
            }
            else
            {
                level->stdev = sqrt(signal_variance(&current->signal[start], level->length - minpoints));
            }
            start += level->length;
            level = level->next;
        }
    }
}


void identify_step_events(event *current, int64_t stepfit_samples, int64_t subevent_minpoints, int attempt_recovery)
{
#ifdef DEBUG
    printf("Step Events\n");
    fflush(stdout);
#endif // DEBUG
    if (current->length < stepfit_samples || current->length < subevent_minpoints)
    {
        if (attempt_recovery || stepfit_samples)
        {
            current->type = STEPRESPONSE;
        }
        else
        {
            current->type = TOOSHORT;
        }
    }
}

void filter_long_events(event *current, int64_t maxpoints)
{
#ifdef DEBUG
    printf("Long Filter\n");
    fflush(stdout);
#endif // DEBUG
    if (current->length > maxpoints)
    {
        current->type = TOOLONG;
    }
}

void filter_short_events(event *current, int64_t minpoints)
{
#ifdef DEBUG
    printf("Short Filter\n");
    fflush(stdout);
#endif // DEBUG
    if (current->length < minpoints)
    {
        current->type = TOOSHORT;
    }
}


int64_t average_cusum_levels(event *current, int64_t subevent_minpoints, double cusum_minstep, int attempt_recovery, int64_t padding_wait)
{
#ifdef DEBUG
    printf("Average CUSUM levels\n");
    fflush(stdout);
#endif // DEBUG
    if (current->type == CUSUM)
    {
        cusum_minstep *= current->local_stdev;
        edge *first_edge = current->first_edge;
        edge *current_edge = first_edge;
        double baseline = signal_average(current->signal, current->padding_before);
        double lastlevel = baseline;
        int64_t j, nStates = 0;
        int64_t typeswitch = 0;
        int passflag = 1;
        int64_t anchor = 0;
        int64_t prev_anchor = anchor;
        double average;
        double residual = 0;
        int64_t totallength = current->length + current->padding_before + current->padding_after;
        int64_t levels = 0;
        while (current_edge)
        {
            levels++;
            current_edge = current_edge->next;
        }
        current_edge = first_edge;
        while (current_edge)
        {
            if (nStates <= levels - 1)
            {
                average = signal_average(&current->signal[anchor + subevent_minpoints], current_edge->location - anchor - subevent_minpoints);
            }
            else
            {
                average = signal_average(&current->signal[anchor + padding_wait + subevent_minpoints], current_edge->location - anchor - padding_wait - subevent_minpoints);
            }
            if (d_abs(average - lastlevel) < cusum_minstep && passflag == 0)
            {
                anchor = prev_anchor;
                passflag = 1;
                continue;
            }
            passflag = 0;
            lastlevel = average;
            for (j=anchor; j<current_edge->location; j++)
            {
                current->filtered_signal[j] = average;
            }
            prev_anchor = anchor;
            anchor = current_edge->location;
            current_edge = current_edge->next;
        }
        residual = (current->signal[0]-current->filtered_signal[0])*(current->signal[0]-current->filtered_signal[0]);
        for (j=1; j<totallength; j++)
        {
            if (signum(current->filtered_signal[j] - current->filtered_signal[j-1]) != 0)
            {
                nStates++;
            }
            residual += (current->signal[j]-current->filtered_signal[j])*(current->signal[j]-current->filtered_signal[j]);
        }
        current->residual = sqrt(residual / totallength);
        nStates++;
        if (nStates < 3)
        {
            if (attempt_recovery)
            {
                current->type = STEPRESPONSE;
                typeswitch = 1;
            }
            else
            {
                current->type = BADLEVELS;
            }
        }
        current_edge = first_edge;
        return typeswitch;
    }
    return 0;
}


void populate_event_levels(event *current)
{
#ifdef DEBUG
    printf("Populate Levels\n");
    fflush(stdout);
#endif // DEBUG'
    if (current->type == CUSUM || current->type == STEPRESPONSE)
    {
        int64_t i;
        int64_t eventlength = current->length;
        int64_t padding = current->padding_before + current->padding_after;
        double *filtered_signal = current->filtered_signal;
        double lastlevel = filtered_signal[0];
        int64_t anchor = 0;
        int64_t numlevels = 0;
        current->first_level = initialize_levels();
        cusumlevel *current_level = current->first_level;


        for (i=0; i<eventlength + padding; i++)
        {
            if (signum(lastlevel - filtered_signal[i]) != 0)
            {
                current_level = add_cusum_level(current_level, lastlevel, i-anchor);
                anchor = i;
                lastlevel = filtered_signal[i];
                numlevels++;
            }
        }
        current_level = add_cusum_level(current_level, lastlevel, i-anchor);
        numlevels++;
        current->numlevels = numlevels;

        if (numlevels < 3)
        {
            current->type = BADLEVELS;
        }
    }
}


void event_max_blockage(event *current, double minstep)
{
#ifdef DEBUG
    printf("Max Blockage\n");
    fflush(stdout);
#endif // DEBUG
    if (current->type == CUSUM || current->type == STEPRESPONSE)
    {
        minstep *= current->local_stdev;
        cusumlevel *current_level = current->first_level;
        double baseline = 0.5*(current->baseline_before + current->baseline_after);
        double blockage;
        double maxblockage = 0;
        double minblockage = DBL_MAX;
        int minfound = 0;
        int64_t maxsamples = 0;
        int64_t minsamples = 0;
        while (current_level)
        {
            blockage = d_abs(baseline - current_level->current);
            if (blockage > maxblockage)
            {
                maxblockage = blockage;
                maxsamples = current_level->length;
            }
            if (blockage < minblockage && blockage > minstep)
            {
                minblockage = blockage;
                minsamples = current_level->length;
                minfound = 1;
            }
            current_level = current_level->next;
        }
        current->max_blockage = maxblockage;
        current->max_length = maxsamples;
        current->min_blockage = minblockage;
        current->min_length = minsamples;
        if (!minfound)
        {
            current->type = BADBASELINE;
        }
    }
}


void refine_event_estimates(event *current)
{
#ifdef DEBUG
    printf("Refine\n");
    fflush(stdout);
#endif // DEBUG
    if (current->type == CUSUM || current->type == STEPRESPONSE)
    {
        cusumlevel *level = current->first_level;
        int64_t newstart = current->start + level->length - current->padding_before;
        current->padding_before = level->length;
        while (level)
        {
            if (level->next)
            {
                level = level->next;
            }
            else
            {
                break;
            }
        }
        int64_t newfinish = current->finish - level->length + current->padding_after;
        current->padding_after = level->length;
        current->start = newstart;
        current->finish = newfinish;
        current->length = newfinish - newstart;
    }
}




double get_cusum_threshold(int64_t length, double minthreshold, double maxthreshold, double sigma, double mun)
{
    length *= 2;
    double arlmin;
    double mindif;
    double h;
    int sign;
    int oldsign;
    double arl;
    double threshold = minthreshold;
    arlmin = ARL(length, sigma, mun, minthreshold);
    oldsign = signum(arlmin);
    mindif = d_abs(arlmin);

    for (h = minthreshold; h<maxthreshold; h += 0.5)
    {
        arl = ARL(length, sigma, mun, h);
        sign = signum(arl);
        if (sign != oldsign)
        {
            threshold = h;
            break;
        }
        else if (d_abs(arl) < mindif)
        {
            mindif = d_abs(arl);
            threshold = h;
        }
    }
    return threshold;
}


void cusum(event *current_event, double delta, double minthreshold, double maxthreshold, int64_t subevent_minpoints)
{
#ifdef DEBUG
    printf("cusum\n");
    fflush(stdout);
#endif // DEBUG
    if (current_event->type == CUSUM)
    {
        delta *= current_event->local_stdev;
        double baseline = current_event->local_baseline;
        double *signal = current_event->signal;
        int64_t length = current_event->length + current_event->padding_before + current_event->padding_after;
        current_event->first_edge = initialize_edges();
        edge *first_edge = current_event->first_edge;
        edge *current_edge = first_edge;
        double threshold = minthreshold;
        int64_t anchor = 0;//initial position
        double mean = signal[anchor];//initial mean value guess
        double variance = 0;//initial variance estimation
        double logp = 0;//log-likelihood ratio for positive jumps
        double logn = 0;//log-likelihood ratio for negative jumps
        double stdev = sqrt(signal_variance(signal, current_event->padding_before));
        double jumpsizestdev = delta/stdev;
        int64_t k = 0;//loop index
        int64_t numjumps = 0;
        int64_t jumppos = 0;
        int64_t jumpneg = 0;
        double varM = signal[0];
        double varS = 0;
        double oldVarM;



        threshold = get_cusum_threshold(length, minthreshold, maxthreshold, jumpsizestdev, -1.0*jumpsizestdev/2.0);
        current_event->threshold = threshold;

        double *cpos;//cumulative log-likelihood for positive jumps
        cpos = calloc_and_check(length, sizeof(double),"Cannot allocate cpos");
        cpos[0] = 0;
        double *cneg;//cumulative log-likelihood for negative jumps
        cneg = calloc_and_check(length, sizeof(double),"Cannot allocate cneg");
        cneg[0] = 0;
        double *gpos;//decision function for positive jumps
        gpos = calloc_and_check(length, sizeof(double),"Cannot allocate gpos");
        gpos[0] = 0;
        double *gneg;//decision function for negative jumps
        gneg = calloc_and_check(length, sizeof(double),"Cannot allocate cneg");
        gneg[0] = 0;



        int64_t detected = 0;
        while (k<length-1)
        {
            k++;
            if (k-anchor < subevent_minpoints && signal[k]*signum(signal[k]) > baseline*signum(baseline))
            {
                break;
            }
            mean = ((k-anchor) * mean + signal[k])/(double) (k+1-anchor);  //mean = signal_average(&signal[anchor], k+1-anchor);
            oldVarM = varM;
            varM = varM + (signal[k] - varM) / (double) (k+1-anchor); //variance = signal_variance(&signal[anchor], k+1-anchor);
            varS = varS + (signal[k] - oldVarM) * (signal[k] - varM);
            variance = varS / (double) (k - anchor);
            logp = delta/variance*(signal[k]-mean-delta/2);
            logn = -delta/variance*(signal[k]-mean+delta/2);
            cpos[k] = cpos[k-1] + logp;
            cneg[k] = cneg[k-1] + logn;
            gpos[k] = my_max(gpos[k-1] + logp, 0);
            gneg[k] = my_max(gneg[k-1] + logn, 0);
            if (gpos[k] > threshold || gneg[k] > threshold)
            {
                if (gpos[k] > threshold)
                {
                    jumppos = anchor + locate_min(&cpos[anchor], k+1-anchor);
                    if (jumppos - current_edge->location > subevent_minpoints && length - jumppos > subevent_minpoints)
                    {
                        current_edge = add_edge(current_edge, jumppos, 1, 0, 0);
                        numjumps++;
                    }
                }
                if (gneg[k] > threshold)
                {
                    jumpneg = anchor + locate_min(&cneg[anchor], k+1-anchor);
                    if (jumpneg - current_edge->location > subevent_minpoints && length - jumpneg > subevent_minpoints)
                    {
                        current_edge = add_edge(current_edge, jumpneg, -1, 0, 0);
                        numjumps++;
                    }
                }
                anchor = k;
                mean = 0;
                variance = 0;
                memset(cpos,'0',length*sizeof(double));
                memset(cneg,'0',length*sizeof(double));
                memset(gpos,'0',length*sizeof(double));
                memset(gneg,'0',length*sizeof(double));
                mean = signal[anchor];
                varM = signal[anchor];
                varS = 0;
            }
        }
        current_edge = add_edge(current_edge, length, 1, 0, 0);
        detected++;
        free(cpos);
        free(cneg);
        free(gpos);
        free(gneg);
    }
}

void event_area(event *current_event, double timestep)
{
#ifdef DEBUG
    printf("Area\n");
    fflush(stdout);
#endif // DEBUG
    if (current_event->type == CUSUM || current_event->type == STEPRESPONSE)
    {
        int64_t i;
        double area = 0;
        int64_t padding = current_event->padding_before;
        int64_t length = current_event->length;
        double *signal = current_event->signal;

        double baseline = 0.5 * (current_event->baseline_before + current_event->baseline_after);

        for (i=padding; i<padding+length; i++)
        {
            area += (signal[i] - baseline) * timestep;
        }
        current_event->area = d_abs(area);
        current_event->average_blockage = current_event->area/(current_event->length*timestep);
    }
}


void event_baseline(event *current_event, double baseline_min, double baseline_max)
{
#ifdef DEBUG
    printf("Baseline\n");
    fflush(stdout);
#endif // DEBUG
    if (current_event->type == CUSUM || current_event->type == STEPRESPONSE)
    {
        double stdev;
        double baseline_before = 0;
        int64_t baseline_before_length = 0;
        double baseline_after = 0;
        double *signal = current_event->signal;
        cusumlevel *current_level = current_event->first_level;


        baseline_before = current_level->current;
        baseline_before_length = current_level->length;
        while (current_level->next)
        {
            current_level = current_level->next;
        }
        baseline_after = current_level->current;


        stdev = sqrt(signal_variance(signal, baseline_before_length));
        if (d_abs(baseline_before - baseline_after) > 2.5 * stdev || baseline_before > baseline_max || baseline_after > baseline_max || baseline_before < baseline_min || baseline_after < baseline_min)
        {
            current_event->type = BADBASELINE;
        }
        current_event->baseline_before = baseline_before;
        current_event->baseline_after = baseline_after;

        int64_t i;
        int64_t length = current_event->length + current_event->padding_before + current_event->padding_after;
        double maxdeviation = 0;
        double deviation;
        double baseline = 0.5 * (baseline_before + baseline_after);
        for (i=0; i<length; i++)
        {
            deviation = d_abs(current_event->signal[i] - baseline);
            if (deviation > maxdeviation)
            {
                maxdeviation = deviation;
            }
        }
        current_event->maxdeviation = maxdeviation;

    }
}



void generate_trace(FILE *input, event *current, int datatype, void *rawsignal, FILE *logfile, bessel *lpfilter, int eventfilter, chimera *daqsetup, edge *current_edge, int64_t last_end, int64_t start, int64_t subevent_minpoints, double savegain, int64_t padding_wait)
{
#ifdef DEBUG
    printf("Generate Trace\n");
    fflush(stdout);
#endif // DEBUG
    if (current->type == CUSUM || current->type == STEPRESPONSE)
    {
        int64_t filter_order = 0;
        int64_t filter_padding = 0;
        int64_t padding = intmin(500, current->length);
        if (lpfilter)
        {
            filter_order = lpfilter->order;
            filter_padding = lpfilter->padding;
            padding = 100 * 2.0 / lpfilter->cutoff;
        }

        int64_t position;
        int64_t read;
        int64_t next_start = get_next_event_start(current_edge);
        current->padding_before = padding;
        current->padding_after = padding + padding_wait;
        if (current->start - start < current->padding_before)
        {
            current->padding_before = current->start - start;
        }
        if (current->start - current->padding_before < last_end)
        {
            current->padding_before = current->start - last_end;
        }
        if (current->padding_after + current->finish > next_start && current->finish != next_start)
        {
            current->padding_after = next_start - current->finish;
        }
        position = current->start - current->padding_before;
        if (position > current->start || current->padding_before < subevent_minpoints || current->padding_after < subevent_minpoints || current->start - current->padding_before < last_end)
        {
            current->type = BADPADDING;
            return;
        }


        current->paddedsignal = calloc_and_check(current->length + current->padding_before + current->padding_after + 2*(filter_order + filter_padding),sizeof(double),"Cannot allocate event signal array");
        current->signal = &current->paddedsignal[filter_order + filter_padding];
        current->filtered_signal = calloc_and_check(current->length + current->padding_before + current->padding_after,sizeof(double),"Cannot allocate event filtered signal array");



        if (fseeko64(input,(off64_t) position*2*sizeof(double),SEEK_SET))
        {
            printf("Cannot locate file position at sample %" PRId64 "\n",position);
            fprintf(logfile,"Cannot locate file position at sample %" PRId64 "\n",position);
            fflush(logfile);
            pause_and_exit(17);
        }

        read = read_current(input, current->signal, rawsignal, position, current->length + current->padding_before + current->padding_after, datatype, daqsetup, savegain);


        if (read != current->length + current->padding_before + current->padding_after)
        {
            printf("Unable to read %" PRId64 " samples for event %" PRId64 ": obtained %" PRId64 "\n",current->length + + current->padding_before + current->padding_after,current->index,read);
            fprintf(logfile,"Unable to read %" PRId64 " samples for event %" PRId64 ": obtained %" PRId64 "\n",current->length + + current->padding_before + current->padding_after,current->index,read);
            fflush(logfile);
            current->type = BADTRACE;
        }
        if (eventfilter)
        {
            filter_signal(current->signal, current->paddedsignal, lpfilter, read);
        }
    }
}




edge *detect_edges(double *signal, double baseline, int64_t length, edge *current, double threshold, double stdev, double hysteresis, int64_t position, int event_direction)
{

    int64_t i = 0;
    double sign;
    double down_threshold;
    double up_threshold;
    int state = 0;

    threshold *= stdev;
    hysteresis *= stdev;

    sign = signum(baseline); //get the sign of the average so that we can properly invert the signal
    baseline *= sign;

    if (length <=0 )
    {
        printf("read length is zero or less, cannot process\n");
        return NULL;
    }

    if (event_direction == 0)
    {
        down_threshold = baseline - threshold;//current thresholds for detection of downspikes and upspikes can be different
        up_threshold = baseline - threshold + hysteresis;

        for (i=0; i<length; i++)
        {
            if (signal[i]*sign < down_threshold && state == 0) //if we are open pore state and detect a downspike
            {
                current = add_edge(current, i+position, state, stdev, baseline);
                state = 1;
            }
            else if (signal[i]*sign > up_threshold && state == 1) //blocked state and detect an upspike
            {
                current = add_edge(current, i+position, state, stdev, baseline);
                state = 0;
            }
        }
    }
    else
    {
        up_threshold = baseline + threshold;
        down_threshold = baseline + threshold - hysteresis;//current thresholds for detection of downspikes and upspikes can be different

        for (i=0; i<length; i++)
        {
            if (signal[i]*sign > up_threshold && state == 0) //if we are open pore state and detect a downspike
            {
                current = add_edge(current, i+position, state, stdev, baseline);
                state = 1;
            }
            else if (signal[i]*sign < down_threshold && state == 1) //blocked state and detect an upspike
            {
                current = add_edge(current, i+position, state, stdev, baseline);
                state = 0;
            }
        }
    }
    return current;
}


void gauss_histogram(double *signal, baseline_struct *baseline, int64_t length)
{
    double *histogram = baseline->histogram;
    double baseline_min = baseline->baseline_min;
    double baseline_max = baseline->baseline_max;
    double delta = baseline->delta;
    int64_t numbins = baseline->numbins;
    int64_t i;
    for (i=0; i<numbins; i++)
    {
        histogram[i] = 0;
    }
    for (i=0; i<length; i++)
    {
        if (signal[i] > baseline_min && signal[i] < baseline_max)
        {
            histogram[(int64_t) ((signal[i]-baseline_min)/delta)] += 1;
        }
    }
    fit_gaussian(baseline);
    //printf("\n\ng\t%g\t%g\n\n",baseline->mean, baseline->stdev);
}

void fit_gaussian(baseline_struct *baseline)
{
    double *x = baseline->current;
    double *y = baseline->histogram;
    int64_t numbins = baseline->numbins;
    double xTx[3][3];
    double xTxinv[3][3];
    double xnlny[3];
    double x0 = 0;
    double x1 = 0;
    double x2 = 0;
    double x3 = 0;
    double x4 = 0;
    double lny = 0;
    double xlny = 0;
    double x2lny = 0;
    int64_t minbin = 0;
    int64_t maxbin = numbins;
    int64_t i,sign,max_location;
    double maxval = signal_max(y,numbins);
    int64_t interval;


    i = locate_max(y,numbins);
    max_location = i;
    sign = signum(baseline->current[i]);

    if (sign < 0)
    {
        while (i >= 0)
        {
            if ((y[i] < my_max(1,exp(-4.5)*maxval) && y[i] > 0))
            {
                minbin = i;
                break;
            }
            i--;
        }
        interval = max_location - minbin;
        maxbin = intmin(max_location + (int64_t) (interval/3.0 * sqrt(2.0 * log(2.0))),numbins-1);
    }
    else if (sign > 0)
    {
        while (i < numbins)
        {
            if (y[i] < my_max(1,exp(-4.5)*maxval) && y[i] > 0)
            {
                maxbin = i;
                break;
            }
            i++;
        }
        interval = maxbin - max_location;
        minbin = intmax(0,max_location - (int64_t) (interval/3.0 * sqrt(2.0 * log(2.0))));
    }
    else
    {
        baseline->stdev = 0;
        baseline->mean = 0;
        baseline->amplitude = 0;
        return;
    }

    double rough_mean = x[max_location];
    double rough_stdev = (double) (interval / 3.0 * (x[1]-x[0]));
    double localx, localy;
    for (i=minbin; i<maxbin; i++)
    {
        localy = y[i] / (double) maxval;
        localx = (x[i] - rough_mean)/rough_stdev;
        x0 += localy;
        x1 += localx*localy;
        x2 += localx*localx*localy;
        x3 += localx*localx*localx*localy;
        x4 += localx*localx*localx*localx*localy;
        if (localy > 0)
        {
            lny += log(localy)*localy;
            xlny += localx*log(localy)*localy;
            x2lny += localx*localx*log(localy)*localy;
        }
    }

    xTx[0][0] = x4;
    xTx[0][1] = x3;
    xTx[0][2] = x2;
    xTx[1][0] = x3;
    xTx[1][1] = x2;
    xTx[1][2] = x1;
    xTx[2][0] = x2;
    xTx[2][1] = x1;
    xTx[2][2] = x0;

    xnlny[0] = x2lny;
    xnlny[1] = xlny;
    xnlny[2] = lny;

    invert_matrix(xTx, xTxinv);

    double params[3];
    params[0] = xTxinv[0][0]*xnlny[0] + xTxinv[0][1]*xnlny[1] + xTxinv[0][2]*xnlny[2];
    params[1] = xTxinv[1][0]*xnlny[0] + xTxinv[1][1]*xnlny[1] + xTxinv[1][2]*xnlny[2];
    params[2] = xTxinv[2][0]*xnlny[0] + xTxinv[2][1]*xnlny[1] + xTxinv[2][2]*xnlny[2];

    baseline->stdev = sqrt(-1.0/(2*params[0]));
    baseline->mean = baseline->stdev*baseline->stdev*params[1];
    baseline->amplitude = exp(params[2] + baseline->mean * baseline->mean/(2.0*baseline->stdev*baseline->stdev));
    baseline->stdev *= rough_stdev;
    baseline->mean += rough_mean;
    baseline->amplitude *= maxval;
}





