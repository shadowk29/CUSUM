#include"detector.h"
#include"io.h"
#include<inttypes.h>
#include<stdint.h>
#include<string.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<float.h>


void filter_event_length(event *current, uint64_t maxpoints, uint64_t minpoints)
{
    uint64_t toolong = 0;
    uint64_t tooshort = 0;
    while (current)
    {
        if (current->length > maxpoints)
        {
            current->type = TOOLONG;
            toolong++;
        }
        else if (current->length < minpoints)
        {
            current->type = TOOSHORT;
            tooshort++;
        }
        current = current->next;
    }
    printf("\n%"PRIu64" events were too long\n%"PRIu64" events were too short\n",toolong,tooshort);
}

void count_all_levels(event *current)
{
    int numlevels;
    while (current)
    {
        if (current->type == 0)
        {
            numlevels = count_levels(current);
            current->numlevels = numlevels;
        }

        current = current->next;
    }
}


int count_levels(event *current)
{
    uint64_t i;
    int numlevels = 1;
    double lastlevel;
    lastlevel = current->filtered_signal[0];
    for (i=0; i<current->length + 2*current->padding; i++)
    {
        if (signum(current->filtered_signal[i]-lastlevel) != 0)
        {
            numlevels++;
            lastlevel = current->filtered_signal[i];
        }
    }
    return numlevels;
}

void assign_cusum_levels(event *current, uint64_t subevent_minpoints)
{
    while (current)
    {
        if (current->type == 0)
        {
            average_cusum_levels(current, subevent_minpoints);
        }
        current = current->next;
    }
}

void average_cusum_levels(event *current, uint64_t subevent_minpoints)
{
    edge *first_edge = current->first_edge;
    edge *current_edge = first_edge;
    double blockage;
    double maxblockage = 0;
    double baseline = 0.5 * (current->baseline_before + current->baseline_after);
    uint64_t j;
    uint64_t anchor = 0;
    double average;
    while (current_edge)
    {
        fflush(stdout);
        average = signal_average(&current->signal[anchor + subevent_minpoints], current_edge->location - anchor - subevent_minpoints);
        fflush(stdout);
        blockage = d_abs(average - baseline);
        if (blockage > maxblockage)
        {
            maxblockage = blockage;
        }
        fflush(stdout);
        for (j=anchor; j<current_edge->location; j++)
        {
            current->filtered_signal[j] = average;
        }
        fflush(stdout);
        anchor = current_edge->location;
        fflush(stdout);
        current_edge = current_edge->next;
        fflush(stdout);
    }
    current_edge = first_edge;
    current->max_blockage = maxblockage;
    current->baseline_before = current->filtered_signal[0];
    current->baseline_after  = current->filtered_signal[current->length + 2 * current->padding - 1];
}

void detect_subevents(event *current_event, double delta, double minthreshold, double maxthreshold, uint64_t subevent_minpoints)
{
    while (current_event)
    {
        if (current_event->type == 0)
        {
            cusum(current_event, delta, minthreshold, maxthreshold, subevent_minpoints);
        }
        current_event = current_event->next;
    }
}

uint64_t locate_min(double *signal, uint64_t length)
{
    double minval = signal[0];
    uint64_t location = 0;
    uint64_t i;
    for (i=0; i<length; i++)
    {
        if (signal[i] < minval)
        {
            minval = signal[i];
            location = i;
        }
    }
    return location;
}



double get_cusum_threshold(uint64_t length, double minthreshold, double maxthreshold, double sigma, double mun)
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


void cusum(event *current_event, double delta, double minthreshold, double maxthreshold, uint64_t subevent_minpoints)
{
    double *signal = current_event->signal;
    uint64_t length = current_event->length + 2 * current_event->padding;
    current_event->first_edge = initialize_edges();
    edge *first_edge = current_event->first_edge;
    edge *current_edge = first_edge;
    double threshold = minthreshold;
    uint64_t anchor = 0;//initial position
    double mean = signal[anchor];//initial mean value guess
    double variance = 0;//initial variance estimation
    double logp = 0;//log-likelihood ratio for positive jumps
    double logn = 0;//log-likelihood ratio for negative jumps
    double stdev = sqrt(signal_variance(signal, current_event->padding));
    double jumpsizestdev = delta/stdev;
    uint64_t k = 0;//loop index
    uint64_t numjumps = 0;
    uint64_t jumppos = 0;
    uint64_t jumpneg = 0;
    double varM = signal[0];
    double varS = 0;
    double oldVarM;

    threshold = get_cusum_threshold(length, minthreshold, maxthreshold, jumpsizestdev, -1.0*jumpsizestdev/2.0);
    current_event->threshold = threshold;

    double *cpos;//cumulative log-likelihood for positive jumps
    if ((cpos = calloc(length, sizeof(double)))==NULL)
    {
        printf("Cannot allocate cpos\n");
        fflush(stdout);
        abort();
    }
    cpos[0] = 0;
    double *cneg;//cumulative log-likelihood for negative jumps
    if ((cneg = calloc(length, sizeof(double)))==NULL)
    {
        printf("Cannot allocate cneg\n");
        fflush(stdout);
        abort();
    }
    cneg[0] = 0;
    double *gpos;//decision function for positive jumps
    if ((gpos = calloc(length, sizeof(double)))==NULL)
    {
        printf("Cannot allocate gpos\n");
        fflush(stdout);
        abort();
    }
    gpos[0] = 0;
    double *gneg;//decision function for negative jumps
    if ((gneg = calloc(length, sizeof(double)))==NULL)
    {
        printf("Cannot allocate gneg\n");
        fflush(stdout);
        abort();
    }
    gneg[0] = 0;




    while (k<length-1)
    {
        k++;
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
                    current_edge = add_edge(current_edge, jumppos, 1);
                    numjumps++;
                }
            }
            if (gneg[k] > threshold)
            {
                jumpneg = anchor + locate_min(&cneg[anchor], k+1-anchor);
                if (jumpneg - current_edge->location > subevent_minpoints && length - jumpneg > subevent_minpoints)
                {
                    current_edge = add_edge(current_edge, jumpneg, -1);
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
    current_edge = add_edge(current_edge, length, 1);
    free(cpos);
    free(cneg);
    free(gpos);
    free(gneg);
}

void event_area(event *current_event, double timestep)
{
    uint64_t i;
    double area = 0;
    uint64_t padding = current_event->padding;
    uint64_t length = current_event->length;
    double *signal = current_event->signal;

    double baseline = 0.5 * (current_event->baseline_before + current_event->baseline_after);

    if (current_event->type != 0)
    {
        current_event->area = 0;
    }
    else
    {
        for (i=padding; i<padding+length; i++)
        {
            area += (signal[i] - baseline) * timestep;
        }
        current_event->area = d_abs(area);
        current_event->average_blockage = current_event->area/(current_event->length*timestep);
    }
}

void assign_event_areas(event *current_event, double timestep)
{
    while (current_event != NULL)
    {
        if (current_event->type == 0)
        {
            event_area(current_event, timestep);
        }
        current_event = current_event->next;
    }
}

void assign_event_baselines(event *current_event)
{
    while (current_event != NULL)
    {
        if (current_event->type == 0)
        {
            event_baseline(current_event);
        }
        current_event = current_event->next;
    }
}

void event_baseline(event *current_event)
{
    uint64_t i;
    double baseline_before = 0;
    double baseline_after = 0;
    uint64_t padding = current_event->padding;
    uint64_t length = current_event->length;
    double *signal = current_event->signal;
    double stdev;

    for (i=0; i<padding; i++)
    {
        baseline_before += signal[i];
        baseline_after += signal[i+length+padding];
    }
    baseline_before = baseline_before / padding;
    baseline_after = baseline_after / padding;

    stdev = 0.5 * (sqrt(signal_variance(signal, padding)) + sqrt(signal_variance(&signal[length+padding], padding)));
    if (d_abs(baseline_before - baseline_after) > 3 * stdev)
    {
        current_event->type = BADBASELINE;
    }
    current_event->baseline_before = baseline_before;
    current_event->baseline_after = baseline_after;
}

void populate_event_traces(FILE *input, event *current_event)
{
    while (current_event != NULL)
    {
        if (current_event->type == 0)
        {
            generate_trace(input, current_event);
        }
        current_event = current_event->next;
    }
}



void generate_trace(FILE *input, event *current)
{
    uint64_t padding;
    uint64_t position;
    uint64_t read;

    if (!current || current->length == 0 || current->index == HEAD)
    {
        return;
    }


    padding = intmax(current->length/2,100);
    if (current->index == 0) //for the first event, we need to make sure we don't overshoot the start of the file
    {
        if (padding > current->start)
        {
            padding = current->start;
        }
    }
    else
    {
        if (current->prev)
        {
            if (padding > current->start - current->prev->finish) //if padding would include some of the previous event, pare it down
            {
                padding = current->start - current->prev->finish;
            }
        }
        if (current->next)
        {
            if (padding > current->next->start - current->finish)
            {
                padding = current->next->start - current->finish;
            }
        }
    }

    position = current->start - padding;
    if (position > current->start)
    {
        printf("Attempting to access negative file index, increase your start time past %" PRIu64 "\n",current->start);
        abort();
    }
    current->padding = padding;


    if ((current->signal = calloc(current->length + 2*padding,sizeof(double)))==NULL)
    {
        printf("Cannot allocate trace array\n");
        return;
    }
    if ((current->filtered_signal = calloc(current->length + 2*padding,sizeof(double)))==NULL)
    {
        printf("Cannot allocate filtered trace array\n");
        return;
    }


    if (fseeko64(input,(off64_t) position*2*sizeof(double),SEEK_SET))
    {
        printf("Cannot location file position at sample %" PRIu64 "\n",position);
        return;
    }

    read = read_current(input, current->signal, position, current->length + 2*padding);
    if (read != current->length + 2*padding)
    {
        printf("Unable to read %" PRIu64 " samples for event %" PRId64 ": obtained %" PRId64 "\n",current->length + 2*padding,current->index,read);
        current->type = -2;
    }
}



event *process_edges(edge *current_edge, event *current_event)
{
    uint64_t start = 0;
    uint64_t finish = 0;
    while (current_edge)
    {
        while (current_edge->type != 0 && current_edge->next) //if for some reason there are multiple of the same type in a row, skip them.
        {
            current_edge = current_edge->next;
        }
        start = current_edge->location;
        while (current_edge->type != 1 && current_edge->next) //if for some reason there are multiple of the same type in a row, skip them.
        {
            current_edge = current_edge->next;
        }
        finish = current_edge->location;
        if (finish > start)
        {
            current_event = add_event(current_event, start, finish);
        }
        current_edge = current_edge->next;
    }
    return current_event;
}



edge *detect_edges(double *signal, double baseline, uint64_t length, edge *current, double threshold, double hysteresis, uint64_t position, int event_direction)
{
    uint64_t i = 0;
    double sign;
    double down_threshold;
    double up_threshold;
    int state = 0;

    sign = signum(baseline); //get the sign of the average so that we can properly invert the signal
    baseline *= sign;

    if (length <=0 )
    {
        printf("read length is zero or less, cannot process\n");
        return NULL;
    }

    if (event_direction == 0)
    {
        down_threshold = baseline*(1 - threshold);//current thresholds for detection of downspikes and upspikes can be different
        up_threshold = baseline*(1 - threshold + hysteresis);

        for (i=0; i<length; i++)
        {
            if (signal[i]*sign < down_threshold && state == 0) //if we are open pore state and detect a downspike
            {
                current = add_edge(current, i+position, state);
                state = 1;
            }
            else if (signal[i]*sign > up_threshold && state == 1) //blocked state and detect an upspike
            {
                current = add_edge(current, i+position, state);
                state = 0;
            }
        }
    }
    else
    {
        up_threshold = baseline*(1 + threshold);
        down_threshold = baseline*(1 + threshold - hysteresis);//current thresholds for detection of downspikes and upspikes can be different
        for (i=0; i<length; i++)
        {
            if (signal[i]*sign > up_threshold && state == 0) //if we are open pore state and detect a downspike
            {
                current = add_edge(current, i+position, state);
                state = 1;
            }
            else if (signal[i]*sign < down_threshold && state == 1) //blocked state and detect an upspike
            {
                current = add_edge(current, i+position, state);
                state = 0;
            }
        }
    }
    return current;
}


double build_histogram(double *signal, histostruct *histogram, uint64_t length, double delta, double baseline_max, double baseline_min)
{
    double maximum = signal_max(signal, length);
    double minimum = signal_min(signal, length);

    double range = maximum - minimum;
    uint64_t numbins = range / delta;
    double baseline = 0;
    uint64_t i,j;
    int sign;
    int newsign;
    double tempbase = 0;
    uint64_t average;

    for (i=0; i<histogram->numbins; i++)
    {
        for (j=0; j<3; j++)
        {
            histogram->histogram[i][j] = 0;
        }
    }

    if (histogram->histogram == NULL || histogram->numbins < numbins)
    {
        if ((histogram->histogram = realloc(histogram->histogram,numbins*sizeof(double *)))==NULL)
        {
            printf("Cannot allocate level 1 with %" PRIu64" using range: %.2lf and delta: %.2lf, min: %.2lf and max: %.2lf\n",numbins, range, delta, minimum, maximum);
            abort();
        }
        for (i=histogram->numbins; i<numbins; i++)
        {
            if ((histogram->histogram[i] = calloc(3,sizeof(double)))==NULL)
            {
                printf("Cannot allocate level 2\n");
                abort();
            }
        }
        histogram->numbins = numbins;
    }


    //populate the current level histogram for the data chunk
    if (length <= 0)
    {
        printf("Final sample interval contains no points\n");
        return 0;
    }
    for (i=0; i<length; i++)
    {
        histogram->histogram[(int) my_min(numbins-1, (int) ((signal[i]-minimum)/delta))][0] += 1;
    }



    //get the first derivative of the histogram
    histogram->histogram[0][FIRST_DERIV] = (histogram->histogram[1][HISTOGRAM]-histogram->histogram[0][HISTOGRAM])/delta; //first derivative uses endpoint only
    histogram->histogram[numbins-1][FIRST_DERIV] = (histogram->histogram[numbins-1][HISTOGRAM]-histogram->histogram[numbins-2][HISTOGRAM])/delta; //first derivative uses endpoint only
    for (i=1; i<numbins-1; i++)
    {
        histogram->histogram[i][FIRST_DERIV] = (histogram->histogram[i+1][HISTOGRAM]-histogram->histogram[i-1][HISTOGRAM])/(2.0*delta);
    }


    //get the second derivative of the histogram
    histogram->histogram[0][SCND_DERIV] = (histogram->histogram[1][HISTOGRAM]-2*histogram->histogram[0][HISTOGRAM])/(delta*delta); //second derivative uses endpoint only
    histogram->histogram[numbins-1][SCND_DERIV] = (-2*histogram->histogram[numbins-1][HISTOGRAM]+histogram->histogram[numbins-2][HISTOGRAM])/(delta*delta); //second derivative uses endpoint only
    for (i=1; i<numbins-1; i++)
    {
        histogram->histogram[i][SCND_DERIV] = (histogram->histogram[i+1][HISTOGRAM]-2*histogram->histogram[i][HISTOGRAM]+histogram->histogram[i-1][HISTOGRAM])/(delta*delta);
    }

    sign = signum(histogram->histogram[0][FIRST_DERIV]);
    average = length/numbins;
    for (i=0; i<numbins; i++) //we ignore the endpoints
    {
        newsign = signum(histogram->histogram[i][FIRST_DERIV]);
        if (newsign != sign)//derivative crosses 0
        {
            sign = newsign;
            if (histogram->histogram[i][SCND_DERIV]<0 && histogram->histogram[i][HISTOGRAM]>average/2)
            {
                tempbase = (i+0.5)*delta+minimum;
            }
            else if (histogram->histogram[intmax(i-1,0)][SCND_DERIV]<0 && histogram->histogram[intmax(i-1,0)][HISTOGRAM]>average/2)//check neighbours to find the negative second derivative.
            {
                tempbase = (intmax(i-1,0)+0.5)*delta+minimum;
            }
            else if (histogram->histogram[intmin(i+1,numbins-1)][SCND_DERIV]<0 && histogram->histogram[intmin(i+1,numbins-1)][HISTOGRAM]>average/2)//check neighbours to find the negative second derivative.
            {
                tempbase = (intmin(i+1,numbins-1)+0.5)*delta+minimum;
            }
            if (tempbase*signum(tempbase) > baseline*signum(baseline) && tempbase < baseline_max && tempbase > baseline_min)
            {
                baseline = tempbase;
            }
        }
    }

    histogram->offset = minimum;
    histogram->delta = delta;

    /*char filename[1024];
    sprintf(filename, "output/histo_%" PRIu64".dat",pos);
    print_histogram(filename, histogram);*/

    printf("Baseline detected was %g\n",baseline);
    return baseline;
}

double signal_max(double *signal, uint64_t length)
{
    uint64_t i;
    double maximum = signal[0];
    for (i=0; i<length; i++)
    {
        if (signal[i] > maximum)
        {
            maximum = signal[i];
        }
    }
    return maximum;
}
//get the smallest (most negative) value for the signal
double signal_min(double *signal, uint64_t length)
{
    uint64_t i;
    double minimum = signal[0];
    for (i=0; i<length; i++)
    {
        if (signal[i] < minimum)
        {
            minimum = signal[i];
        }
    }
    return minimum;
}

double signal_average(double *signal, uint64_t length)
{
    uint64_t i;
    double average;
    average = 0;
    for (i=0; i<length; i++)
    {
        average += signal[i];
    }
    return average/length;
}

double signal_extreme(double *signal, uint64_t length, double sign)
{
    uint64_t i;
    double tempmax;
    tempmax = 0;
    for (i=0; i<length; i++)
    {
        if (signal[i]*sign > tempmax)
        tempmax = signal[i]*sign;
    }
    return tempmax;
}

double signal_variance(double *signal, uint64_t length)
{
    if (length < 2)
    {
        printf("Cannot calculate variance with less than 2 samples\n");
        return 0;
    }
    uint64_t i;
    double variance;
    variance = 0;
    double average;
    average = signal_average(signal, length);
    for (i=0; i<length; i++)
    {
        variance += (signal[i] - average)*(signal[i] - average);
    }
    variance = variance / (length - 1);
    return variance;
}
