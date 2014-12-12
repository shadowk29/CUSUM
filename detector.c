#include"detector.h"
#include"io.h"
#include<inttypes.h>
#include<stdint.h>
#include<string.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<float.h>

void event_area(event *current_event, double timestep)
{
    uint64_t i;
    double area = 0;
    uint64_t padding = current_event->padding;
    uint64_t length = current_event->length;
    double *signal = current_event->signal;

    double baseline = 0.5 * (current_event->baseline_before + current_event->baseline_after);

    if (current_event->type == -1)
    {
        current_event->area = 0;
    }
    else
    {
        for (i=padding; i<padding+length; i++)
        {
            area += (signal[i] - baseline) * timestep;
        }
        current_event->area = area;
    }
}

void assign_event_areas(event *current_event, double timestep)
{
    while (current_event != NULL)
    {
        event_area(current_event, timestep);
        current_event = current_event->next;
    }
}

void assign_event_baselines(event *current_event)
{
    while (current_event != NULL)
    {
        event_baseline(current_event);
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
        current_event->type = -1;
    }
    current_event->baseline_before = baseline_before;
    current_event->baseline_after = baseline_after;
}

void populate_event_traces(FILE *input, event *current_event, uint64_t order)
{
    while (current_event != NULL)
    {
        generate_trace(input, current_event, order);
        current_event = current_event->next;
    }
}



void generate_trace(FILE *input, event *current, uint64_t order)
{
    uint64_t padding;
    uint64_t position;
    int64_t read;

    if (!current || current->length == 0 || current->index == -1000)
    {
        return;
    }


    padding = intmax(current->length/2,100);
    if (padding < order) //padding needs to be at least big enough for the filter to work
    {
        padding = order;
    }
    if (current->index == 0) //for the first event, we need to make sure we don't overshoot the start of the file
    {
        if (padding > current->start)
        {
            printf("Reducing padding because of proximity to file start\n");
            padding = current->start;
        }
    }
    else
    {
        if (!current->prev) //make sure there is a previous one before we do anything with it
        {
            printf("previous struct is null\n");
        }
        else
        {
            if (padding > current->start - current->prev->finish) //if padding would include some of the previous event, pare it down
            {
                printf("Reducing padding for event %" PRId64" because of proxmity to previous event\n",current->index);
                padding = current->start - current->prev->finish;
            }
        }
        if (!current->next)
        {
            printf("Next struct is null\n");
        }
        else if (padding > current->next->start - current->finish) // if padding would include some of the next event, pare it down
        {
            printf("Reducing padding for event %" PRId64" because of proxmity to next event\n",current->index);
            padding = current->next->start - current->finish;
        }
    }
    if (padding < order) //padding needs to be at least big enough for the filter to work
    {
        printf("Cannot pad event %" PRId64" with enough points to filter\n Reduce the filter order to less than %" PRIu64 " and try again\n", current->index, padding);
        abort();
    }

    position = current->start - padding;
    if (position < 0)
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



edge *detect_edges(double *signal, uint64_t padding, double baseline, uint64_t length, edge *current, double threshold, double hysteresis, uint64_t position, int event_direction)
{
    uint64_t i = 0;
    double sign;
    double down_threshold;
    double up_threshold;
    int state = 0;

    sign = signum(baseline); //get the sign of the average so that we can properly invert the signal
    baseline *= sign;

    if (length < 2*padding)
    {
        printf("read length < 2*padding, cannot process edges over negative number of samples\n");
        return NULL;
    }

    if (event_direction == 0)
    {
        down_threshold = baseline*(1 - threshold);//current thresholds for detection of downspikes and upspikes can be different
        up_threshold = baseline*(1 - threshold + hysteresis);

        for (i=padding; i<length-padding; i++)
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
        for (i=padding; i<length-padding; i++)
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


double build_histogram(double *signal, histostruct *histogram, uint64_t length, double delta, uint64_t padding, uint64_t pos, double baseline_max, double baseline_min)
{
    double maximum = signal_max(signal, length, padding);
    double minimum = signal_min(signal, length, padding);

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
    if (length <= 2*padding)
    {
        printf("Final sample interval contains no points: increase your read length to more than twice your filter padding\n");
        return 0;
    }
    for (i=padding; i<length-padding; i++)
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

double signal_max(double *signal, uint64_t length, uint64_t padding)
{
    uint64_t i;
    double maximum = signal[padding];
    for (i=padding; i<length-padding; i++)
    {
        if (signal[i] > maximum)
        {
            maximum = signal[i];
        }
    }
    return maximum;
}
//get the smallest (most negative) value for the signal
double signal_min(double *signal, uint64_t length, uint64_t padding)
{
    uint64_t i;
    double minimum = signal[padding];
    for (i=padding; i<length-padding; i++)
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
