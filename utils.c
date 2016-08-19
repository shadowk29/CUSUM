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
#include"utils.h"

double signal_max(double *signal, int64_t length)
{
    int64_t i;
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
double signal_min(double *signal, int64_t length)
{
    int64_t i;
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

double signal_average(double *signal, int64_t length)
{
    int64_t i;
    double average;
    average = 0;
    for (i=0; i<length; i++)
    {
        average += signal[i];
    }
    return average/length;
}

double signal_extreme(double *signal, int64_t length, double sign)
{
    int64_t i;
    double tempmax;
    tempmax = 0;
    for (i=0; i<length; i++)
    {
        if (signal[i]*sign > tempmax)
        tempmax = signal[i]*sign;
    }
    return tempmax;
}

double signal_variance(double *signal, int64_t length)
{
    if (length < 2)
    {
        printf("Cannot calculate variance with less than 2 samples\n");
        return 0;
    }
    int64_t i;
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

FILE *fopen64_and_check(const char *fname, const char *mode, int error)
{
    FILE *buffer;
    if ((buffer=fopen64(fname,mode))==NULL)
    {
        printf("Cannot open file %s\n",fname);
        exit(error);
    }
    return buffer;
}

void *calloc_and_check(size_t num, size_t size, char *msg)
{
    void *block;
    if ((block=calloc(num,size))==NULL)
    {
        printf("%s\n",msg);
        exit(1);
    }
    return block;
}


int64_t count_edges(edge *current_edge)
{
    int64_t count = 0;
    while (current_edge)
    {
        count++;
        current_edge = current_edge->next;
    }
    return count;
}

inline void progressbar(int64_t pos, int64_t finish, const char *msg, double elapsed)
{
    double ratio = pos/(double)finish;
    double remaining;

    if (pos == 0)
    {
        remaining = 0;
    }
    else
    {
        remaining = elapsed * (double) (finish - pos) / (double) pos;
    }


    int64_t hours = (int64_t) remaining / 3600;
    int64_t rhours = (int64_t) remaining % 3600;
    int64_t minutes = rhours / 60;
    int64_t seconds = rhours % 60;
    printf("%3d%%\t", (int)(ratio*100) );
    printf("%02"PRId64":%02"PRId64":%02"PRId64" remaining\t%s       \r",hours,minutes,seconds,msg);
    fflush(stdout);
}

int64_t get_filesize(FILE *input, int datatype)
{
    int64_t length;
    if (datatype != 0)
    {
        fseeko64(input, 0, SEEK_END);
        length = ftello64(input);
        fseeko64(input, 0, SEEK_SET);
        return length / (datatype / 8 * 2);
    }
    else
    {
        fseeko64(input, 0, SEEK_END);
        length = ftello64(input);
        fseeko64(input, 0, SEEK_SET);
        return length / 2;
    }

}


inline int signum(double num)
{
    return (EPS<num)-(num<-EPS);
}
inline double my_min(double a, double b)
{
    return a < b ? a : b;
}

inline double my_max(double a, double b)
{
    return a > b ? a : b;
}

inline double d_abs(double num)
{
    return num > 0 ? num : -num;
}

inline int64_t intmin(int64_t a, int64_t b)
{
    return a < b ? a : b;
}
inline int64_t intmax(int64_t a, int64_t b)
{
    return a > b ? a : b;
}

edge *initialize_edges(void)
{
    edge *head;
    head = calloc_and_check(1,sizeof(edge),"Cannot initialize edge list");
    head->next = NULL;
    head->location = 0;
    head->type = HEAD;
    return head;
}


edge *add_edge(edge *current, int64_t location, int type)
{
    if (current->type == HEAD) //if we are adding a new node to the head node that hasn't been filled yet
    {
        current->location = location;
        current->type = type;
    }
    else
    {//if the current node is filled with useful information and we actually need more memory
        current->next = calloc_and_check(1,sizeof(edge),"Cannot allocate next edge");
        current->next->location = location;
        current->next->type = type;
        current->next->next = NULL;
        current = current->next;
    }
    return current;
}

void free_edges(edge *current)
{
    edge *temp;
    while (current)
    {
        temp = current->next;
        free(current);
        current = temp;
    }
}

void free_levels(cusumlevel *current)
{
    cusumlevel *temp;
    while (current)
    {
        temp = current->next;
        free(current);
        current = temp;
    }
}


event *initialize_events(void)
{
    event *head;
    head = calloc_and_check(1,sizeof(event),"Cannot allocate head event");
    head->type = 0;
    head->threshold = 0;
    head->rc1 = 0;
    head->rc2 = 0;
    head->index = HEAD;
    head->signal = NULL;
    head->filtered_signal = NULL;
    head->first_edge = NULL;
    head->first_level = NULL;
    return head;
}

event *add_event(event *current, int64_t start, int64_t finish, int64_t index)
{
    current->type = 0;
    current->index = index;
    current->start = start;
    current->finish = finish;
    current->length = finish-start;
    current->threshold = 0;
    current->rc1 = 0;
    current->rc2 = 0;
    current->first_edge = NULL;
    current->first_level = NULL;
    current->signal = NULL;
    current->filtered_signal = NULL;
    return current;
}


void free_single_event(event *current)
{
#ifdef DEBUG
    printf("Free event\n");
    fflush(stdout);
#endif // DEBUG
    if (current->signal)
    {
        free(current->signal);
    }
    if (current->filtered_signal)
    {
        free(current->filtered_signal);
    }
    if (current->first_edge)
    {
        free_edges(current->first_edge);
    }
    if (current->first_level)
    {
        free_levels(current->first_level);
    }
}


cusumlevel *initialize_levels(void)
{
    cusumlevel *head;
    head=calloc_and_check(1,sizeof(cusumlevel),"Cannot allocate head level");
    head->current = 0;
    head->length = 0;
    head->next = NULL;
    return head;
}

cusumlevel *add_cusum_level(cusumlevel *lastlevel, double current, int64_t length)
{
    cusumlevel *temp;
    if (lastlevel && lastlevel->length > 0)
    {
        lastlevel->next=calloc_and_check(1,sizeof(cusumlevel),"Cannot allocate next level");
        lastlevel->next->current = current;
        lastlevel->next->length = length;
        lastlevel->next->next = NULL;
        temp = lastlevel->next;
    }
    else
    {
        lastlevel->current = current;
        lastlevel->length = length;
        lastlevel->next = NULL;
        temp = lastlevel;
    }
    return temp;
}


double ARL(int64_t length, double sigma, double mun, double h)
{
    return (exp(-2.0*mun*(h/sigma+1.166))-1.0+2.0*mun*(h/sigma+1.166))/(2.0*mun*mun)-(double) length;
}



