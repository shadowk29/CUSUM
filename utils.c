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

inline void progressbar(uint64_t pos, uint64_t finish)
{
    // Calculuate the ratio of complete-to-incomplete.
    double ratio = pos/(double)finish;
    int   c     = (int) (ratio * 40)+1;

    // Show the percentage complete.
    printf("%3d%% [", (int)(ratio*100+1) );

    // Show the load bar.
    int i;
    for (i=0; i<c; i++)
       printf("=");

    for (i=c; i<40; i++)
       printf(" ");

    // ANSI Control codes to go back to the
    // previous line and clear it.
    printf("]\r");
    fflush(stdout);
}

uint64_t get_filesize(FILE *input, int datatype)
{
    uint64_t length;
    fseeko64(input, 0, SEEK_END);
    length = ftello64(input);
    fseeko64(input, 0, SEEK_SET);
    return length / (datatype / 8 * 2);
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
    if ((head = malloc(sizeof(edge)))==NULL)
    {
        printf("Cannot allocate head node\n");
        abort();
    }
    head->next = NULL;
    head->location = 0;
    head->type = HEAD;
    return head;
}

event *delete_bad_events(event *head)
{
    event *newhead;
    event *current;
    event *tempnext;
    event *tempprev;
    event *temp;
    current = head;
    newhead = head;
    while (current)
    {
        if (current->type != 0)
        {
           tempnext = current->next;
           tempprev = current->prev;
           temp = current;
           current = current->next;
           free_single_event(temp);
           if (tempprev && tempnext)
           {
               tempprev->next = tempnext;
               tempnext->prev = tempprev;
           }
           else if (!tempprev && tempnext) //head of the list
           {
               newhead = tempnext;
               newhead->index = 0;
               tempnext->prev = NULL;
           }
           else if (tempprev && !tempnext) //end of the list
           {
               tempprev->next = NULL;
           }
           else if (!tempprev && !tempnext) //singleton list, odd case
           {
               newhead = NULL;
           }
        }
        else
        {
            current = current->next;
        }
    }
    return newhead;
}

edge *add_edge(edge *current, uint64_t location, int type)
{
    if (current->type == HEAD) //if we are adding a new node to the head node that hasn't been filled yet
    {
        current->location = location;
        current->type = type;
    }
    else
    {//if the current node is filled with useful information and we actually need more memory
        if ((current->next = malloc(sizeof(edge)))==NULL)
        {
            printf("Cannot allocate new node after %" PRIu64 "\n",current->location);
            abort();
        }
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
    if ((head = malloc(sizeof(event)))==NULL)
    {
        printf("Cannot allocate head node for event list\n");
        abort();
    }
    head->type = 0;
    head->threshold = 0;
    head->rc1 = 0;
    head->rc2 = 0;
    head->index = HEAD;
    head->signal = NULL;
    head->filtered_signal = NULL;
    head->first_edge = NULL;
    head->first_level = NULL;
    head->next = NULL;
    head->prev = NULL;
    return head;
}

event *add_event(event *current, uint64_t start, uint64_t finish)
{
    if (current->index == HEAD) //if we are adding a new node to the head node that hasn't been filled yet
    {
        current->index = 0;
        current->start = start;
        current->finish = finish;
        current->length = finish-start;
        current->threshold = 0;
        current->rc1 = 0;
        current->rc2 = 0;
        current->first_edge = NULL;
        current->first_level = NULL;
        current->next = NULL;
        current->prev = NULL;
    }
    else
    {//if the current node is filled with useful information and we actually need more memory
        if ((current->next = malloc(sizeof(event)))==NULL)
        {
            printf("Cannot allocate new event node after %" PRIu64 "\n",current->index);
            abort();
        }
        current->next->type = 0;
        current->next->first_edge = NULL;
        current->next->first_level = NULL;
        current->next->next = NULL;
        current->next->signal = NULL;
        current->filtered_signal = NULL;
        current->next->prev = current;
        current->next->index = current->index + 1;
        current->next->start = start;
        current->next->finish = finish;
        current->next->length = finish-start;
        current = current->next;
        current->threshold = 0;
        current->rc1 = 0;
        current->rc2 = 0;
    }
    return current;
}

void free_events(event *current)
{
    while (current->next)
    {
        current = current->next;
        free_single_event(current->prev);
    }
    free_single_event(current);
}

void free_single_event(event *current)
{
    if (current)
    {
        if (current->signal)
        {
            free(current->signal);
        }
        if (current->filtered_signal)
        {
            free(current->filtered_signal);
        }
        free_edges(current->first_edge);
        free_levels(current->first_level);
    }
    free(current);
}


cusumlevel *initialize_levels(void)
{
    cusumlevel *head;
    if ((head=malloc(sizeof(cusumlevel)))==NULL)
    {
        printf("Cannot allocate head level node\n");
        abort();
    }
    head->current = 0;
    head->length = 0;
    head->next = NULL;
    return head;
}

cusumlevel *add_cusum_level(cusumlevel *lastlevel, double current, uint64_t length)
{
    cusumlevel *temp;
    if (lastlevel && lastlevel->length > 0)
    {
        if ((lastlevel->next=malloc(sizeof(cusumlevel)))==NULL)
        {
            printf("Cannot allocate level node\n");
            abort();
        }
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


double ARL(uint64_t length, double sigma, double mun, double h)
{
    return (exp(-2.0*mun*(h/sigma+1.166))-1.0+2.0*mun*(h/sigma+1.166))/(2.0*mun*mun)-(double) length;
}

butterworth *initialize_filter(butterworth *lpfilter, uint64_t order, double cutoff, uint64_t length)
{
    if ((lpfilter = malloc(sizeof(butterworth)))==NULL)
    {
        printf("Cannot allocate filter memory\n");
        abort();
    }
    lpfilter->dcof = NULL;
    lpfilter->ccof = NULL;
    lpfilter->temp = NULL;
    lpfilter->tempback = NULL;
    lpfilter->paddedsignal = NULL;
    lpfilter->dcof = dcof_bwlp(order, cutoff);
    lpfilter->ccof = ccof_bwlp(order);
    lpfilter->scale = sf_bwlp(order,cutoff);
    lpfilter->order = order;
    if ((lpfilter->temp = calloc(length+2*order, sizeof(double)))==NULL)
    {
        printf("Cannot allocate temp filter array\n");
        abort();
    }
    if ((lpfilter->tempback = calloc(length+2*order, sizeof(double)))==NULL)
    {
        printf("Cannot allocate tempback array\n");
        abort();
    }
    if ((lpfilter->paddedsignal = calloc(length + 2*order,sizeof(double)))==NULL)
    {
        printf("Cannot allocate padded signal\n");
        abort();
    }
    return lpfilter;
}

void free_filter(butterworth *lpfilter)
{
    free(lpfilter->dcof);
    free(lpfilter->ccof);
    free(lpfilter->temp);
    free(lpfilter->tempback);
    free(lpfilter->paddedsignal);
    free(lpfilter);
}

