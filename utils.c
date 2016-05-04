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
    if ((head = malloc(sizeof(edge)))==NULL)
    {
        printf("Cannot allocate head node\n");
        exit(28);
    }
    head->next = NULL;
    head->location = 0;
    head->type = HEAD;
    return head;
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
            exit(29);
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
        exit(30);
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
    return head;
}

event *add_event(event *current, uint64_t start, uint64_t finish, uint64_t index)
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
    if ((head=malloc(sizeof(cusumlevel)))==NULL)
    {
        printf("Cannot allocate head level node\n");
        exit(32);
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
            exit(33);
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



