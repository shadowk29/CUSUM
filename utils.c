#include"utils.h"
#include<stdio.h>
#include<stdlib.h>


int signum(double num)
{
    return (EPS<num)-(num<-EPS);
}
double my_min(double a, double b)
{
    return a < b ? a : b;
}

double my_max(double a, double b)
{
    return a > b ? a : b;
}

double d_abs(double num)
{
    return num > 0 ? num : -num;
}

int64_t intmin(int64_t a, int64_t b)
{
    return a < b ? a : b;
}
int64_t intmax(int64_t a, int64_t b)
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
    if (current)
    {
        edge *temp;
        while (current->next)
        {
            temp = current;
            current = current->next;
            free(temp);
        }
        free(current);
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
    head->index = HEAD;
    head->first_edge = NULL;
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
        current->first_edge = NULL;
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
        current->next->next = NULL;
        current->next->signal = NULL;
        current->next->prev = current;
        current->next->index = current->index + 1;
        current->next->start = start;
        current->next->finish = finish;
        current->next->length = finish-start;
        current = current->next;
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
        free(current->signal);
        free(current->filtered_signal);
        if (current->first_edge)
        {
            free(current->first_edge);
        }
        free(current);
    }
}
