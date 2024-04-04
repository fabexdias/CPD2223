#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include "debug.h"
#include "matrix.h"

#include "lib/nqueue/queue.h"

#define DELTA 4

typedef struct
{
    bool valid;
    unsigned int ncities;
    double *graph;
    double *short1;
    double *short2;
} tsp_repr;

typedef struct
{
    unsigned int *tour;
    double cost;
    double bound;
    unsigned int length;
    unsigned int index;
} tsp_node;

typedef struct
{
    unsigned int *tour;
    double cost;
} tsp_result;

void tsp_delrepr(tsp_repr t)
{
    if (t.graph)
    {
        free(t.graph);
    }
    if (t.short1)
    {
        free(t.short1);
    }
    if (t.short2)
    {
        free(t.short2);
    }
}

tsp_repr tsp_mkrepr(FILE *input)
{
    tsp_repr t;
    t.graph = NULL;
    t.short1 = NULL;
    t.short2 = NULL;

    int nroutes = -1;
    int from = -1, to = -1;
    double cost = -1;

    if (fscanf(input, "%u %d", &(t.ncities), &nroutes) == 2)
    {
        t.graph = matrix_alloc(t.ncities);
        t.short1 = array_alloc(t.ncities);
        t.short2 = array_alloc(t.ncities);

        if (!t.graph || !t.short1 || !t.short2)
        {
            error("Failed to allocate memory; aborting!\n");
            t.valid = false;
            return t;
        }
    }

    for (int i = 0; i < nroutes; i++)
    {
        if (fscanf(input, "%d %d %lf", &from, &to, &cost) == 3)
        {
            matrix_write(t.graph, t.ncities, from, to, cost);
            matrix_write(t.graph, t.ncities, to, from, cost);

            if (t.short1[from] > cost)
            {
                t.short2[from] = t.short1[from];
                t.short1[from] = cost;
            }
            else if (t.short2[from] > cost)
            {
                t.short2[from] = cost;
            }

            if (t.short1[to] > cost)
            {
                t.short2[to] = t.short1[to];
                t.short1[to] = cost;
            }
            else if (t.short2[to] > cost)
            {
                t.short2[to] = cost;
            }
        }
        else
        {
            error("Error reading file.\n");
            t.valid = false;
            return t;
        }
    }

    fclose(input);

    t.valid = true;
    return t;
}

tsp_node *tsp_mknode(unsigned int length)
{
    tsp_node *node = malloc(sizeof(tsp_node));
    node->tour = arrayi_alloc(length);
    return node;
}

void tsp_delnode(tsp_node *node)
{
    if (node->tour)
    {
        free(node->tour);
    }

    free(node);
}

char tsp_queue_cmp(void *a, void *b)
{
    // Lowest lower-bound goes first; if both happen to be tied, the one with the lowest index goes first.
    if (((tsp_node *)a)->bound == ((tsp_node *)b)->bound)
    {
        return (((tsp_node *)a)->index > ((tsp_node *)b)->index);
    }

    return (((tsp_node *)a)->bound > ((tsp_node *)b)->bound);
}

tsp_result tsp_exe(tsp_repr rep, double lowerbound, double limit)
{
    double *graph = rep.graph;
    double *short1 = rep.short1;
    double *short2 = rep.short2;
    unsigned int ncities = rep.ncities;

    unsigned int *btour = arrayi_alloc(ncities);
    double btourcost = INFINITY;

    tsp_result result;

    priority_queue_t *queue = queue_create(tsp_queue_cmp);
    tsp_node *current = tsp_mknode(ncities);

    btour[0] = 0;
    current->tour[0] = 0;
    current->cost = 0;
    current->bound = lowerbound;
    current->length = 1;
    current->index = 0;

    queue_push(queue, current);

    while (queue->size)
    {
        current = queue_pop(queue);

        if (current->bound >= btourcost)
        {
            tsp_delnode(current);
            break;
        }

        if (current->length == ncities)
        {
            if (current->cost + matrix_read(graph, ncities, current->index, 0) < btourcost)
            {
                btourcost = current->cost + matrix_read(graph, ncities, current->index, 0);
                for (unsigned int i = 0; i < ncities; i++)
                {
                    btour[i] = current->tour[i];
                }
            }
        }
        else
        {
            debug("Level: %u\n", current->length);
            for (unsigned int i = 0; i < ncities; i++)
            {
                bool ontour = false;
                double cost = matrix_read(graph, ncities, current->index, i);
                if (cost != INFINITY && current->index != i)
                {
                    for (unsigned int j = 0; j < current->length; j++)
                    {
                        if (current->tour[j] == i)
                        {
                            ontour = true;
                            debug("Visiting node %u\n", current->index);
                            break;
                        }
                    }
                    if (!ontour)
                    {
                        double update = (cost >= short2[i] ? short2[i] : short1[i]) + (cost >= short2[current->index] ? short2[current->index] : short1[current->index]);

                        double newBound = current->bound + cost - (update / 2);
                        if (newBound > btourcost || newBound > limit)
                        {
                            continue;
                        }

                        tsp_node *new = tsp_mknode(current->length + 1);
                        for (unsigned int j = 0; j < current->length; j++)
                        {
                            new->tour[j] = current->tour[j];
                        }
                        new->tour[current->length] = i;
                        new->cost = current->cost + cost;
                        new->bound = newBound;
                        new->length = current->length + 1;
                        new->index = i;
                        queue_push(queue, new);
                    }
                }
            }
        }
        tsp_delnode(current);
        debug("Queue size: %lu\n", queue->size);
    }

    result.tour = btour;
    result.cost = btourcost;

    while (queue->size > 0)
    {
        current = queue_pop(queue);
        tsp_delnode(current);
    }
    queue_delete(queue);
    free(queue);
    return result;
}

void help(char *me)
{
    printf("USAGE: %s inputfile lowerbound\n * Where inputfile is a file;\n * Where lowerbound is a number.\n", me);
}

int main(int argc, char *argv[])
{
    double limit = INFINITY;
    double exec_time;

    // Argument validation
    if (argc <= 2)
    {
        error("No arguments provided.\n");
        help(argv[0]);
        return 1;
    }
    else if (argc > 3)
    {
        error("Too many arguments.\n");
        help(argv[0]);
        return 1;
    }

    // Make sure that the lowerbound arg is a number
    if (argc != 3 || atof(argv[2]) <= 0)
    {
        // Either the value is 0 (not allowed) or it is not a number.
        // If that's the case ignore this.
        warn("The lowerbound argument doesn't seem to be a positive number. Ignoring.\n");
    }
    else
    {
        limit = (double)atof(argv[2]);
    }
    info("Cost must be <= %f\n", limit);
    info("File target: %s\n", argv[1]);
    FILE *input = fopen(argv[1], "r");
    if (!input)
    {
        syserr("Could not open the input file");
        return 1;
    }

    tsp_repr t = tsp_mkrepr(input);
    if (!t.valid)
    {
        // An earlier error happened preventing us from carrying on
        tsp_delrepr(t);
        return 1;
    }

    double lowerbound = 0;
    for (size_t i = 0; i < t.ncities; ++i)
    {
        lowerbound += t.short1[i] + t.short2[i];
    }
    lowerbound /= 2;
    info("Lowerbound at root = %f\n", lowerbound);

    exec_time = -omp_get_wtime();

    tsp_result result = tsp_exe(t, lowerbound, limit);

    exec_time += omp_get_wtime();

    fprintf(stderr, "%.1fs\n", exec_time);
    if (lowerbound > limit)
    {
        info("Lowerbound %f is higher than the desired limit %f.\n", lowerbound, limit);
        printf("NO SOLUTION\n");
    }
    else if (result.cost > limit)
    {
        info("Either the graph is not connected OR the lowerbound limit is too low!\n");
        printf("NO SOLUTION\n");
    }
    else
    {
        // Print the tour
        printf("%.1f\n", result.cost);
        printf("%d", result.tour[0]);
        for (size_t i = 1; i < t.ncities; ++i)
        {
            printf(" %d", result.tour[i]);
        }
        printf(" 0\n");
    }

    // Cleanup
    free(result.tour);
    tsp_delrepr(t);
    return 0;
}
