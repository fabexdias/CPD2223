#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include "debug.h"
#include "matrix.h"

#include "lib/nqueue/queue.h"

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

#define LOCK_QUEUE(i) omp_set_lock(locks + i)
#define UNLOCK_QUEUE(i) omp_unset_lock(locks + i)

tsp_result tsp_exe(tsp_repr rep, double lowerbound, double limit)
{
    double *graph = rep.graph;
    double *short1 = rep.short1;
    double *short2 = rep.short2;
    unsigned int ncities = rep.ncities;
    const unsigned int thread_num = omp_get_max_threads();
    info("Running with numthreads = %d\n", thread_num);

    unsigned int *btour = arrayi_alloc(ncities);
    double btourcost = limit;
    btour[0] = 0;
    tsp_result result;

    /**
        Each thread will have it's own priority queue
        To ensure orderly and relatively fine-grained load balancing, each queue will have it's own lock
    */
    priority_queue_t **queues = malloc(thread_num * sizeof(priority_queue_t *));
    omp_lock_t *locks = malloc((thread_num) * sizeof(omp_lock_t));
    bool *waiting = malloc(thread_num * sizeof(bool));

    unsigned int finish = 0;
    for (size_t k = 0; k < thread_num; k++)
    {
        queues[k] = queue_create(tsp_queue_cmp);
        waiting[k] = false;
        omp_init_lock(locks + k); // omp lock functions use pointers, hence why we're doing this
    }

    info("Starting parallel\n");
#pragma omp parallel default(none) \
    shared(stderr, queues, locks, btour, btourcost, graph, short1, short2, ncities, thread_num, waiting, finish, lowerbound)
    {
        int idx = omp_get_thread_num();
        tsp_node *current = NULL, *new = NULL;
        debug("Running preamble, %d\n", idx);
        LOCK_QUEUE(idx);
        double cost;
        double update;
        double newBound;

#pragma omp for nowait
        for (size_t i = 0; i < ncities; i++)
        {
            cost = matrix_read(graph, ncities, 0, i);
            if (cost != INFINITY && 0 != i)
            {
                update = (cost >= short2[i] ? short2[i] : short1[i]) + (cost >= short2[0] ? short2[0] : short1[0]);
                newBound = lowerbound + cost - (update / 2);
                if (newBound > btourcost)
                {
                    continue;
                }
                new = tsp_mknode(2);
                new->tour[0] = 0;
                new->tour[1] = i;
                new->cost = cost;
                new->bound = newBound;
                new->length = 2;
                new->index = i;
                debug("Pushing for thread = %d\n", idx);
                queue_push(queues[idx], new);
            }
        }
        UNLOCK_QUEUE(idx);

        do
        {
            // Each thread pops the top element of the queue (if available)
            if (queues[idx]->size > 0)
            {
                debug("Popping!\n");
                LOCK_QUEUE(idx);
                current = queue_pop(queues[idx]);
                // There needs to be a counter that is reset every iteration, if all threads enter this if exiting needs to be activated
                if (current->bound > btourcost)
                {
                    // If the best node doesn't work, then the others in the queue probably don't work either
                    // Delete all the nodes because they are 100% garbage
                    debug("Only low quality nodes at thread %d. Clearing the queue.\n", idx);
                    while (queues[idx]->size > 0)
                    {
                        // NOTE: There's some optimization potential here because the bubbe-down step is unnecessary
                        tsp_delnode(current);
                        current = queue_pop(queues[idx]);
                    }
#pragma omp atomic write
                    waiting[idx] = true;
#pragma omp atomic update
                    finish++;
                    UNLOCK_QUEUE(idx);
                }
                else if (current->length == ncities)
                {
                    // This node represents a complete loop, so a tour cost can be computed
                    // It becomes the new solution if it's better than the solution computed so far
                    UNLOCK_QUEUE(idx);
#pragma omp critical(update_btour)
                    {
                        double newcost = current->cost + matrix_read(graph, ncities, current->index, 0);
                        if (newcost < btourcost)
                        {
                            for (size_t i = 1; i < ncities; i++)
                            {
                                btour[i] = current->tour[i];
                            }

                            btourcost = newcost;
                        }
                        else if (newcost == btourcost)
                        {
                            // If the cost is equal, prefer the one going through lower-numbered nodes first.
                            bool check = true;
                            for (size_t i = 1; i < ncities; i++)
                            {
                                if (check && btour[i] > current->tour[i])
                                {
                                    break;
                                }
                                else if (check && btour[i] < current->tour[i])
                                {
                                    check = false;
                                }
                                btour[i] = current->tour[i];
                            }
                        }
                    }
                }
                else
                {
                    // Visit this node and generate all children nodes for it
                    UNLOCK_QUEUE(idx);
                    debug("Level: %u\n", current->length);
                    debug("READ! %p: length %d tid %d\n", (void *)current->tour, current->length, idx);
                    for (size_t c = 0; c < ncities; c++)
                    {
                        bool alreadyAdded = false;
                        bool pushed = true;
                        size_t here = current->index;
                        double cost = matrix_read(graph, ncities, here, c);
                        if (cost != INFINITY && c != here)
                        {

                            // The path to the city exists AND it's not to the same city we were already
                            for (unsigned int j = 0; j < current->length; j++)
                            {
                                debug("Checking a given tour.\n");
                                if (current->tour[j] == c)
                                {
                                    // We've already gone to this city, though
                                    alreadyAdded = true;
                                    break;
                                }
                            }
                            if (alreadyAdded)
                            {
                                // We'd be visiting a node we've already visited, so ignore this node
                                continue;
                            }

                            update = (cost >= short2[c] ? short2[c] : short1[c]) + (cost >= short2[here] ? short2[here] : short1[here]);
                            newBound = current->bound + cost - (update / 2);

                            // Make sure that this path is decent enough. Otherwise skip
                            if (newBound > btourcost)
                            {
                                continue;
                            }

                            // This node is good!
                            new = tsp_mknode(current->length + 1);
                            for (unsigned int j = 0; j < current->length; j++)
                            {
                                new->tour[j] = current->tour[j];
                            }
                            new->tour[current->length] = c;
                            new->cost = current->cost + cost;
                            new->bound = newBound;
                            new->length = current->length + 1;
                            new->index = c;

                            if (!pushed || finish == 0)
                            {
                                LOCK_QUEUE(idx);
                                queue_push(queues[idx], new);
                                UNLOCK_QUEUE(idx);
                                pushed = true;
                            }
                            else
                            {
                                for (size_t ii = 0; ii < thread_num; ii++)
                                {
                                    if (waiting[ii]) // send work to waiting threads
                                    {
                                        LOCK_QUEUE(ii);
                                        if (!waiting[ii]) // send work to waiting threads
                                        {
                                            UNLOCK_QUEUE(ii);
                                            continue;
                                        }
                                        debug("ADD! %p: length %d tid %d\n", (void *)new->tour, new->length, idx);
                                        queue_push(queues[ii], new);
#pragma omp atomic update
                                        finish--;
                                        waiting[ii] = false;

                                        UNLOCK_QUEUE(ii);
                                        pushed = false;
                                        break;
                                    }
                                    else if (finish == 0)
                                    {
                                        break;
                                    }
                                }

                                if (pushed)
                                {
                                    LOCK_QUEUE(idx);
                                    queue_push(queues[idx], new);
                                    UNLOCK_QUEUE(idx);
                                    pushed = false;
                                }
                            }
                            debug("Done pushing!\n");
                        }
                    }
                }
                if (current != NULL)
                {
                    tsp_delnode(current);
                }
                debug("Queue size: %lu\n", queues[idx]->size);
            }
            else
            {
                waiting[idx] = true;
#pragma omp atomic update
                finish++;
                debug("Queue is empty\n");
            }

            if (finish == thread_num)
            {
                debug("%u %u\n", finish, thread_num);
                for (size_t k = 0; k < thread_num; k++)
                {
                    waiting[k] = false;
                    debug("%u\n", k);
                }
            }

            int i;
            while (waiting[idx])
            {
#pragma omp atomic update
                i++;
            }
        } while (finish != thread_num);
    }

    result.tour = btour;
    result.cost = btourcost;

    for (size_t k = 0; k < thread_num; k++)
    {
        while (queues[k]->size > 0)
        {
            tsp_delnode(queue_pop(queues[k]));
        }
        queue_delete(queues[k]);
        free(queues[k]);
        omp_destroy_lock(locks + k);
    }
    free(queues);
    free(waiting);

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
