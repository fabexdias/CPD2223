#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>

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

char *packnode(tsp_node *node, unsigned int ncities)
{
    char *buffer = malloc((sizeof(tsp_node) - sizeof(unsigned int *) + ncities * sizeof(unsigned int)));
    *(unsigned int *)&buffer[0] = node->length;
    *(unsigned int *)&buffer[4] = node->index;
    *(double *)&buffer[8] = node->cost;
    *(double *)&buffer[16] = node->bound;

    for (size_t i = 0; i < node->length; i++)
    {
        *(unsigned int *)&buffer[24 + i * 4] = node->tour[i];
    }
    return buffer;
}

tsp_node *unpacknode(char *buffer)
{
    unsigned int size = *(unsigned int *)&buffer[0];
    tsp_node *new = tsp_mknode(size);

    new->length = size;
    new->index = *(unsigned int *)&buffer[4];
    new->cost = *(double *)&buffer[8];
    new->bound = *(double *)&buffer[16];
    debug("%u %u %f %f\n", size, new->index, new->cost, new->bound);
    for (size_t i = 0; i < size; i++)
    {
        new->tour[i] = *(unsigned int *)&buffer[24 + i * 4];
    }

    return new;
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

bool array_contains(unsigned int *arr, int size, unsigned int value)
{
    for (int i = 0; i < size; i++)
    {
        if (arr[i] == value)
        {
            return true;
        }
    }
    return false;
}

tsp_result tsp_exe(int rank, int size, tsp_repr rep, double lowerbound, double limit)
{
    double *graph = rep.graph;
    double *short1 = rep.short1;
    double *short2 = rep.short2;
    unsigned int ncities = rep.ncities;

    unsigned int *btour = arrayi_alloc(ncities);
    double btourcost = limit;
    btour[0] = 0;
    tsp_result result;

    MPI_Request sendrequest;
    MPI_Status status;

    int flag = 0;
    int paused = 1;
    int noted = -1, prevnoted = -1;
    char *recvbuff = NULL, *sendbuff = NULL;

    priority_queue_t *queue = queue_create(tsp_queue_cmp);

    tsp_node *current = NULL, *new = NULL, *recv = NULL;
    double cost;
    double update;
    double newBound;

    // Preamble - push all node 0 neighbours to the various process queues
    for (size_t i = rank; i < ncities; i = i + size)
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
            debug("%d) Pushing node %d with bound %f and cost %f to queue\n", rank, (int)i, new->bound, new->cost);
            queue_push(queue, new);
            paused = 0;
        }
    }

    if (rank == 0 && size > 1)
    {
        sendbuff = malloc((size + 1) * sizeof(char) + sizeof(double));

        for (int i = 0; i < size; i++)
        {
            sendbuff[i] = 0;
        }

        sendbuff[size] = paused;

        *(double *)&sendbuff[size + 1] = limit;

        MPI_Send(sendbuff, ((size + 1) * sizeof(char) + sizeof(double)), MPI_CHAR, 1, 1, MPI_COMM_WORLD);
        free(sendbuff);
        sendbuff = NULL;
    }

    int pops = 0;
    srand(time(NULL));
    while (1)
    {
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);

        if (flag == 1 && size > 1)
        {
            flag = 0;
            if (status.MPI_TAG == 1) // token -> recvbuff[0...size-1] = global_paused[0...size-1], recvbuff[size] = to stop (0 not, 1 maybe, 2 stop), recvbuff[size+1...size+9] = bestcost
            {
                recvbuff = malloc((size + 1) * sizeof(char) + sizeof(double));
                MPI_Recv(recvbuff, ((size + 1) * sizeof(char) + sizeof(double)), MPI_CHAR, status.MPI_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                debug("%d) Recieved the token from %d\n", rank, status.MPI_SOURCE);

                if (rank == 0 && recvbuff[size] == 1) // Break all
                    recvbuff[size] = 2;

                if (rank == 0 && paused == 1 && recvbuff[size] == 0) // Start break chain
                    recvbuff[size] = 1;

                if (recvbuff[size] == 1 && paused == 0) // Stop break chain
                {
                    recvbuff[size] = 0;
                }

                if (recvbuff[size] == 2) // Breaking
                {
                    // All threads are paused, we can stop the execution
                    MPI_Send(recvbuff, ((size + 1) * sizeof(char) + sizeof(double)), MPI_CHAR, rank == size - 1 ? 0 : rank + 1, 1, MPI_COMM_WORLD);
                    debug("%d) All threads are paused, we can stop the execution\n", rank);
                    queue_delete(queue);
                    free(queue);
                    free(recvbuff);
                    result.tour = btour;
                    result.cost = btourcost;
                    debug("%d) Returning result with cost %f\n", rank, result.cost);
                    break;
                }

                if (*(double *)&recvbuff[size + 1] < limit)
                    limit = *(double *)&recvbuff[size + 1];

                // SETUP SENDING
                if (btourcost < *(double *)&recvbuff[size + 1])
                    *(double *)&recvbuff[size + 1] = btourcost;

                noted = -1;

                if (queue->size > ncities)
                    for (int i = 0; i < size; i++)
                    {
                        if (recvbuff[i] == 1 && paused == 0 && prevnoted != i && rank != i)
                        {
                            noted = i;
                            recvbuff[i] = 0;
                            break;
                        }
                    }

                recvbuff[rank] = paused;

                MPI_Send(recvbuff, ((size + 1) * sizeof(char) + sizeof(double)), MPI_CHAR, rank == size - 1 ? 0 : rank + 1, 1, MPI_COMM_WORLD);

                free(recvbuff);
                recvbuff = NULL;
            }
            else if (status.MPI_TAG == 2)
            {
                paused = 0;
                recvbuff = malloc((sizeof(tsp_node) - sizeof(unsigned int *) + ncities * sizeof(unsigned int)));
                debug("%d) Packdge recieved!\n", rank);
                MPI_Recv(recvbuff, (sizeof(tsp_node) - sizeof(unsigned int *) + (ncities) * sizeof(unsigned int)), MPI_CHAR, status.MPI_SOURCE, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                recv = unpacknode(recvbuff);
                debug("%d) Packdge recieved: %d %f\n", rank, recv->length, recv->bound);
                free(recvbuff);
                recvbuff = NULL;
                queue_push(queue, recv);
            }
        }

        if (queue->size > 0)
        {
            current = queue_pop(queue);
            pops++;
            // debug("%d) Popping.\n", rank);

            if (current->bound > btourcost || current->bound >= limit)
            {
                // If the best node doesn't work, then the others in the queue probably don't work either
                // Delete all the nodes because they are 100% garbage
                debug("%d) Clearing the queue.\n", rank);
                while (queue->size > 0)
                {
                    // NOTE: There's some optimization potential here because the bubbe-down step is unnecessary
                    tsp_delnode(current);
                    current = queue_pop(queue);
                }
            }
            else if (current->length == ncities)
            {
                double newcost = current->cost + matrix_read(graph, ncities, current->index, 0);
                if (newcost < btourcost && newcost < limit)
                {
                    for (size_t i = 1; i < ncities; i++)
                    {
                        btour[i] = current->tour[i];
                    }
                    btourcost = newcost;
                }
            }
            else if ((pops > 20000 && size > 1 && size < 16) || (pops > 7500 && size > 15) || ((noted > -1) && prevnoted != noted))
            {
                int index = 0;

                if (noted > -1)
                    index = noted;
                else
                    index = rand() % size;

                if (index == rank)
                {
                    if (rank == 0)
                        index = size - 1;
                    else
                        index--;
                }

                debug("%d) Sending pops to %d %d, %f, (%f)\n", rank, index, current->length, current->bound, limit);
                sendbuff = packnode(current, ncities);
                MPI_Isend(sendbuff, (sizeof(tsp_node) - sizeof(unsigned int *) + (ncities) * sizeof(unsigned int)), MPI_CHAR, index, 2, MPI_COMM_WORLD, &sendrequest);

                if (noted > -1)
                    prevnoted = noted;
                else
                    pops = 0;
            }
            else
            {
                for (size_t i = 0; i < ncities; i++)
                {
                    cost = matrix_read(graph, ncities, current->index, i);
                    if (cost != INFINITY && i != current->index && !array_contains(current->tour, current->length, (unsigned int)i))
                    {
                        update = (cost >= short2[i] ? short2[i] : short1[i]) + (cost >= short2[current->index] ? short2[current->index] : short1[current->index]);
                        newBound = current->bound + cost - (update / 2);
                        if (newBound > btourcost || newBound > limit)
                        {
                            continue;
                        }
                        new = tsp_mknode(current->length + 1);
                        for (size_t j = 0; j < current->length; j++)
                        {
                            new->tour[j] = current->tour[j];
                        }
                        new->tour[current->length] = i;
                        new->cost = current->cost + cost;
                        new->bound = newBound;
                        new->length = current->length + 1;
                        new->index = i;

                        // Distribute the new node to another random process if own queue still has nodes and other process is paused
                        queue_push(queue, new);
                    }
                }
            }
            tsp_delnode(current);
            if ((noted > -1 && prevnoted == noted) || (pops == 0 && size > 1))
            {
                debug("%d) Waiting for send.\n", rank);
                MPI_Wait(&sendrequest, MPI_STATUSES_IGNORE);
                free(sendbuff);
                sendbuff = NULL;
                noted = -1;
            }

            if (queue->size == 0)
            {
                debug("%d) Queue ces't finni!\n", rank);
                paused = 1;
                if (size == 1)
                {
                    queue_delete(queue);
                    free(queue);
                    result.tour = btour;
                    result.cost = btourcost;
                    debug("%d) Returning result with cost %f\n", rank, result.cost);
                    break;
                }
            }
        }
    }
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

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Argument validation
    if (argc <= 2 || argc > 3)
    {
        if (rank == 0)
        {
            error("Too many arguments.\n");
            help(argv[0]);
        }
        MPI_Finalize();
        return 0;
    }

    // Make sure that the lowerbound arg is a number
    if (argc != 3 || atof(argv[2]) <= 0)
    {
        // Either the value is 0 (not allowed) or it is not a number.
        // If that's the case ignore this.
        if (rank == 0)
        {
            error("The lowerbound argument doesn't seem to be a positive number. Ignoring.\n");
        }
        MPI_Finalize();
        return 0;
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
        if (rank == 0)
            syserr("Could not open the input file");
        MPI_Finalize();
        return 0;
    }

    tsp_repr t = tsp_mkrepr(input);
    if (!t.valid)
    {
        // An earlier error happened preventing us from carrying on
        tsp_delrepr(t);
        MPI_Finalize();
        return 0;
    }

    double lowerbound = 0;
    for (size_t i = 0; i < t.ncities; ++i)
    {
        lowerbound += t.short1[i] + t.short2[i];
    }
    lowerbound /= 2;
    info("Lowerbound at root = %f\n", lowerbound);
    double *overallbest = malloc(size * sizeof(double));

    exec_time = -MPI_Wtime();

    tsp_result result = tsp_exe(rank, size, t, lowerbound, limit);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allgather(&result.cost, 1, MPI_DOUBLE, overallbest, 1, MPI_DOUBLE, MPI_COMM_WORLD);

    exec_time += MPI_Wtime();

    if (rank == 0)
        fprintf(stderr, "%.1fs\n", exec_time);

    int index = -1;
    double min = limit;
    for (int i = 0; i < size; i++)
    {
        if (overallbest[i] < min)
        {
            index = i;
            min = overallbest[i];
        }
    }
    free(overallbest);
    // printf("%f %f %d %d\n", overallbest, result.cost, result.tour[1], result.tour[t.ncities - 1]);

    if (index == rank)
    {
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
    }

    // Cleanup
    free(result.tour);
    tsp_delrepr(t);
    MPI_Finalize();
    return 0;
}
