#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define INF 1000000 // An arbitrarily large number for infinity

typedef struct
{
    int id;
    int parent_id;
    int level;
    int *path;
    int path_len;
    int bound;
} node_t;

// Compare function for the priority queue
int cmp_node(const void *a, const void *b)
{
    node_t *na = (node_t *)a;
    node_t *nb = (node_t *)b;
    return na->bound - nb->bound;
}

int *create_graph(int n)
{
    // Allocate memory for a complete graph with n nodes
    int *graph = (int *)malloc(sizeof(int) * n * n);

    // Fill in the graph with random distances between nodes
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j)
            {
                graph[i * n + j] = 0;
            }
            else
            {
                graph[i * n + j] = rand() % 100;
            }
        }
    }

    return graph;
}

int calc_bound(node_t *node, int *graph, int n)
{
    // Calculate the bound for a given node

    // Add up the distance so far
    int bound = 0;
    for (int i = 0; i < node->path_len - 1; i++)
    {
        bound += graph[node->path[i] * n + node->path[i + 1]];
    }

    // Add in the minimum outgoing edge from each node not in the path
    for (int i = 0; i < n; i++)
    {
        if (i == node->id || i == node->path[0] || i == node->path[node->path_len - 1])
        {
            continue;
        }

        int min_edge = INF;
        for (int j = 0; j < n; j++)
        {
            if (j == i || j == node->id || j == node->path[0] || j == node->path[node->path_len - 1])
            {
                continue;
            }
            if (graph[i * n + j] < min_edge)
            {
                min_edge = graph[i * n + j];
            }
        }
        bound += min_edge;
    }

    return bound;
}

void expand_node(node_t *node, node_t *children, int *graph, int n, int *num_children)
{
    // Generate all children of a given node

    *num_children = 0;

    for (int i = 0; i < n; i++)
    {
        // Check if i is already in the path
        int in_path = 0;
        for (int j = 0; j < node->path_len; j++)
        {
            if (i == node->path[j])
            {
                in_path = 1;
                break;
            }
        }

        if (!in_path)
        {
            node_t child;
            child.id = i;
            child.parent_id = node->id;
            child.level = node->level + 1;
            child.path_len = node->path_len + 1;
            child.path = (int *)malloc(sizeof(int) * child.path_len);
            for (int j = 0; j < node->path_len; j++)
            {
                child.path[j] = node->path[j];
            }
            child.path[node->path_len] = i;
            // Calculate the bound for the child
            child.bound = calc_bound(&child, graph, n);

            children[*num_children] = child;
            (*num_children)++;
        }
    }
}

void TSP(int rank, int size, int *graph, int n)
{
    int num_nodes_expanded = 0;
    // Initialize the priority queue with the root node
    node_t root;
    root.id = 0;
    root.parent_id = -1;
    root.level = 0;
    root.path_len = 1;
    root.path = (int *)malloc(sizeof(int));
    root.path[0] = 0;
    root.bound = calc_bound(&root, graph, n);

    node_t *queue = (node_t *)malloc(sizeof(node_t) * size);
    int queue_size = 1;
    queue[0] = root;

    int *min_tour = (int *)malloc(sizeof(int) * n);
    int min_tour_len = INF;
    int global_min_tour_len;

    while (queue_size > 0)
    {
        // Get the next node from the queue
        node_t node = queue[0];
        for (int i = 1; i < queue_size; i++)
        {
            queue[i - 1] = queue[i];
        }
        queue_size--;

        if (node.level == n - 1)
        {
            // We've reached a leaf node, so check if it forms a better tour
            int tour_len = node.bound + graph[node.path[node.path_len - 1] * n];
            if (tour_len < min_tour_len)
            {
                min_tour_len = tour_len;
                for (int i = 0; i < n; i++)
                {
                    min_tour[i] = node.path[i];
                }
            }
        }
        else
        {
            // Expand the node
            node_t children[n];
            int num_children;
            expand_node(&node, children, graph, n, &num_children);
            num_nodes_expanded++;

            // Add the children to the queue
            for (int i = 0; i < num_children; i++)
            {
                // Check if this child can possibly produce a better tour
                if (children[i].bound < min_tour_len)
                {
                    // Add the child to the queue
                    queue[queue_size] = children[i];
                    queue_size++;
                }
            }

            // Sort the queue by bound
            qsort(queue, queue_size, sizeof(node_t), cmp_node);
        }

        // Check if we need to send or receive nodes
        if (queue_size < size && rank == 0)
        {
            // Send a node to the first idle process
            int dest_rank = -1;
            for (int i = 1; i < size; i++)
            {
                int is_idle;
                MPI_Recv(&is_idle, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if (is_idle)
                {
                    dest_rank = i;
                    break;
                }
            }
            if (dest_rank != -1)
            {
                MPI_Send(&queue[queue_size], sizeof(node_t), MPI_BYTE, dest_rank, 0, MPI_COMM_WORLD);
                queue_size++;
            }
        }
        else if (queue_size > 0 && rank != 0)
        {
            // Receive a node from process 0
            int is_idle = 0;
            MPI_Send(&is_idle, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            node_t child;
            MPI_Recv(&child, sizeof(node_t), MPI_BYTE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Expand the received node
            node_t children[n];
            int num_children;
            expand_node(&child, children, graph, n, &num_children);
            num_nodes_expanded++;

            // Add the children to the queue
            for (int i = 0; i < num_children; i++)
            {
                // Check if this child can possibly produce a better tour
                if (children[i].bound < min_tour_len)
                {
                    // Add the child to the queue
                    queue[queue_size] = children[i];
                    queue_size++;
                }
            }

            // Sort the queue by bound
            qsort(queue, queue_size, sizeof(node_t), cmp_node);
        }

        // Check if we need to broadcast the current min tour length
        if (queue_size == 0)
        {
            MPI_Allreduce(&min_tour_len, &global_min_tour_len, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
            if (global_min_tour_len < min_tour_len)
            {
                min_tour_len = global_min_tour_len;
            }
            else
            {
                break;
            }
        }
    }

    // Print the results
    if (rank == 0)
    {
        printf("Min tour length: %d\n", min_tour_len);
        printf("Min tour: ");
        for (int i = 0; i < n; i++)
        {
            printf("%d ", min_tour[i]);
        }
        printf("\n");
        printf("Nodes expanded: %d\n", num_nodes_expanded);
    }
}
int main(int argc, char **argv)
{
    int rank, size;
    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Read the input file
    char *input_filename = argv[1];
    int n;
    int *graph = read_input(input_filename, &n);

    // Run the TSP algorithm
    TSP(rank, size, graph, n);

    // Finalize MPI
    MPI_Finalize();

    return 0;
}
