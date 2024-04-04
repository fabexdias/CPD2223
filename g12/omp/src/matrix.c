#include "matrix.h"
#include <limits.h>
#include <math.h>
#include <stdlib.h>

unsigned int* arrayi_alloc(const size_t size) {
    unsigned int* m = malloc(sizeof(unsigned int) * size);
    for (size_t s = 0; s < size; ++s) {
        m[s] = UINT_MAX;
    }

    return m;
}

double* array_alloc(const size_t size) {
    double* m = malloc(sizeof(double) * size);
    for (size_t s = 0; s < size; ++s) {
        m[s] = INFINITY;
    }

    return m;
}

double* matrix_alloc(const size_t size) {
    double* m = malloc(sizeof(double) * size * size);
    for (size_t s = 0; s < size * size; ++s) {
        m[s] = INFINITY;
    }

    return m;
}

inline double matrix_read(double* matrix, const size_t size, const size_t x, const size_t y) {
    return matrix[x * size + y];
}

inline void matrix_write(double* matrix, const size_t size, const size_t x, const size_t y, double val) {
    matrix[x * size + y] = val;
}
