/*
    Matrix heap data structure. Optimized for as least memory accesses as possible
*/

#pragma once
#include <stdlib.h>
#include <string.h>

unsigned int* arrayi_alloc(const size_t size);
double* array_alloc(const size_t size);
double* matrix_alloc(const size_t size);

extern double matrix_read(double* matrix, const size_t size, const size_t x, const size_t y);
extern void matrix_write(double* matrix, const size_t size, const size_t x, const size_t y, double val);
