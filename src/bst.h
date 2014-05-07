#ifndef BST_H
#define BST_H

#include <stdint.h>

void* bst_alloc( size_t n );

/**
 * @brief bst_compute
 * @param root
 * @param p
 * @param q
 * @param n
 * @return
 */
double bst_compute( int* root, double* p, double* q, size_t n, void*_mem );

double bst_root_get( double* root, size_t i, size_t j );

void bst_free( void* _mem );

#endif // BST_H
