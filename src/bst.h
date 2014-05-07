#ifndef BST_H
#define BST_H

#include <stdint.h>

void* bst_alloc( size_t n );

/**
 * @brief bst_compute
 * @param _bst_obj an internal data structure
 * @param p array containing probabilities of the key (lengh n)
 * @param q array containing probabilities of the dummy keys (length n+1)
 * @param n size of p
 * @return expected propability of the entire tree
 */
double bst_compute( void *_bst_obj, double* p, double* q, size_t n );

int bst_get_root( void* _bst_obj, size_t i, size_t j );

void bst_free( void* _mem );

#endif // BST_H
