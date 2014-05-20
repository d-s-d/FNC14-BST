#ifndef BST_H
#define BST_H

#include <stdlib.h>

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

size_t bst_get_root( void* _bst_obj, size_t i, size_t j );

void bst_free( void* _bst_obj );

size_t bst_flops( size_t n );

// Function types

typedef void*  (*bst_alloc_fn)( size_t n );
typedef double (*bst_compute_fn)( void *_bst_obj, double* p, double* q, size_t n );
typedef size_t (*bst_get_root_fn)( void* _bst_obj, size_t i, size_t j );
typedef void   (*bst_free_fn)( void* _bst_obj );
typedef size_t (*bst_flops_fn) ( size_t n );

typedef struct {
    char            *name;
    bst_alloc_fn    alloc;
    bst_compute_fn  compute;
    bst_get_root_fn root;
    bst_free_fn     free;
    bst_flops_fn    flops;
} bst_impl_t;

#endif // BST_H
