#include<bst.h>
#include<math.h>
#include<string.h>


/*
 * Disclaimer: The indices differ from the high-level pseudo-code description
 * found in the book [1].
 *
 * Here, e[i,j] is the expected value for the tree containing keys i to j-1,
 * since the first key is p_0 and not p_1 as in the book. w[.,.] and root[.,.]
 * are defined similarly.
 *
 */

/*
 * This implementation differs from the reference by traversing the e table
 * row by row bottom up rather than diagonal by diagonal.
 */

#define STRIDE (n+1)
#define IDX(i,j) ((i)*STRIDE + j)

typedef struct {
    double* e;
    double* w;
    int* r;
    size_t n;
} segments_t;

void* bst_alloc_100_ref_bottomup( size_t n ) {
    segments_t* mem = (segments_t*) malloc( sizeof(segments_t) );
    size_t sz = (n+1)*(n+1);
    mem->e = (double*) malloc( sz * sizeof(double) );
    mem->w = (double*) malloc( sz * sizeof(double) );
    mem->r = (int*)    malloc( sz * sizeof(int) );
    memset( mem->r, -1, sz * sizeof(int) );
    return mem;
}

#include <stdio.h>
double bst_compute_100_ref_bottomup( void*_bst_obj, double* p, double* q, size_t n ) {
    segments_t* mem = (segments_t*) _bst_obj;
    int i, l, r, j;
    double t;
    double* e = mem->e, *w = mem->w;
    int* root = mem->r;
    // initialization
    mem->n = n;
    for( i = 0; i < n+1; i++ ) {
        e[IDX(i,i)] = q[i];
        w[IDX(i,i)] = q[i];
    }

    for (i = n-1; i >= 0; --i) {
        for (j = i+1; j < n+1; ++j) {
            e[IDX(i,j)] = INFINITY;
            w[IDX(i,j)] = w[IDX(i,j-1)] + p[j-1] + q[j];
            for (r=i; r<j; ++r) {
                t = e[IDX(i,r)] + e[IDX(r+1,j)] + w[IDX(i,j)];
                if (t < e[IDX(i,j)]) {
                    e[IDX(i,j)] = t;
                    root[IDX(i,j)] = r;
                }
            }
        }
    }

    return e[IDX(0,n)];
}

size_t bst_get_root_100_ref_bottomup( void* _bst_obj, size_t i, size_t j )
{
    // [i,j], in table: [i-1, j]+1
    segments_t *mem = _bst_obj;
    size_t n = mem->n;
    int *root = mem->r;
    return (size_t) root[(i-1)*(n+1)+j]+1;
}

void bst_free_100_ref_bottomup( void* _mem ) {
    segments_t* mem = (segments_t*) _mem;
    free( mem->e );
    free( mem->w );
    free( mem->r );
    free( mem );
}
