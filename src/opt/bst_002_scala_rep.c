#include<bst.h>
#include<math.h>
#include<string.h>

//#define DIAG_NUM(j,i) (j-i+1)
//#define EULER_SUM(n,d) ((n-d+1)*(n-d+2)/2)
//#define IDX(n,j,i) (EULER_SUM(n,0) - EULER_SUM(n,DIAG_NUM(j,i)) - 1 + i)

/*
 * Disclaimer: The indices differ from the high-level pseudo-code description
 * found in the book [1].
 *
 * Here, e[i,j] is the expected value for the tree containing keys i to j-1,
 * since the first key is p_0 and not p_1 as in the book. w[.,.] and root[.,.]
 * are defined similarly.
 *
 */

#define STRIDE (n+1)
#define IDX(i,j) ((i)*STRIDE + j)

typedef struct {
    double* e;
    double* w;
    int* r;
    size_t n;
} segments_t;

void* bst_alloc_002_scala_rep( size_t n ) {
    segments_t* mem = (segments_t*) malloc( sizeof(segments_t) );
    size_t sz = (n+1)*(n+1);
    mem->e = (double*) malloc( sz * sizeof(double) );
    mem->w = (double*) malloc( sz * sizeof(double) );
    mem->r = (int*)    malloc( sz * sizeof(int) );
    memset( mem->r, -1, sz * sizeof(int) );
    return mem;
}

double bst_compute_002_scala_rep( void*_bst_obj, double* p, double* q,
                                  size_t n ) {
    segments_t* mem = (segments_t*) _bst_obj;
    size_t i, l, r, j;
    double t, t_min, w_cur;
    int r_min;
    double* e = mem->e, *w = mem->w;
    int* root = mem->r;
    // initialization
    mem->n = n;
    for( i = 0; i < n+1; i++ ) {
        e[IDX(i,i)] = q[i];
        w[IDX(i,i)] = q[i];
    }

    for( l = 1; l < n+1; l++ ) {
        for( i = 0; i < n-l+1; i++ ) {
            j = i+l;
            t_min = INFINITY;
            w_cur = w[IDX(i,j-1)] + p[j-1] + q[j];
            for( r = i; r < j; r++ ) { // l many iterations
                // TODO: check these indices
                t = e[IDX(i,r)] + e[IDX(r+1,j)] + w_cur;
                if( t < t_min ) {
                    t_min = t;
                    r_min = r;
                }
            }
            e[IDX(i,j)] = t_min;
            w[IDX(i,j)] = w_cur;
            root[IDX(i,j)] = r_min;
        }
    }
    return e[IDX(0,n)];
}

size_t bst_get_root_002_scala_rep( void* _bst_obj, size_t i, size_t j )
{
    // [i,j], in table: [i-1, j]+1
    segments_t *mem = _bst_obj;
    size_t n = mem->n;
    int *root = mem->r;
    return (size_t) root[(i-1)*(n+1)+j]+1;
}

void bst_free_002_scala_rep( void* _mem ) {
    segments_t* mem = (segments_t*) _mem;
    free( mem->e );
    free( mem->w );
    free( mem->r );
    free( mem );
}

size_t bst_flops_002_scala_rep( size_t n ) {
    size_t n3 = n*n*n;
    size_t n2 = n*n;
    return (n3 + 5*n2)/2 + 2*n;
}
