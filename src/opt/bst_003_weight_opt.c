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

/*
 * In this optimization, the memory layout of the weight matrix w was changed.
 * The first "line" of the w-array now contains the total weight of the tree
 * from 0 to 0, from 1 to 0, ... from n to n. The second line contains the
 * respective values for the tree from 0 to 1, 1 to 2, etc.
 */

#define STRIDE (n+1)
#define IDX(i,j) ((i)*STRIDE + j)

typedef struct {
    double* e;
    double* w;
    int* r;
    size_t n;
} segments_t;

void* bst_alloc_003_weight_opt( size_t n ) {
    segments_t* mem = (segments_t*) malloc( sizeof(segments_t) );
    size_t sz = (n+1)*(n+1);
    size_t sz2 = ((n+1)*(n+2))/2;
    mem->e = (double*) malloc( sz  * sizeof(double) );
    mem->w = (double*) malloc( sz2 * sizeof(double) );
    mem->r = (int*)    malloc( sz2 * sizeof(int) );
    memset( mem->r, -1, sz2 * sizeof(int) );
    return mem;
}

double bst_compute_003_weight_opt( void*_bst_obj, double* p, double* q,
                                   size_t n ) {
    segments_t* mem = (segments_t*) _bst_obj;
    size_t i, l, r, j;
    double t, t_min, w_cur;
    int r_min;
    double* e = mem->e, *w = mem->w;
    int* root = mem->r;
    int w_idx_prev, w_idx_cur;
    // initialization
    mem->n = n;
    for( i = 0; i < n+1; i++ ) {
        e[IDX(i,i)] = q[i];
        w[i]        = q[i]; // weight optimization
        // note that the value w[n] is not going to be read since it is
        // synonymous with q[n]
    }
    w_idx_cur  = n+1; // in the first round, we set the values on the second row
    w_idx_prev = 0; // and we read the values from the first row
    for( l = 1; l < n+1; l++ ) {
        for( i = 0; i < n-l+1; i++ ) {
            j = i+l;
            t_min = INFINITY;
            w_cur = w[w_idx_prev] + p[j-1] + q[j];
            for( r = i; r < j; r++ ) { // l many iterations
                t = e[IDX(i,r)] + e[IDX(j,r+1)] + w_cur;
                if( t < t_min ) {
                    t_min = t; // e[IDX(i,j)] = t;
                    r_min = r; // root[IDX(i,j)] = r;
                }
            }
            w[w_idx_cur]    = w_cur;
            e[IDX(i,j)]     = t_min;
            e[IDX(j,i)]     = t_min;
            root[w_idx_cur] = r_min;
            w_idx_cur  += 1;
            w_idx_prev += 1;
        }
        // the previous row in the weight matrix is longer by one entry
        w_idx_prev += 1;
    }
    return e[IDX(0,n)];
}


// TODO: This is INCORRECT
size_t bst_get_root_003_weight_opt( void* _bst_obj, size_t i, size_t j )
{
    // [i,j], in table: [i-1, j]+1
    segments_t *mem = _bst_obj;
    size_t n = mem->n;
    int *root = mem->r;
    return (size_t) root[(i-1)*(n+1)+j]+1;
}

void bst_free_003_weight_opt( void* _mem ) {
    segments_t* mem = (segments_t*) _mem;
    free( mem->e );
    free( mem->w );
    free( mem->r );
    free( mem );
}

size_t bst_flops_003_weight_opt( size_t n ) {
    size_t n3 = n*n*n;
    size_t n2 = n*n;
    return (n3 + 5*n2)/2 + 2*n;
}
