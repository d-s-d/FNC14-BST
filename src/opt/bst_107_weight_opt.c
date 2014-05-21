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
 * This version applies scalar replacement to bst_100_ref_bottomup.
 */

#define STRIDE (n+1)
#define IDX(i,j) ((i)*STRIDE + j)

typedef struct {
    double* e;
    double* w;
    int* r;
    size_t n;
} segments_t;

void* bst_alloc_107_weight_opt( size_t n ) {
    segments_t* mem = (segments_t*) malloc( sizeof(segments_t) );
    size_t sz = (n+1)*(n+1);
    size_t sz2 = ((n+1)*(n+2))/2;
    mem->e = (double*) malloc( sz * sizeof(double) );
    // mem->w = (double*) malloc( sz2 * sizeof(double) );
    mem->w = NULL;
    mem->r = (int*)    malloc( sz * sizeof(int) );
    memset( mem->r, -1, sz * sizeof(int) );
    return mem;
}

double bst_compute_107_weight_opt( void*_bst_obj, double* p, double* q, size_t n ) {
    segments_t* mem = (segments_t*) _bst_obj;
    int i, l, r, j;
    double t, t_min; // , w_cur;
    int r_min;
    double* e = mem->e; //, *w = mem->w;
    int* root = mem->r;
    double w_cur;
    // initialization
    // mem->n = n;
    for( i = 0; i < n+1; i++ ) {
        e[IDX(i,i)]  = q[i];
    }

    for (i = n-1; i >= 0; --i) {
        w_cur = q[i];
        for (j = i+1; j < n+1; ++j) {
            t_min = INFINITY;
            w_cur = w_cur + p[j-1] + q[j];
            for (r=i; r<j; ++r) {
                t = e[IDX(i,r)] + e[IDX(j,r+1)] + w_cur;
                if (t < t_min) {
                    t_min = t;
                    r_min = r;
                }
            }
            e[IDX(i,j)]    = t_min;
            e[IDX(j,i)]    = t_min;
            root[IDX(i,j)] = r_min;
        }
    }

    return e[IDX(0,n)];
}

size_t bst_get_root_107_weight_opt( void* _bst_obj, size_t i, size_t j )
{
    // [i,j], in table: [i-1, j]+1
    segments_t *mem = _bst_obj;
    size_t n = mem->n;
    int *root = mem->r;
    return (size_t) root[(i-1)*(n+1)+j]+1;
}

void bst_free_107_weight_opt( void* _mem ) {
    segments_t* mem = (segments_t*) _mem;
    free( mem->e );
    // free( mem->w );
    free( mem->r );
    free( mem );
}

size_t bst_flops_107_weight_opt( size_t n ) {
    size_t n3 = n*n*n;
    size_t n2 = n*n;
    return (size_t) ( n3/3.0 + 2*n2 + 5.0*n/3 );
}
