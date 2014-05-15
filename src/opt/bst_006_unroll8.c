#include<bst.h>
#include<math.h>
#include<string.h>
#include <assert.h>

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
 * This optimization is based on bst_005. It unrolls the inner-most r-loop
 * by 8.
 */

#define STRIDE (n+1)
#define IDX(i,j) ((i)*STRIDE + j)

typedef struct {
    double* e;
    double* w;
    int* r;
    size_t n;
} segments_t;

void* bst_alloc_006_unroll8( size_t n ) {
    segments_t* mem = (segments_t*) malloc( sizeof(segments_t) );
    size_t sz = (n+1)*(n+1);
    size_t sz2 = ((n+1)*(n+2))/2;
    mem->e = (double*) malloc( sz  * sizeof(double) );
    mem->w = (double*) malloc( sz2 * sizeof(double) );
    mem->r = (int*)    malloc( sz2 * sizeof(int) );
    memset( mem->r, -1, sz2 * sizeof(int) );
    return mem;
}

double bst_compute_006_unroll8( void*_bst_obj, double* p, double* q,
                                   size_t n ) {
    segments_t* mem = (segments_t*) _bst_obj;
    size_t i, l, r, j;
    size_t ue; // unroll_end
    double t, t_min, w_cur;
    double t1, t2, t3, t4, t5, t6, t7, t8;
    double t1_min, t2_min, t3_min, t4_min, t5_min, t6_min, t7_min, t8_min;
    int r_min;
    int r1_min, r2_min, r3_min, r4_min, r5_min, r6_min, r7_min, r8_min;
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
            w_cur = w[w_idx_prev] + p[j-1] + q[j];

            t_min  = INFINITY;
            t1_min = INFINITY;
            t2_min = INFINITY;
            t3_min = INFINITY;
            t4_min = INFINITY;
            t5_min = INFINITY;
            t6_min = INFINITY;
            t7_min = INFINITY;
            t8_min = INFINITY;

            r  = i;
            ue = j-(l%8);
            if (l/8 > 0) { // can unroll by 8
                for (; r < ue; r+=8) {
                    t1 = e[IDX(i,r)]   + e[IDX(j,r+1)] + w_cur;
                    t2 = e[IDX(i,r+1)] + e[IDX(j,r+2)] + w_cur;
                    t3 = e[IDX(i,r+2)] + e[IDX(j,r+3)] + w_cur;
                    t4 = e[IDX(i,r+3)] + e[IDX(j,r+4)] + w_cur;
                    t5 = e[IDX(i,r+4)] + e[IDX(j,r+5)] + w_cur;
                    t6 = e[IDX(i,r+5)] + e[IDX(j,r+6)] + w_cur;
                    t7 = e[IDX(i,r+6)] + e[IDX(j,r+7)] + w_cur;
                    t8 = e[IDX(i,r+7)] + e[IDX(j,r+8)] + w_cur;

                    if (t1 < t1_min ) { t1_min = t1; r1_min = r;   }
                    if (t2 < t2_min ) { t2_min = t2; r2_min = r+1; }
                    if (t3 < t3_min ) { t3_min = t3; r3_min = r+2; }
                    if (t4 < t4_min ) { t4_min = t4; r4_min = r+3; }
                    if (t5 < t5_min ) { t5_min = t5; r5_min = r+4; }
                    if (t6 < t6_min ) { t6_min = t6; r6_min = r+5; }
                    if (t7 < t7_min ) { t7_min = t7; r7_min = r+6; }
                    if (t8 < t8_min ) { t8_min = t8; r8_min = r+7; }
                }

                if (t2_min < t1_min) { t1_min = t2_min; r1_min = r2_min; }
                if (t4_min < t3_min) { t3_min = t4_min; r3_min = r4_min; }
                if (t6_min < t5_min) { t5_min = t6_min; r5_min = r6_min; }
                if (t8_min < t7_min) { t7_min = t8_min; r7_min = r8_min; }
                if (t3_min < t1_min) { t1_min = t3_min; r1_min = r3_min; }
                if (t7_min < t5_min) { t5_min = t7_min; r5_min = r7_min; }
                if (t5_min < t1_min) {
                    t_min = t5_min;
                    r_min = r5_min;
                } else {
                    t_min = t1_min;
                    r_min = r1_min;
                }
            }
            // cleanup of unrolling
            for (; r < j; ++r) {
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
size_t bst_get_root_006_unroll8( void* _bst_obj, size_t i, size_t j )
{
    // [i,j], in table: [i-1, j]+1
    segments_t *mem = _bst_obj;
    size_t n = mem->n;
    int *root = mem->r;
    return (size_t) root[(i-1)*(n+1)+j]+1;
}

void bst_free_006_unroll8( void* _mem ) {
    segments_t* mem = (segments_t*) _mem;
    free( mem->e );
    free( mem->w );
    free( mem->r );
    free( mem );
}
