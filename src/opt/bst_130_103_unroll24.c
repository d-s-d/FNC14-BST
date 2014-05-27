#include <bst.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <x86intrin.h>


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
 * loop unrolling by 8 to 103
 */

#define STRIDE (n+1)
#define IDX(i,j) ((i)*STRIDE + j)

typedef struct {
    double* e;
    double* w;
    int* r;
    size_t n;
} segments_t;

void* bst_alloc_130_103_unroll24( size_t n ) {
    segments_t* mem = (segments_t*) malloc( sizeof(segments_t) );
    size_t sz = (n+1)*(n+1);
    // 32B alignment for AVX.
    mem->e = _mm_malloc(sz * sizeof(double), 32); assert(mem->e);
    mem->w = _mm_malloc(sz * sizeof(double), 32); assert(mem->w);
    mem->r = _mm_malloc(sz * sizeof(int), 32);    assert(mem->r);
    mem->n = n;
    memset( mem->r, -1, sz * sizeof(int) );
    return mem;
}

double bst_compute_130_103_unroll24( void*_bst_obj, double* p, double* q, size_t nn ) {
    segments_t* mem = (segments_t*) _bst_obj;
    int n;
    int i, l, r, j;
    double t, t_min, w_cur;
    int r_min;
    double* e = mem->e, *w = mem->w;
    int* root = mem->r;
    // initialization
    // mem->n = nn;
    n = nn; // subtractions with n potentially negative. say hello to all the bugs

    e[IDX(n,n)] = q[n];
    for (i = n-1; i >= 0; --i) {
        e[IDX(i,i)] = q[i];
        w[IDX(i,i)] = q[i];
        for (j = i+1; j < n+1; ++j) {
            e[IDX(i,j)] = INFINITY;
            w[IDX(i,j)] = w[IDX(i,j-1)] + p[j-1] + q[j];
        }

        for (r = i; r < n; ++r) {
            // we use 4-way AVX packed double instructions.
            // we require 32B aligned loads, i.e. j%4==0
            for (j = r+1; (j < n+1) && (j%4 != 0); ++j ) {
                t = e[IDX(i,r)] + e[IDX(r+1,j)] + w[IDX(i,j)];
                if (t < e[IDX(i,j)]) {
                    e[IDX(i,j)] = t;
                    root[IDX(i,j)] = r;
                }
            }

            double ir = e[IDX(i,r)];
            int idx_rj = IDX(r+1, j);
            int idx_ij = IDX(i, j);
#define DJ 24
#define T_BLOCK(dj) double t##dj = ir + e[idx_rj+dj] + w[idx_ij+dj]
#define M_BLOCK(dj) if (t##dj < e[idx_ij+dj]) { e[idx_ij+dj] = t##dj; root[idx_ij+dj] = r; }
            for (; j < (n+1-(DJ-1)); j += DJ, idx_rj += DJ, idx_ij += DJ) {
                // for (int jx = j; jx < (j+8); ++jx) {
                //     t = e[IDX(i,r)] + e[IDX(r+1,jx)] + w[IDX(i,jx)];
                //     if (t < e[IDX(i,jx)]) {
                //         e[IDX(i,jx)] = t;
                //         root[IDX(i,jx)] = r;
                //     }
                // }

                // double t1 = ir + e[idx_rj+0] + w[idx_ij+0];
                // double t2 = ir + e[idx_rj+1] + w[idx_ij+1];
                // double t3 = ir + e[idx_rj+2] + w[idx_ij+2];
                // double t4 = ir + e[idx_rj+3] + w[idx_ij+3];
                // double t5 = ir + e[idx_rj+4] + w[idx_ij+4];
                // double t6 = ir + e[idx_rj+5] + w[idx_ij+5];
                // double t7 = ir + e[idx_rj+6] + w[idx_ij+6];
                // double t8 = ir + e[idx_rj+7] + w[idx_ij+7];

                // if (t1 < e[idx_ij+0]) { e[idx_ij+0] = t1; root[idx_ij+0] = r; }
                // if (t2 < e[idx_ij+1]) { e[idx_ij+1] = t2; root[idx_ij+1] = r; }
                // if (t3 < e[idx_ij+2]) { e[idx_ij+2] = t3; root[idx_ij+2] = r; }
                // if (t4 < e[idx_ij+3]) { e[idx_ij+3] = t4; root[idx_ij+3] = r; }
                // if (t5 < e[idx_ij+4]) { e[idx_ij+4] = t5; root[idx_ij+4] = r; }
                // if (t6 < e[idx_ij+5]) { e[idx_ij+5] = t6; root[idx_ij+5] = r; }
                // if (t7 < e[idx_ij+6]) { e[idx_ij+6] = t7; root[idx_ij+6] = r; }
                // if (t8 < e[idx_ij+7]) { e[idx_ij+7] = t8; root[idx_ij+7] = r; }

                T_BLOCK(0);
                T_BLOCK(1);
                T_BLOCK(2);
                T_BLOCK(3);
                T_BLOCK(4);
                T_BLOCK(5);
                T_BLOCK(6);
                T_BLOCK(7);
                T_BLOCK(8);
                T_BLOCK(9);
                T_BLOCK(10);
                T_BLOCK(11);
                T_BLOCK(12);
                T_BLOCK(13);
                T_BLOCK(14);
                T_BLOCK(15);
                T_BLOCK(16);
                T_BLOCK(17);
                T_BLOCK(18);
                T_BLOCK(19);
                T_BLOCK(20);
                T_BLOCK(21);
                T_BLOCK(22);
                T_BLOCK(23);

                M_BLOCK(0);
                M_BLOCK(1);
                M_BLOCK(2);
                M_BLOCK(3);
                M_BLOCK(4);
                M_BLOCK(5);
                M_BLOCK(6);
                M_BLOCK(7);
                M_BLOCK(8);
                M_BLOCK(9);
                M_BLOCK(10);
                M_BLOCK(11);
                M_BLOCK(12);
                M_BLOCK(13);
                M_BLOCK(14);
                M_BLOCK(15);
                M_BLOCK(16);
                M_BLOCK(17);
                M_BLOCK(18);
                M_BLOCK(19);
                M_BLOCK(20);
                M_BLOCK(21);
                M_BLOCK(22);
                M_BLOCK(23);
            }

            // and do the non-8 cleanup
            for (; j < n+1; ++j) {
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

size_t bst_get_root_130_103_unroll24( void* _bst_obj, size_t i, size_t j )
{
    // [i,j], in table: [i-1, j]+1
    segments_t *mem = _bst_obj;
    size_t n = mem->n;
    int *root = mem->r;
    return (size_t) root[(i-1)*(n+1)+j]+1;
}

void bst_free_130_103_unroll24( void* _mem ) {
    segments_t* mem = (segments_t*) _mem;

    /*
    size_t n = mem->n;
    for (size_t i=0; i<=n; ++i) {
        for (size_t j=0; j<=n; ++j) {
            printf(" %.3lf", mem->e[i*(n+1)+j]);
        }
        printf("\n");
    }
    */

    _mm_free( mem->e );
    _mm_free( mem->w );
    _mm_free( mem->r );
    free( mem );
}

size_t bst_flops_130_103_unroll24( size_t n ) {
    size_t n3 = n*n*n;
    size_t n2 = n*n;
    return (size_t) ( n3/3.0 + 2*n2 + 5.0*n/3 );
}
