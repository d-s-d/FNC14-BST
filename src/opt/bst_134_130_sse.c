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
 * AVX x 2 to 131
 */

#define STRIDE (n+1)
#define IDX(i,j) ((i)*STRIDE + j)

typedef struct {
    double* e;
    double* w;
    int* r;
    size_t n;
} segments_t;

void* bst_alloc_134_130_sse( size_t n ) {
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

double bst_compute_134_130_sse( void*_bst_obj, double* p, double* q, size_t nn ) {
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


            // now we can work in junk of 8s
            int idx_rj = IDX(r+1, j);
            int idx_ij = IDX(i, j);
            
            double e_ir_s = e[IDX(i,r)];
            __m128d e_ir  = _mm_set_pd(e_ir_s, e_ir_s);
            __m128  rs    = _mm_castsi128_ps(_mm_set_epi32(r, r, r, r));

            for (; j < (n+1-7); j += 8, idx_rj += 8, idx_ij += 8) {

                __m128d w_ij1 = _mm_load_pd( &w[idx_ij+0] );
                __m128d w_ij2 = _mm_load_pd( &w[idx_ij+2] );
                __m128d w_ij3 = _mm_load_pd( &w[idx_ij+4] );
                __m128d w_ij4 = _mm_load_pd( &w[idx_ij+6] );

                __m128d t1_1  = _mm_add_pd(e_ir, w_ij1);
                __m128d e_rj1 = _mm_load_pd( &e[idx_rj+0] );
                __m128d t2_1  = _mm_add_pd(e_ir, w_ij2);
                __m128d e_rj2 = _mm_load_pd( &e[idx_rj+2] );
                __m128d t3_1  = _mm_add_pd(e_ir, w_ij3);
                __m128d e_rj3 = _mm_load_pd( &e[idx_rj+4] );
                __m128d t4_1  = _mm_add_pd(e_ir, w_ij4);
                __m128d e_rj4 = _mm_load_pd( &e[idx_rj+6] );

                __m128d t1    = _mm_add_pd(t1_1, e_rj1);
                __m128d e_ij1 = _mm_load_pd( &e[idx_ij+0] );
                __m128d t2    = _mm_add_pd(t2_1, e_rj2);
                __m128d e_ij2 = _mm_load_pd( &e[idx_ij+2] );
                __m128d t3    = _mm_add_pd(t3_1, e_rj3);
                __m128d e_ij3 = _mm_load_pd( &e[idx_ij+4] );
                __m128d t4    = _mm_add_pd(t4_1, e_rj4);
                __m128d e_ij4 = _mm_load_pd( &e[idx_ij+6] );

                __m128d cmp1  = _mm_cmplt_pd(t1, e_ij1);
                __m128d cmp2  = _mm_cmplt_pd(t2, e_ij2);
                __m128d cmp3  = _mm_cmplt_pd(t3, e_ij3);
                __m128d cmp4  = _mm_cmplt_pd(t4, e_ij4);

                /*
                __m128d res1  = _mm_or_pd(_mm_and_pd(cmp1, t1), _mm_andnot_pd(cmp1, e_ij1));
                __m128d res2  = _mm_or_pd(_mm_and_pd(cmp2, t2), _mm_andnot_pd(cmp2, e_ij2));
                __m128d res3  = _mm_or_pd(_mm_and_pd(cmp3, t3), _mm_andnot_pd(cmp3, e_ij3));
                __m128d res4  = _mm_or_pd(_mm_and_pd(cmp4, t4), _mm_andnot_pd(cmp4, e_ij4));

                _mm_store_pd( &e[idx_ij+0], res1);
                _mm_store_pd( &e[idx_ij+2], res2);
                _mm_store_pd( &e[idx_ij+4], res3);
                _mm_store_pd( &e[idx_ij+6], res4);
                */

                _mm_maskstore_pd( &e[idx_ij+0], _mm_castpd_si128(cmp1), t1);
                _mm_maskstore_pd( &e[idx_ij+2], _mm_castpd_si128(cmp2), t2);
                _mm_maskstore_pd( &e[idx_ij+4], _mm_castpd_si128(cmp3), t3);
                _mm_maskstore_pd( &e[idx_ij+6], _mm_castpd_si128(cmp4), t4);

                __m128i rm1 = _mm_castps_si128(_mm_shuffle_ps(
                        _mm_castpd_ps(cmp1), _mm_castpd_ps(cmp2),
                        _MM_SHUFFLE(0,2,0,2)));

                __m128i rm2 = _mm_castps_si128(_mm_shuffle_ps(
                        _mm_castpd_ps(cmp3), _mm_castpd_ps(cmp4),
                        _MM_SHUFFLE(0,2,0,2)));

                _mm_maskstore_ps( &root[idx_ij+0], rm1, rs);
                _mm_maskstore_ps( &root[idx_ij+4], rm2, rs);
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

size_t bst_get_root_134_130_sse( void* _bst_obj, size_t i, size_t j )
{
    // [i,j], in table: [i-1, j]+1
    segments_t *mem = _bst_obj;
    size_t n = mem->n;
    int *root = mem->r;
    return (size_t) root[(i-1)*(n+1)+j]+1;
}

void bst_free_134_130_sse( void* _mem ) {
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

size_t bst_flops_134_130_sse( size_t n ) {
    size_t n3 = n*n*n;
    size_t n2 = n*n;
    return (size_t) ( n3/3.0 + 2*n2 + 5.0*n/3 );
}
