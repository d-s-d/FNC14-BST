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

void* bst_alloc_132_131_double( size_t n ) {
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

double bst_compute_132_131_double( void*_bst_obj, double* p, double* q, size_t nn ) {
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
            double ir = e[IDX(i,r)];
            int idx_rj = IDX(r+1, j);
            int idx_ij = IDX(i, j);
            
            double e_ir_s = e[IDX(i,r)];
            __m256d e_ir  = _mm256_set_pd(e_ir_s, e_ir_s, e_ir_s, e_ir_s);
            __m256  rs    = _mm256_castsi256_ps(_mm256_set_epi32(r, r, r, r, r, r, r, r));

            for (; j < (n+1-15); j += 16, idx_rj += 16, idx_ij += 16) {
                // 8-way
                __m256d u1_e_rj1 = _mm256_load_pd(e + idx_rj + 0);
                __m256d u1_e_ij1 = _mm256_load_pd(e + idx_ij + 0);
                __m256d u1_w_ij1 = _mm256_load_pd(w + idx_ij + 0);
                __m256d u1_e_rj2 = _mm256_load_pd(e + idx_rj + 4);
                __m256d u1_e_ij2 = _mm256_load_pd(e + idx_ij + 4);
                __m256d u1_w_ij2 = _mm256_load_pd(w + idx_ij + 4);

                __m256d u1_t1_1 = _mm256_add_pd(e_ir, u1_e_rj1);
                __m256d u1_t1   = _mm256_add_pd(u1_t1_1, u1_w_ij1);
                __m256d u1_t2_1 = _mm256_add_pd(e_ir, u1_e_rj2);
                __m256d u1_t2   = _mm256_add_pd(u1_t2_1, u1_w_ij2);

                __m256d u1_cmp1 = _mm256_cmp_pd(u1_t1, u1_e_ij1, _CMP_LT_OQ);
                __m256d u1_cmp2 = _mm256_cmp_pd(u1_t2, u1_e_ij2, _CMP_LT_OQ);

                _mm256_maskstore_pd(e + idx_ij + 0, _mm256_castpd_si256(u1_cmp1), u1_t1);
                _mm256_maskstore_pd(e + idx_ij + 4, _mm256_castpd_si256(u1_cmp2), u1_t2);

                __m256 u1_rootmask = _mm256_insertf128_ps(
                        _mm256_castps128_ps256(_mm256_cvtpd_ps(u1_cmp1)),
                        _mm256_cvtpd_ps(u1_cmp2), 1);
                _mm256_maskstore_ps(root + idx_ij + 0, _mm256_castps_si256(u1_rootmask), rs);

                // 8-way second
                __m256d u2_e_rj1 = _mm256_load_pd(e + idx_rj + 8);
                __m256d u2_e_ij1 = _mm256_load_pd(e + idx_ij + 8);
                __m256d u2_w_ij1 = _mm256_load_pd(w + idx_ij + 8);
                __m256d u2_e_rj2 = _mm256_load_pd(e + idx_rj + 12);
                __m256d u2_e_ij2 = _mm256_load_pd(e + idx_ij + 12);
                __m256d u2_w_ij2 = _mm256_load_pd(w + idx_ij + 12);

                __m256d u2_t1_1 = _mm256_add_pd(e_ir, u2_e_rj1);
                __m256d u2_t1   = _mm256_add_pd(u2_t1_1, u2_w_ij1);
                __m256d u2_t2_1 = _mm256_add_pd(e_ir, u2_e_rj2);
                __m256d u2_t2   = _mm256_add_pd(u2_t2_1, u2_w_ij2);

                __m256d u2_cmp1 = _mm256_cmp_pd(u2_t1, u2_e_ij1, _CMP_LT_OQ);
                __m256d u2_cmp2 = _mm256_cmp_pd(u2_t2, u2_e_ij2, _CMP_LT_OQ);

                _mm256_maskstore_pd(e + idx_ij + 8, _mm256_castpd_si256(u2_cmp1), u2_t1);
                _mm256_maskstore_pd(e + idx_ij + 12, _mm256_castpd_si256(u2_cmp2), u2_t2);

                __m256 u2_rootmask = _mm256_insertf128_ps(
                        _mm256_castps128_ps256(_mm256_cvtpd_ps(u2_cmp1)),
                        _mm256_cvtpd_ps(u2_cmp2), 1);
                _mm256_maskstore_ps(root + idx_ij + 8, _mm256_castps_si256(u2_rootmask), rs);
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

size_t bst_get_root_132_131_double( void* _bst_obj, size_t i, size_t j )
{
    // [i,j], in table: [i-1, j]+1
    segments_t *mem = _bst_obj;
    size_t n = mem->n;
    int *root = mem->r;
    return (size_t) root[(i-1)*(n+1)+j]+1;
}

void bst_free_132_131_double( void* _mem ) {
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

size_t bst_flops_132_131_double( size_t n ) {
    size_t n3 = n*n*n;
    size_t n2 = n*n;
    return (size_t) ( n3/3.0 + 2*n2 + 5.0*n/3 );
}
