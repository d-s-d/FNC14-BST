#include<bst.h>
#include<math.h>
#include<string.h>
#include <assert.h>
#include <x86intrin.h>
#include <stdint.h>


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
    int64_t* r;
    size_t n;
} segments_t;

void* bst_alloc_144_141_avx_8x4( size_t n ) {
    segments_t* mem = (segments_t*) malloc( sizeof(segments_t) );
    size_t sz = (n+1)*(n+1);
    // 32B alignment for AVX.
    mem->e = _mm_malloc(sz * sizeof(double), 32); assert(mem->e);
    mem->w = _mm_malloc(sz * sizeof(double), 32); assert(mem->w);
    mem->r = _mm_malloc(sz * sizeof(int64_t), 32);    assert(mem->r);
    mem->n = n;
    memset( mem->r, -1, sz * sizeof(int64_t) );
    return mem;
}

#define NB 8
#if (NB!=8)
#error this version requires NB=8
#endif
double bst_compute_144_141_avx_8x4( void*_bst_obj, double* p, double* q, size_t nn ) {
    segments_t* mem = (segments_t*) _bst_obj;
    int n;
    int i, l, j;
    int64_t r;
    double t, t_min, w_cur;
    int64_t r_min;
    double* e = mem->e, *w = mem->w;
    int64_t* root = mem->r;
    // initialization
    mem->n = nn;
    n = nn; // subtractions with n potentially negative. say hello to all the bugs
    int idx = 0;
    for( i = 0; i < n+1; i++, idx += (n+2)) {
        e[idx] = q[i];
        w[idx] = q[i];
    }

    int ib;

    // compute bottom right triangle NBxNB (i.e. bottom NB rows)
    i = n-1;
    int idx_ii = IDX(i,i);
    for (; (i >= 0) && (i > (n-NB)); --i, idx_ii -= (n+2) ) {
        int idx_ij = idx_ii+1;
        for (j = i+1; j < n+1; ++j, idx_ij += 1) {
            w[idx_ij] = w[idx_ij-1] + p[j-1] + q[j];
            t_min = INFINITY;
            int idx_ir = idx_ii;
            int idx_r1j = idx_ii + STRIDE + (j-i);
            for (r=i; r<j; ++r, idx_ir += 1, idx_r1j += STRIDE) {
                t = e[idx_ir] + e[idx_r1j] + w[idx_ij];
                if (t < t_min) {
                    t_min = t;
                    r_min = r;
                }
            }
            e[idx_ij]    = t_min;
            root[idx_ij] = r_min;
        }
    }

    // compute the remaining rows
    // printf("Rest of rows: i=%d\n", i);
    for (; i >= 0; --i, idx_ii -= (n+2)) {
        // First, the starting NB values in this row are computed.
        // This corresponds to completing the NBxNB triangle right down from (i,i)
        int idx_ij = idx_ii+1;
        for (j = i+1; j < (i+NB); ++j, idx_ij += 1) {
            w[idx_ij] = w[idx_ij-1] + p[j-1] + q[j];
            t_min = INFINITY;
            int idx_ir = idx_ii;
            int idx_r1j = idx_ii + STRIDE + (j-i);
            for (r=i; r<j; ++r, idx_ir += 1, idx_r1j += STRIDE) {
                t = e[idx_ir] + e[idx_r1j] + w[idx_ij];
                if (t < t_min) {
                    t_min = t;
                    r_min = r;
                }
            }
            e[idx_ij]    = t_min;
            root[idx_ij] = r_min;
        }

        // Now we compute the rest of the row, but do updates in chunk of NB's.
        // Since we now have the first NB values in this row, we can compute
        // the first NB iterations of the r-loop for all the remaining values
        // in this row.
        // (this needs to be seperated from the loop afterwards since we also do
        //  initialization to INFINITY))
        for (; j < (n+1); ++j, idx_ij += 1) {
            w[idx_ij] = w[idx_ij-1] + p[j-1] + q[j];
            t_min = INFINITY;
            int idx_ir = idx_ii;
            int idx_r1j = idx_ii + STRIDE + (j-i);
            for (r=i; r < (i+NB); ++r, idx_ir +=1, idx_r1j += STRIDE) {
                t = e[idx_ir] + e[idx_r1j] + w[idx_ij];
                if (t < t_min) {
                    t_min = t;
                    r_min = r;
                }
            }
            e[idx_ij]    = t_min;
            root[idx_ij] = r_min;
        }

        // We now continue to update the values in this row in chunks of NB
        // as long as possible.
        const int NB_STRIDE = NB * (STRIDE+1);
        int idx_ib_x = idx_ii + NB_STRIDE;
        for (ib = i+NB; (ib+NB) < (n+1); ib += NB, idx_ib_x += NB_STRIDE) {
            //printf("got in here for i=%d\n", i);

            // Again we start by finishing computing the next NB values of
            // row 'i'. The last values needed for that are from the triangle
            // in down right from (ib,ib).
            idx_ij = idx_ii + (ib+1) - i;
            for (j = (ib+1); j < (ib+NB); ++j, idx_ij += 1) {
                t_min = e[idx_ij];
                r_min = root[idx_ij];
                int idx_ir = idx_ii + ib - i;
                int idx_r1j = idx_ib_x + STRIDE + j - ib;
                for (r=ib; r<j; ++r, idx_ir += 1, idx_r1j += STRIDE) {
                    t = e[idx_ir] + e[idx_r1j] + w[idx_ij];
                    if (t < t_min) {
                        t_min = t;
                        r_min = r;
                    }
                }
                e[idx_ij]    = t_min;
                root[idx_ij] = r_min;
            }

            // Now, having NB new values in row 'i', we compute the next NB
            // r-iterations for the remaining values in row 'i'.
            
            // make this 4-aligned
            for (; (j < (n+1)) && (j%4); ++j, idx_ij += 1) {
                t_min = e[idx_ij];
                r_min = root[idx_ij];
                int idx_ir = idx_ii + ib - i;
                int idx_r1j = idx_ib_x + STRIDE + j - ib;
                for (r=ib; r<(ib+NB); ++r, idx_ir += 1, idx_r1j += STRIDE) {
                    t = e[idx_ir] + e[idx_r1j] + w[idx_ij];
                    if (t < t_min) {
                        t_min = t;
                        r_min = r;
                    }
                }
                e[idx_ij]    = t_min;
                root[idx_ij] = r_min;
            }

            // run in chunks of 8x4
            int idx_ir_here = idx_ii+ib-i;
            __m256d v_e_ir0 = _mm256_set1_pd(e[idx_ir_here+0]);
            __m256d v_e_ir1 = _mm256_set1_pd(e[idx_ir_here+1]);
            __m256d v_e_ir2 = _mm256_set1_pd(e[idx_ir_here+2]);
            __m256d v_e_ir3 = _mm256_set1_pd(e[idx_ir_here+3]);
            __m256d v_e_ir4 = _mm256_set1_pd(e[idx_ir_here+4]);
            __m256d v_e_ir5 = _mm256_set1_pd(e[idx_ir_here+5]);
            __m256d v_e_ir6 = _mm256_set1_pd(e[idx_ir_here+6]);
            __m256d v_e_ir7 = _mm256_set1_pd(e[idx_ir_here+7]);
            __m256i v_r0    = _mm256_set1_epi64x(ib+0);
            __m256i v_r1    = _mm256_set1_epi64x(ib+1);
            __m256i v_r2    = _mm256_set1_epi64x(ib+2);
            __m256i v_r3    = _mm256_set1_epi64x(ib+3);
            __m256i v_r4    = _mm256_set1_epi64x(ib+4);
            __m256i v_r5    = _mm256_set1_epi64x(ib+5);
            __m256i v_r6    = _mm256_set1_epi64x(ib+6);
            __m256i v_r7    = _mm256_set1_epi64x(ib+7);
            for (; j < (n+1-3); j += 4, idx_ij += 4) {
                // all the loads
                __m256d t_min0 = _mm256_loadu_pd( &e[idx_ij+0] );
                __m256i r_min0 = _mm256_lddqu_si256( &root[idx_ij+0] );
                __m256d w0     = _mm256_loadu_pd( &w[idx_ij+0] );

                int ir_start_0 = idx_ib_x + STRIDE + j - ib;
                int ir_start_1 = ir_start_0 + STRIDE;
                int ir_start_2 = ir_start_1 + STRIDE;
                int ir_start_3 = ir_start_2 + STRIDE;
                int ir_start_4 = ir_start_3 + STRIDE;
                int ir_start_5 = ir_start_4 + STRIDE;
                int ir_start_6 = ir_start_5 + STRIDE;
                int ir_start_7 = ir_start_6 + STRIDE;
                __m256d e_ir00 = _mm256_load_pd( &e[ir_start_0+0] );
                __m256d e_ir10 = _mm256_load_pd( &e[ir_start_1+0] );
                __m256d e_ir20 = _mm256_load_pd( &e[ir_start_2+0] );
                __m256d e_ir30 = _mm256_load_pd( &e[ir_start_3+0] );
                __m256d e_ir40 = _mm256_load_pd( &e[ir_start_4+0] );
                __m256d e_ir50 = _mm256_load_pd( &e[ir_start_5+0] );
                __m256d e_ir60 = _mm256_load_pd( &e[ir_start_6+0] );
                __m256d e_ir70 = _mm256_load_pd( &e[ir_start_7+0] );

                // do the comps
                __m256d t00_1 = _mm256_add_pd(v_e_ir0, w0);
                __m256d t00_2 = _mm256_add_pd(t00_1, e_ir00);

                __m256d t10_1 = _mm256_add_pd(v_e_ir1, w0);
                __m256d t10_2 = _mm256_add_pd(t10_1, e_ir10);

                __m256d t20_1 = _mm256_add_pd(v_e_ir2, w0);
                __m256d t20_2 = _mm256_add_pd(t20_1, e_ir20);

                __m256d t30_1 = _mm256_add_pd(v_e_ir3, w0);
                __m256d t30_2 = _mm256_add_pd(t30_1, e_ir30);

                __m256d t40_1 = _mm256_add_pd(v_e_ir4, w0);
                __m256d t40_2 = _mm256_add_pd(t40_1, e_ir40);

                __m256d t50_1 = _mm256_add_pd(v_e_ir5, w0);
                __m256d t50_2 = _mm256_add_pd(t50_1, e_ir50);

                __m256d t60_1 = _mm256_add_pd(v_e_ir6, w0);
                __m256d t60_2 = _mm256_add_pd(t60_1, e_ir60);

                __m256d t70_1 = _mm256_add_pd(v_e_ir7, w0);
                __m256d t70_2 = _mm256_add_pd(t70_1, e_ir70);

                // mins
                __m256d cmp0    = _mm256_cmp_pd(t00_2, t_min0 , _CMP_LT_OQ);
                __m256d t_min01 = _mm256_blendv_pd(t_min0,  t00_2, cmp0);
                __m256i r_min01 = _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(r_min0),  _mm256_castsi256_pd(v_r0), cmp0));
                __m256d cmp1    = _mm256_cmp_pd(t10_2, t_min01, _CMP_LT_OQ);
                __m256d t_min02 = _mm256_blendv_pd(t_min01, t10_2, cmp1);
                __m256i r_min02 = _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(r_min01),  _mm256_castsi256_pd(v_r1), cmp1));
                __m256d cmp2    = _mm256_cmp_pd(t20_2, t_min02, _CMP_LT_OQ);
                __m256d t_min03 = _mm256_blendv_pd(t_min02, t20_2, cmp2);
                __m256i r_min03 = _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(r_min02),  _mm256_castsi256_pd(v_r2), cmp2));
                __m256d cmp3    = _mm256_cmp_pd(t30_2, t_min03, _CMP_LT_OQ);
                __m256d t_min04 = _mm256_blendv_pd(t_min03, t30_2, cmp3);
                __m256i r_min04 = _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(r_min03),  _mm256_castsi256_pd(v_r3), cmp3));
                __m256d cmp4    = _mm256_cmp_pd(t40_2, t_min04, _CMP_LT_OQ);
                __m256d t_min05 = _mm256_blendv_pd(t_min04, t40_2, cmp4);
                __m256i r_min05 = _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(r_min04),  _mm256_castsi256_pd(v_r4), cmp4));
                __m256d cmp5    = _mm256_cmp_pd(t50_2, t_min05, _CMP_LT_OQ);
                __m256d t_min06 = _mm256_blendv_pd(t_min05, t50_2, cmp5);
                __m256i r_min06 = _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(r_min05),  _mm256_castsi256_pd(v_r5), cmp5));
                __m256d cmp6    = _mm256_cmp_pd(t60_2, t_min06, _CMP_LT_OQ);
                __m256d t_min07 = _mm256_blendv_pd(t_min06, t60_2, cmp6);
                __m256i r_min07 = _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(r_min06),  _mm256_castsi256_pd(v_r6), cmp6));
                __m256d cmp7    = _mm256_cmp_pd(t70_2, t_min07, _CMP_LT_OQ);
                __m256d t_min08 = _mm256_blendv_pd(t_min07, t70_2, cmp7);
                __m256i r_min08 = _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(r_min07),  _mm256_castsi256_pd(v_r7), cmp7));

                // writeback
                _mm256_storeu_pd( &e[idx_ij+0], t_min08);
                _mm256_storeu_si256( &root[idx_ij+0], r_min08);
            }


            // non-4 cleanup
            for (; j < (n+1); ++j, idx_ij += 1) {
                t_min = e[idx_ij];
                r_min = root[idx_ij];
                int idx_ir = idx_ii + ib - i;
                int idx_r1j = idx_ib_x + STRIDE + j - ib;
                for (r=ib; r<(ib+NB); ++r, idx_ir += 1, idx_r1j += STRIDE) {
                    t = e[idx_ir] + e[idx_r1j] + w[idx_ij];
                    if (t < t_min) {
                        t_min = t;
                        r_min = r;
                    }
                }
                e[idx_ij]    = t_min;
                root[idx_ij] = r_min;
            }
        }

        // There are less than NB elements remaining in row 'i'. The values
        // missing come from the triangle down right of (ib,ib)
        idx_ij = idx_ii + (ib+1) - i;
        for (j = (ib+1); j < (n+1); ++j, idx_ij += 1) {
            t_min = e[idx_ij];
            r_min = root[idx_ij];
            int idx_ir = idx_ii + ib - i;
            int idx_r1j = idx_ib_x + STRIDE + j - ib;
            for (r=ib; r<j; ++r, idx_ir += 1, idx_r1j += STRIDE) {
                t = e[idx_ir] + e[idx_r1j] + w[idx_ij];
                if (t < t_min) {
                    t_min = t;
                    r_min = r;
                }
            }
            e[idx_ij]    = t_min;
            root[idx_ij] = r_min;
        }
    }

    return e[n]; //e[IDX(0,n)];
}

size_t bst_get_root_144_141_avx_8x4( void* _bst_obj, size_t i, size_t j )
{
    // [i,j], in table: [i-1, j]+1
    segments_t *mem = _bst_obj;
    size_t n = mem->n;
    int64_t *root = mem->r;
    return (size_t) root[(i-1)*(n+1)+j]+1;
}

void bst_free_144_141_avx_8x4( void* _mem ) {
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

size_t bst_flops_144_141_avx_8x4( size_t n ) {
    size_t n3 = n*n*n;
    size_t n2 = n*n;
    return (size_t) ( n3/3.0 + 2*n2 + 5.0*n/3 );
}
