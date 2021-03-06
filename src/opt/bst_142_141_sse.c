#include<bst.h>
#include<math.h>
#include<string.h>
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

void* bst_alloc_142_141_sse( size_t n ) {
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

#define NB 4
#if (NB!=4)
#error this version requires NB=4
#endif
double bst_compute_142_141_sse( void*_bst_obj, double* p, double* q, size_t nn ) {
    segments_t* mem = (segments_t*) _bst_obj;
    int n;
    int i, l, r, j;
    double t, t_min, w_cur;
    int r_min;
    double* e = mem->e, *w = mem->w;
    int* root = mem->r;
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

            // run in chunks of 4
            int idx_ir_here = idx_ii+ib-i;
            double e_ir0 = e[idx_ir_here+0];
            double e_ir1 = e[idx_ir_here+1];
            double e_ir2 = e[idx_ir_here+2];
            double e_ir3 = e[idx_ir_here+3];
            __m128d v_e_ir0 = _mm_set_pd(e_ir0, e_ir0);
            __m128d v_e_ir1 = _mm_set_pd(e_ir1, e_ir1);
            __m128d v_e_ir2 = _mm_set_pd(e_ir2, e_ir2);
            __m128d v_e_ir3 = _mm_set_pd(e_ir3, e_ir3);
            for (; j < (n+1-3); j += 4, idx_ij += 4) {
                // for (int jx = 0; jx < 4; ++jx) {
                //     t_min = e[idx_ij+jx];
                //     r_min = root[idx_ij+jx];
                //     int idx_ir = idx_ii + ib - i;
                //     int idx_r1j = idx_ib_x + STRIDE + j - ib + jx;
                //     for (r=ib; r<(ib+NB); ++r, idx_ir += 1, idx_r1j += STRIDE) {
                //         t = e[idx_ir] + e[idx_r1j] + w[idx_ij +jx];
                //         if (t < t_min) {
                //             t_min = t;
                //             r_min = r;
                //         }
                //     }
                //     e[idx_ij +jx]    = t_min;
                //     root[idx_ij +jx] = r_min;
                // }

                // double t_min0 = e[idx_ij+0];
                // double t_min1 = e[idx_ij+1];
                // double t_min2 = e[idx_ij+2];
                // double t_min3 = e[idx_ij+3];
                // int    r_min0 = root[idx_ij+0];
                // int    r_min1 = root[idx_ij+1];
                // int    r_min2 = root[idx_ij+2];
                // int    r_min3 = root[idx_ij+3];
                // double w0     = w[idx_ij+0];
                // double w1     = w[idx_ij+1];
                // double w2     = w[idx_ij+2];
                // double w3     = w[idx_ij+3];

                // int ir_start_0 = idx_ib_x + STRIDE + j - ib;
                // int ir_start_1 = ir_start_0 + STRIDE;
                // int ir_start_2 = ir_start_1 + STRIDE;
                // int ir_start_3 = ir_start_2 + STRIDE;
                // double e_ir00 = e[ir_start_0+0];
                // double e_ir01 = e[ir_start_0+1];
                // double e_ir02 = e[ir_start_0+2];
                // double e_ir03 = e[ir_start_0+3];
                // double e_ir10 = e[ir_start_1+0];
                // double e_ir11 = e[ir_start_1+1];
                // double e_ir12 = e[ir_start_1+2];
                // double e_ir13 = e[ir_start_1+3];
                // double e_ir20 = e[ir_start_2+0];
                // double e_ir21 = e[ir_start_2+1];
                // double e_ir22 = e[ir_start_2+2];
                // double e_ir23 = e[ir_start_2+3];
                // double e_ir30 = e[ir_start_3+0];
                // double e_ir31 = e[ir_start_3+1];
                // double e_ir32 = e[ir_start_3+2];
                // double e_ir33 = e[ir_start_3+3];

                // double t1_1 = e_ir0 + e_ir00 + w0;
                // double t2_1 = e_ir0 + e_ir01 + w1;
                // double t3_1 = e_ir0 + e_ir02 + w2;
                // double t4_1 = e_ir0 + e_ir03 + w3;
                // t_min0 = (t1_1 < t_min0) ? t1_1 : t_min0;
                // t_min1 = (t2_1 < t_min1) ? t2_1 : t_min1;
                // t_min2 = (t3_1 < t_min2) ? t3_1 : t_min2;
                // t_min3 = (t4_1 < t_min3) ? t4_1 : t_min3;

                // double t1_2 = e_ir1 + e_ir10 + w0;
                // double t2_2 = e_ir1 + e_ir11 + w1;
                // double t3_2 = e_ir1 + e_ir12 + w2;
                // double t4_2 = e_ir1 + e_ir13 + w3;
                // t_min0 = (t1_2 < t_min0) ? t1_2 : t_min0;
                // t_min1 = (t2_2 < t_min1) ? t2_2 : t_min1;
                // t_min2 = (t3_2 < t_min2) ? t3_2 : t_min2;
                // t_min3 = (t4_2 < t_min3) ? t4_2 : t_min3;

                // double t1_3 = e_ir2 + e_ir20 + w0;
                // double t2_3 = e_ir2 + e_ir21 + w1;
                // double t3_3 = e_ir2 + e_ir22 + w2;
                // double t4_3 = e_ir2 + e_ir23 + w3;
                // t_min0 = (t1_3 < t_min0) ? t1_3 : t_min0;
                // t_min1 = (t2_3 < t_min1) ? t2_3 : t_min1;
                // t_min2 = (t3_3 < t_min2) ? t3_3 : t_min2;
                // t_min3 = (t4_3 < t_min3) ? t4_3 : t_min3;

                // double t1_4 = e_ir3 + e_ir30 + w0;
                // double t2_4 = e_ir3 + e_ir31 + w1;
                // double t3_4 = e_ir3 + e_ir32 + w2;
                // double t4_4 = e_ir3 + e_ir33 + w3;
                // t_min0 = (t1_4 < t_min0) ? t1_4 : t_min0;
                // t_min1 = (t2_4 < t_min1) ? t2_4 : t_min1;
                // t_min2 = (t3_4 < t_min2) ? t3_4 : t_min2;
                // t_min3 = (t4_4 < t_min3) ? t4_4 : t_min3;
                
                // e[idx_ij+0] = t_min0;
                // e[idx_ij+1] = t_min1;
                // e[idx_ij+2] = t_min2;
                // e[idx_ij+3] = t_min3;
                // // XXX don't care about r for now

                // all the loads
                __m128d t_min0 = _mm_loadu_pd( &e[idx_ij+0] );
                __m128d t_min2 = _mm_loadu_pd( &e[idx_ij+2] );
                __m128d w0     = _mm_loadu_pd( &w[idx_ij+0] );
                __m128d w2     = _mm_loadu_pd( &w[idx_ij+2] );

                int ir_start_0 = idx_ib_x + STRIDE + j - ib;
                int ir_start_1 = ir_start_0 + STRIDE;
                int ir_start_2 = ir_start_1 + STRIDE;
                int ir_start_3 = ir_start_2 + STRIDE;
                __m128d e_ir00 = _mm_load_pd( &e[ir_start_0+0] );
                __m128d e_ir02 = _mm_load_pd( &e[ir_start_0+2] );
                __m128d e_ir10 = _mm_load_pd( &e[ir_start_1+0] );
                __m128d e_ir12 = _mm_load_pd( &e[ir_start_1+2] );
                __m128d e_ir20 = _mm_load_pd( &e[ir_start_2+0] );
                __m128d e_ir22 = _mm_load_pd( &e[ir_start_2+2] );
                __m128d e_ir30 = _mm_load_pd( &e[ir_start_3+0] );
                __m128d e_ir32 = _mm_load_pd( &e[ir_start_3+2] );

                // do the comps
                __m128d t00_1 = _mm_add_pd(v_e_ir0, w0);
                __m128d t02_1 = _mm_add_pd(v_e_ir0, w2);
                __m128d t00_2 = _mm_add_pd(t00_1, e_ir00);
                __m128d t02_2 = _mm_add_pd(t02_1, e_ir02);

                __m128d t10_1 = _mm_add_pd(v_e_ir1, w0);
                __m128d t12_1 = _mm_add_pd(v_e_ir1, w2);
                __m128d t10_2 = _mm_add_pd(t10_1, e_ir10);
                __m128d t12_2 = _mm_add_pd(t12_1, e_ir12);

                __m128d t20_1 = _mm_add_pd(v_e_ir2, w0);
                __m128d t22_1 = _mm_add_pd(v_e_ir2, w2);
                __m128d t20_2 = _mm_add_pd(t20_1, e_ir20);
                __m128d t22_2 = _mm_add_pd(t22_1, e_ir22);

                __m128d t30_1 = _mm_add_pd(v_e_ir3, w0);
                __m128d t32_1 = _mm_add_pd(v_e_ir3, w2);
                __m128d t30_2 = _mm_add_pd(t30_1, e_ir30);
                __m128d t32_2 = _mm_add_pd(t32_1, e_ir32);

                // mins
                __m128d t_min01 = _mm_blendv_pd(t_min0, t00_2, _mm_cmplt_pd(t00_2, t_min0));
                __m128d t_min21 = _mm_blendv_pd(t_min2, t02_2, _mm_cmplt_pd(t02_2, t_min2));
                __m128d t_min02 = _mm_blendv_pd(t_min01, t10_2, _mm_cmplt_pd(t10_2, t_min01));
                __m128d t_min22 = _mm_blendv_pd(t_min21, t12_2, _mm_cmplt_pd(t12_2, t_min21));
                __m128d t_min03 = _mm_blendv_pd(t_min02, t20_2, _mm_cmplt_pd(t20_2, t_min02));
                __m128d t_min23 = _mm_blendv_pd(t_min22, t22_2, _mm_cmplt_pd(t22_2, t_min22));
                __m128d t_min04 = _mm_blendv_pd(t_min03, t30_2, _mm_cmplt_pd(t30_2, t_min03));
                __m128d t_min24 = _mm_blendv_pd(t_min23, t32_2, _mm_cmplt_pd(t32_2, t_min23));

                // writeback
                _mm_storeu_pd( &e[idx_ij+0], t_min04);
                _mm_storeu_pd( &e[idx_ij+2], t_min24);
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

size_t bst_get_root_142_141_sse( void* _bst_obj, size_t i, size_t j )
{
    // [i,j], in table: [i-1, j]+1
    segments_t *mem = _bst_obj;
    size_t n = mem->n;
    int *root = mem->r;
    return (size_t) root[(i-1)*(n+1)+j]+1;
}

void bst_free_142_141_sse( void* _mem ) {
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

size_t bst_flops_142_141_sse( size_t n ) {
    size_t n3 = n*n*n;
    size_t n2 = n*n;
    return (size_t) ( n3/3.0 + 2*n2 + 5.0*n/3 );
}
