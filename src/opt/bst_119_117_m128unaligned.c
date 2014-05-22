#include<bst.h>
#include<math.h>
#include<string.h>
#include<x86intrin.h>

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
 * This version swaps the r and j loop form 100_ref_bottomup.c
 */

#define STRIDE (n+1)
#define IDX(i,j) ((n+1)*(n+2)/2 - (n-(i)+1)*(n-(i)+2)/2 + (j) - (i))

typedef struct {
    double* e;
    double* w;
    int* r;
    size_t n;
} segments_t;

void* bst_alloc_119_117_m128unaligned( size_t n ) {
    segments_t* mem = (segments_t*) malloc( sizeof(segments_t) );
    size_t sz2 = (n+1)*(n+2)/2;
    // XXX: for testing: calloc
    mem->e = (double*) calloc(1,  sz2 * sizeof(double) );
    mem->w = (double*) calloc(1,  sz2 * sizeof(double) );
    mem->r = (int*)    calloc(1,  sz2 * sizeof(int) );
    mem->n = n;
    memset( mem->r, -1, sz2 * sizeof(int) );
    return mem;
}

double bst_compute_119_117_m128unaligned( void*_bst_obj, double* p, double* q, size_t nn ) {
    segments_t* mem = (segments_t*) _bst_obj;
    int n, i, r, l_end, j, l_end_pre;
    double t, e_tmp;
    double* e = mem->e, *w = mem->w;
    int* root = mem->r;
    __m128d v_tmp;
    __m128d v00, v01, v02;
    __m128d v10, v11, v12;
    __m128i v_cur_roots, v_new_roots, v_old_roots, v_root_inc;
    __m128 v_rootmask;
    // initialization
    // mem->n = nn;
    n = nn; // subtractions with n potentially negative. say hello to all the bugs

    int idx1, idx2, idx3;
    
    idx1 = IDX(n,n);
    e[idx1] = q[n];
    idx1++;
    for (i = n-1; i >= 0; --i) {
        idx1 -= 2*(n-i)+1;
        idx2 = idx1 + 1;
        e[idx1] = q[i];
        w[idx1] = q[i];
        for (j = i+1; j < n+1; ++j,++idx2) {
            e[idx2] = INFINITY;
            w[idx2] = w[idx2-1] + p[j-1] + q[j];
        }
        idx3 = idx1; 
        for (r = i; r < n; ++r) {
            // idx2 = IDX(r+1, r+1);
            idx1 = idx3;
            l_end = idx2 + (n-r);
            // l_end points to the first entry after the current row
            e_tmp = e[idx1++];
            // calculate until a multiple of 8 doubles is left
            // 8 = 4 * 2 128-bit vectors
            l_end_pre = idx2 + ((n-r)&3);
            for( ; (idx2 < l_end_pre) && (idx2 < l_end); ++idx2 ) {
                t = e_tmp + e[idx2] + w[idx1];
                if (t < e[idx1]) {
                    e[idx1] = t;
                    root[idx1] = r;
                }
                idx1++;
            }
            
            v_tmp = _mm_set_pd( e_tmp, e_tmp );
            // execute the shit for 4 vectors of size 2
            v_cur_roots = _mm_set_epi32(r, r+1, r+2, r+3);
            v_root_inc = _mm_set_epi32(4,4,4,4);
            for( ; idx2 < l_end; idx2 += 4 ) {
                v00 = _mm_loadu_pd( &e[idx2  ] );
                v10 = _mm_loadu_pd( &e[idx2+2] );

                v01 = _mm_loadu_pd( &w[idx1  ] );
                v11 = _mm_loadu_pd( &w[idx1+2] );

                v01 = _mm_add_pd( v01, v_tmp ); // supoptimal for raw-dependency
                v11 = _mm_add_pd( v11, v_tmp );

                v01 = _mm_add_pd( v01, v00 );
                v11 = _mm_add_pd( v11, v10 );

                v00 = _mm_loadu_pd( &e[idx1  ] );
                v10 = _mm_loadu_pd( &e[idx1+2] );

                v02 = _mm_cmplt_pd( v01, v00 );
                v12 = _mm_cmplt_pd( v11, v10 );

                v00 = _mm_or_pd( _mm_and_pd( v02, v01 ), _mm_andnot_pd( v02, v00 ));
                v10 = _mm_or_pd( _mm_and_pd( v12, v11 ), _mm_andnot_pd( v12, v10 ));

                _mm_storeu_pd( &e[idx1  ], v00 );
                _mm_storeu_pd( &e[idx1+2], v10 );

                v_rootmask = _mm_shuffle_ps(
                        _mm_castpd_ps( v02 ),
                        _mm_castpd_ps( v12 ),
                        _MM_SHUFFLE(0,2,0,2) );

                v_old_roots = _mm_lddqu_si128( &root[idx1] );
                v_new_roots = _mm_or_si128(
                        _mm_and_si128(    v_cur_roots, 
                            _mm_castps_si128( v_rootmask ) ),
                        _mm_andnot_si128( v_old_roots,
                            _mm_castps_si128( v_rootmask ) )
                        );

                _mm_storeu_si128( &root[idx1], v_new_roots );

                v_cur_roots = _mm_add_epi32( v_cur_roots, v_root_inc );

                idx1 += 4;
            }
            idx3++;
        }
    }


    return e[IDX(0,n)];
}

size_t bst_get_root_119_117_m128unaligned( void* _bst_obj, size_t i, size_t j )
{
    // [i,j], in table: [i-1, j]+1
    segments_t *mem = _bst_obj;
    size_t n = mem->n;
    int *root = mem->r;
    return (size_t) root[(i-1)*(n+1)+j]+1;
}

void bst_free_119_117_m128unaligned( void* _mem ) {
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

    free( mem->e );
    free( mem->w );
    free( mem->r );
    free( mem );
}

size_t bst_flops_119_117_m128unaligned( size_t n ) {
    double n3 = n*n*n;
    double n2 = n*n;
    return (size_t) ( n3/3.0 + 2*n2 + 5.0*n/3 );
}
