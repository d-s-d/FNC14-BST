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

void* bst_alloc_140_102_opt( size_t n ) {
    segments_t* mem = (segments_t*) malloc( sizeof(segments_t) );
    size_t sz = (n+1)*(n+1);
    // XXX: for testing: calloc
    mem->e = (double*) calloc(1,  sz * sizeof(double) );
    mem->w = (double*) calloc(1,  sz * sizeof(double) );
    mem->r = (int*)    calloc(1,  sz * sizeof(int) );
    memset( mem->r, -1, sz * sizeof(int) );
    return mem;
}

#define NB 4
double bst_compute_140_102_opt( void*_bst_obj, double* p, double* q, size_t nn ) {
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

size_t bst_get_root_140_102_opt( void* _bst_obj, size_t i, size_t j )
{
    // [i,j], in table: [i-1, j]+1
    segments_t *mem = _bst_obj;
    size_t n = mem->n;
    int *root = mem->r;
    return (size_t) root[(i-1)*(n+1)+j]+1;
}

void bst_free_140_102_opt( void* _mem ) {
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

size_t bst_flops_140_102_opt( size_t n ) {
    size_t n3 = n*n*n;
    size_t n2 = n*n;
    return (size_t) ( n3/3.0 + 2*n2 + 5.0*n/3 );
}
