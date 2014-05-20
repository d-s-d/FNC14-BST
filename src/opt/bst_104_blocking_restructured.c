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

extern int nb1;

typedef struct {
    double* e;
    double* w;
    int* r;
    size_t n;
} segments_t;

void* bst_alloc_104_blocking_restructured( size_t n ) {
    segments_t* mem = (segments_t*) malloc( sizeof(segments_t) );
    size_t sz = (n+1)*(n+1);
    // XXX: for testing: calloc
    mem->e = (double*) calloc(1,  sz * sizeof(double) );
    mem->w = (double*) calloc(1,  sz * sizeof(double) );
    mem->r = (int*)    calloc(1,  sz * sizeof(int) );
    memset( mem->r, -1, sz * sizeof(int) );
    return mem;
}

inline void compute_min( int i, int j, int n, int r1, int r2, double* e,
        double* w, int* root ) {
    double t;
    int r;
    for( r = r1; r < r2; ++r ) {
        t = e[IDX(i,r)] + e[IDX(r+1,j)] + w[IDX(i,j)];
        if (t < e[IDX(i,j)]) {
            e[IDX(i,j)] = t;
            root[IDX(i,j)] = r;
        }
    }
}

double bst_compute_104_blocking_restructured( void*_bst_obj, double* p, double* q, size_t nn ) {
    segments_t* mem = (segments_t*) _bst_obj;
    int n;
    int i, l, j;
    double t, t_min, w_cur;
    int r_min;
    double* e = mem->e, *w = mem->w;
    int* root = mem->r;
    // initialization
    mem->n = nn;
    n = nn; // subtractions with n potentially negative. say hello to all the bugs
    for( i = 0; i < n+1; i++ ) {
        e[IDX(i,i)] = q[i];
        w[IDX(i,i)] = q[i];
    }

    int ib;

    // compute bottom right triangle nb1xnb1 (i.e. bottom nb1 rows)
    for (i = n-1; (i >= 0) && (i > (n-nb1)); --i) {
        for (j = i+1; j < n+1; ++j) {
            e[IDX(i,j)] = INFINITY;
            w[IDX(i,j)] = w[IDX(i,j-1)] + p[j-1] + q[j];
            compute_min( i, j, n, i, j, e, w, root );
       }
    }

    // compute the remaining rows
    // printf("Rest of rows: i=%d\n", i);
    for (; i >= 0; --i) {
        // First, the starting nb1 values in this row are computed.
        // This corresponds to completing the nb1xnb1 triangle right down from (i,i)
        for (j = i+1; j < (i+nb1); ++j) {
            e[IDX(i,j)] = INFINITY;
            w[IDX(i,j)] = w[IDX(i,j-1)] + p[j-1] + q[j];
            compute_min( i, j, n, i, j, e, w, root );
       }

        // Now we compute the rest of the row, but do updates in chunk of nb1's.
        // Since we now have the first nb1 values in this row, we can compute
        // the first nb1 iterations of the r-loop for all the remaining values
        // in this row.
        // (this needs to be seperated from the loop afterwards since we also do
        //  initialization to INFINITY))
        for (; j < (n+1); ++j) {
            e[IDX(i,j)] = INFINITY;
            w[IDX(i,j)] = w[IDX(i,j-1)] + p[j-1] + q[j];
            compute_min( i, j, n, i, i+nb1, e, w, root );
       }

        // We now continue to update the values in this row in chunks of nb1
        // as long as possible.
        for (ib = i+nb1; (ib+nb1) < (n+1); ib += nb1) {
            //printf("got in here for i=%d\n", i);

            // Again we start by finishing computing the next nb1 values of
            // row 'i'. The last values needed for that are from the triangle
            // in down right from (ib,ib).
            for (j = (ib+1); j < (ib+nb1); ++j) {
                compute_min( i, j, n, ib, j, e, w, root );
            }

            // Now, having nb1 new values in row 'i', we compute the next nb1
            // r-iterations for the remaining values in row 'i'.
            for (; j < (n+1); ++j) {
                compute_min( i, j, n, ib, ib + nb1, e, w, root );
            }
        }

        // There are less than nb1 elements remaining in row 'i'. The values
        // missing come from the triangle down right of (ib,ib)
        for (j = (ib+1); j < (n+1); ++j) {
            compute_min( i, j, n, ib, j, e, w, root );
       }
    }

    return e[IDX(0,n)];
}

size_t bst_get_root_104_blocking_restructured( void* _bst_obj, size_t i, size_t j )
{
    // [i,j], in table: [i-1, j]+1
    segments_t *mem = _bst_obj;
    size_t n = mem->n;
    int *root = mem->r;
    return (size_t) root[(i-1)*(n+1)+j]+1;
}

void bst_free_104_blocking_restructured( void* _mem ) {
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

size_t bst_flops_104_blocking_restructured( size_t n ) {
    size_t n3 = n*n*n;
    size_t n2 = n*n;
    return (n3 + 5*n2)/2 + 2*n;
}
