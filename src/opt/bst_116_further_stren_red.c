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

void* bst_alloc_116_further_stren_red( size_t n ) {
    segments_t* mem = (segments_t*) malloc( sizeof(segments_t) );
    size_t sz = (n+1)*(n+1);
    size_t sz2 = (n+1)*(n+2)/2;
    // XXX: for testing: calloc
    mem->e = (double*) calloc(1,  sz2 * sizeof(double) );
    mem->w = (double*) calloc(1,  sz2 * sizeof(double) );
    mem->r = (int*)    calloc(1,  sz2 * sizeof(int) );
    mem->n = n;
    memset( mem->r, -1, sz2 * sizeof(int) );
    return mem;
}

double bst_compute_116_further_stren_red( void*_bst_obj, double* p, double* q, size_t nn ) {
    segments_t* mem = (segments_t*) _bst_obj;
    int n, i, r, j;
    double t;
    double* e = mem->e, *w = mem->w;
    int* root = mem->r;
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
        for (j = i+1; j < n+1; ++j) {
            e[idx2] = INFINITY;
            w[idx2] = w[idx2-1] + p[j-1] + q[j];
            idx2++;
        }
        idx3 = idx1; 
        for (r = i; r < n; ++r) {
            // idx2 = IDX(r+1, r+1);
            idx1 = idx3;
            double e_tmp = e[idx1++];
            for (j = r+1; j < n+1; ++j) {
                printf("idx1: %d, idx2: %d\n", idx1, idx2);
                t = e_tmp + e[idx2++] + w[idx1];
                if (t < e[idx1]) {
                    e[idx1] = t;
                    root[idx1] = r;
                }
                idx1++;
            }
            idx3++;
        }
    }


    return e[IDX(0,n)];
}

size_t bst_get_root_116_further_stren_red( void* _bst_obj, size_t i, size_t j )
{
    // [i,j], in table: [i-1, j]+1
    segments_t *mem = _bst_obj;
    size_t n = mem->n;
    int *root = mem->r;
    return (size_t) root[(i-1)*(n+1)+j]+1;
}

void bst_free_116_further_stren_red( void* _mem ) {
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

size_t bst_flops_116_further_stren_red( size_t n ) {
    double n3 = n*n*n;
    double n2 = n*n;
    return (size_t) ( n3/3.0 + 2*n2 + 5.0*n/3 );
}
