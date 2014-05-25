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
 * This is an extension of bst_117. It implements the padding semantics that
 * enable proper alignment of vectors.
 *
 */

#define STRIDE (n+1)
#define IDX(i,j) ((n+1)*(n+2)/2 - (n-(i)+1)*(n-(i)+2)/2 + (j) - (i))

typedef struct {
    double* e;
    double* w;
    int* r;
    size_t n;
    size_t e_sz;
} segments_t;

inline size_t compute_size( size_t n ) {
    size_t m = (n+1)%4;
    size_t k = (3-m)*(3-m+1)/2;
    size_t rempad = 6 - k;
    return (n+1)*(n+2)/2 + ((n+1)/4)*6 + rempad;
}

void* bst_alloc_126_m256_padding( size_t n ) {
    segments_t* mem = (segments_t*) malloc( sizeof(segments_t) );
    mem->e_sz = compute_size( n );
    // XXX: for testing: calloc
    mem->e = (double*) calloc(1, mem->e_sz * sizeof(double) );
    mem->w = (double*) calloc(1, mem->e_sz * sizeof(double) );
    mem->r = (int*)    calloc(1, mem->e_sz * sizeof(int) );
    mem->n = n;
    memset( mem->r, -1, mem->e_sz * sizeof(int) );
    return mem;
}

double bst_compute_126_m256_padding( void*_bst_obj, double* p, double* q, size_t nn ) {
    segments_t* mem = (segments_t*) _bst_obj;
    int n, i, r, l_end, j;
    double t, e_tmp;
    double* e = mem->e, *w = mem->w;
    int* root = mem->r;
    // initialization
    // mem->n = nn;
    n = nn; // subtractions with n potentially negative. say hello to all the bugs

    int idx1, idx2, idx3, pad, pad_r;
    
    idx1 = ((int) mem->e_sz) - 1;
    // the conventio is that iteration i, idx1 points to the first element of line i+1
    e[idx1++] = q[n];
    
    // pad contains the padding for row i+1
    // for row n it's always 3
    pad = 3;
    for (i = n-1; i >= 0; --i) {
        idx1 -= 2*(n-i)+1 + pad;
        idx2 = idx1 + 1;
        e[idx1] = q[i];
        w[idx1] = q[i];
        for (j = i+1; j < n+1; ++j,++idx2) {
            e[idx2] = INFINITY;
            w[idx2] = w[idx2-1] + p[j-1] + q[j];
        }
        // idx2 now points to the beginning of the next line.
        idx2 += pad; // padding of line i+1

        idx3 = idx1;
        pad_r = pad; // padding of line r
        for (r = i; r < n; ++r) {
            pad_r = (pad_r+1)&3; // padding of line r+1
            // idx2 = IDX(r+1, r+1);
            idx1 = idx3;
            l_end = idx2 + (n-r);
            e_tmp = e[idx1++];
            for( ; idx2 < l_end; ++idx2 ) {
                // printf("idx1: %d, idx2: %d\n", idx1, idx2);
                t = e_tmp + e[idx2] + w[idx1];
                if (t < e[idx1]) {
                    e[idx1] = t;
                    root[idx1] = r;
                }
                idx1++;
            }
            idx2 += pad_r;
            idx3++;
        }
        pad = (pad-1)&3;
        // pad is 3,2,1,0,3,2,1,0, etc.
    }

    // the index of the last item of the first row is ((n/4)+1)*4-1, due to the padding
    // if n is even, the total number of entries in the first
    // row of the table is odd, so we need padding
    return e[ ((n/4)+1)*4 - 1 ];
}

size_t bst_get_root_126_m256_padding( void* _bst_obj, size_t i, size_t j )
{
    // [i,j], in table: [i-1, j]+1
    segments_t *mem = _bst_obj;
    size_t n = mem->n;
    int *root = mem->r;
    return (size_t) root[(i-1)*(n+1)+j]+1;
}

void bst_free_126_m256_padding( void* _mem ) {
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

size_t bst_flops_126_m256_padding( size_t n ) {
    double n3 = n*n*n;
    double n2 = n*n;
    return (size_t) ( n3/3.0 + 2*n2 + 5.0*n/3 );
}
