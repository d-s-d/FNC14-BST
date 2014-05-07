#include<bst.h>

typedef struct {
    double* e;
    double* w;
    int* r;
} segments_t;

bool simple_bst_test() {
    void* bst_obj;
    double p[] = {     .15, .1 , .05, .1 , .2};
    double q[] = {.05, .1 , .05, .05, .05, .1};
    double ref_res = 2.75;

    size_t n = sizeof(p)/sizeof(p[0]);
    bst_obj = bst_alloc( n );
    if( !bst_obj ) {
        printf("ERROR: bst_alloc returned NULL.\n");
        return false;
    }

    double res = bst_compute( bst_obj, p, q, n );
    if( abs( res - ref_res ) > 0.000001 ) {
        printf("ERROR: Result sould be %lf, but is %lf.\n", ref_res, res);
        return false;
    } else {
        printf("Result is: %lf\n", res);
        printf("Root structure: \n");
        size_t i,j;
        int* root = ((segments_t*)bst_obj)->r;
        for( i = 0; i < n+1; i++ ) {
            for( j = 0; j < n+1; j++ ) {
                printf("%.3d ", root[i*(n+1)+j]);
            }
            printf("\n");
        }
    }
    return true;
}
