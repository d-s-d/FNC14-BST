#include <stdlib.h>
#include <stdio.h>

#include "bst.h"

void traverse(size_t i, size_t j, size_t indent, void *_bst_obj)
{
    size_t root = bst_get_root(_bst_obj, i, j);

    for (size_t k=0; k<indent; ++k) printf(" ");
    printf("%-2zu\n", root);

    if (i != j) {
        if (i == root) {
            traverse(i+1, j, indent+3, _bst_obj);
        } else if (j == root) {
            traverse(i, j-1, indent+3, _bst_obj);
        } else {
            traverse(i, root-1, indent+3, _bst_obj);
            traverse(root+1, j, indent+3, _bst_obj);
        }
    }
}

int main(int argc, char *argv[])
{
    // size_t n = 10;
    // double p[n];
    // double q[n+1];

    size_t n = 5;
    double p[] = {     .15, .1 , .05, .1 , .2};
    double q[] = {.05, .1 , .05, .05, .05, .1};

    void *bst_data = bst_alloc(n);
    if (!bst_data) {
        fprintf(stderr, "Error allocating memroy\n");
        return -1;
    }

    double e = bst_compute(bst_data, p, q, n);

    printf("Cost is: %lf. Root is %d.\n", e, bst_get_root(bst_data, 1, n));
    traverse(1, n, 0, bst_data);

    bst_free(bst_data);

    return 0;
}
