#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "bst.h"

struct {
    char   *input_file; ///< name of input file for data
    char   *log_name;   ///< name of logfile

    // derived values
    double *p;
    double *q;
    size_t max_n;     ///< max N supported by p, q
    FILE   log_fd;    ///< handle for logfile
    int    log_valid; ///< true iff logfile open for writes
} config;

bst_impl_t reference = {
    .name    = "reference",
    .alloc   = bst_alloc,
    .compute = bst_compute,
    .root    = bst_get_root,
    .free    = bst_free
};

void traverse(size_t i, size_t j, size_t indent, void *_bst_obj,
        bst_get_root_fn root_fn)
{
    // printf("DBG: i=%zu, j=%zu\n", i, j);
    size_t root = root_fn(_bst_obj, i, j);

    for (size_t k=0; k<indent; ++k) printf(" ");
    printf("%-2zu\n", root);

    if (i != j) {
        if (i == root) {
            traverse(i+1, j, indent+3, _bst_obj, root_fn);
        } else if (j == root) {
            traverse(i, j-1, indent+3, _bst_obj, root_fn);
        } else {
            traverse(i, root-1, indent+3, _bst_obj, root_fn);
            traverse(root+1, j, indent+3, _bst_obj, root_fn);
        }
    }
}

#define RUN_TEST_ERR(...)                                                \
    fprintf(stderr, "ERROR[run_test, impl=%s, n=%zu]: ", impl->name, n); \
    fprintf(stderr, __VA_ARGS__);

int run_test(size_t n, bst_impl_t *impl,
        double ref, int verify)
{
    assert(impl);

    void *bst_data = impl->alloc(n);
    if (!bst_data) {
        RUN_TEST_ERR("Could not allocate memory.\n");
        return -1;
    }

    double e = impl->compute(bst_data, config.p, config.q, n);

    printf("Cost is: %lf. Root is %d.\n", e, impl->root(bst_data, 1, n));
    traverse(1, n, 0, bst_data, impl->root);

    int pass = 0;

    if (verify) {
        if (abs(e - ref) > 1e-5) {
            pass = -2;
            RUN_TEST_ERR("Wrong result: %lf (diff: %lf)\n", e, e-ref);
        }
    }

    impl->free(bst_data);
    return pass;
}

void sweep_tests(size_t from, size_t to, size_t step)
{
    for (size_t n=from; n<=to; n+=step) {
        if (run_test(n, &reference, 2.75, 1) < 0) {
            printf("Test failed.\n");
        } else {
            printf("Test passed.\n");
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

    config.max_n = n;
    config.p     = p;
    config.q     = q;

    sweep_tests(1, 5, 1);

    return 0;
}
