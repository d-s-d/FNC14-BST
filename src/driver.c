#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdarg.h>
#include <getopt.h>
#include <stdbool.h>

#include "urange.h"
#include "variants.h"
#include "bst.h"
#include "perf/perfmon_wrapper.h"
#include "driver_log.h"

#define def_str(s) _def_str((s)) // for some reason, need (s)
#define _def_str(s) #s

struct {
    // Test Setup
    urange_t N;
    char   **tests;  ///< NULL-terminated list of implementations to test
    unsigned int seed; ///< seed for random number generator

    // Logging Setup
    char   *logfile; ///< name of logfile

    // Determined at runtime
    double *p;
    double *q;
    FILE   *fd;      ///< handle for logfile
    struct perf_data *perf_data;
} config;

bst_impl_t implementations[] = {
    {
        .name    = "ref/bst_ref.c",
        .alloc   = bst_alloc,
        .compute = bst_compute,
        .root    = bst_get_root,
        .free    = bst_free
    }, {
        .name    = "ref/bst_ref.c_dummy",
        .alloc   = bst_alloc,
        .compute = bst_compute,
        .root    = bst_get_root,
        .free    = bst_free
    }, {
        .name    = "opt/bst_001_transposed.c",
        .alloc   = bst_alloc_001_transposed,
        .compute = bst_compute_001_transposed,
        .root    = bst_get_root_001_transposed,
        .free    = bst_free_001_transposed
    }
};
#define impl_size (sizeof(implementations)/sizeof(bst_impl_t))

#define LOG(...) printf(__VA_ARGS__)
#define ERROR_MSG(...) fprintf(stderr, __VA_ARGS__)

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

    perf_reset(config.perf_data);
    perf_start(config.perf_data);
    double e = impl->compute(bst_data, config.p, config.q, n);
    perf_stop(config.perf_data);

#ifdef DEBUG
    printf("Cost is: %lf. Root is %d.\n", e, impl->root(bst_data, 1, n));
    traverse(1, n, 0, bst_data, impl->root);
#endif

    int pass = 0;

    if (verify) {
        if (abs(e - ref) > 1e-5) {
            pass = -2;
            RUN_TEST_ERR("Wrong result: %lf (diff: %lf)\n", e, e-ref);
        }
    }

    // log performance
    perf_update_values(config.perf_data);
    log_idouble("cycles",           config.perf_data[0].value);
    log_idouble("cache-references", config.perf_data[1].value);
    log_idouble("cache-misses",     config.perf_data[2].value);

    impl->free(bst_data);
    return pass;
}

void run_configuration()
{
    // some output log
    log_size("from", config.N.start);
    log_size("to",   config.N.stop);
    log_size("step", config.N.step);
    log_fmt("seed", "%u", config.seed);
    log_str("userflags", def_str(M_ENV_USERFLAGS));
    log_str("git-revision", def_str(M_ENV_GITREV));

    log_struct("runs");

    // debug output
    LOG("%zd implementations available:\n", impl_size);
    for (size_t i=0; i<impl_size; ++i) {
        LOG("  %s\n", implementations[i].name);
    }

    // run over all implementations in config.tests
    for (char **impl_name=config.tests; *impl_name; ++impl_name) {

        // See if implementation available
        LOG("Evaluating implementation '%s'.\n", *impl_name);

        bst_impl_t *impl;
        int found = 0;
        for (size_t i=0; i<impl_size && !found; ++i) {
            if (strncmp(implementations[i].name, *impl_name, 100) == 0) {
                found = 1;
                impl = &implementations[i];
            }
        }

        if (!found) {
            LOG("ERROR: Implementation not found.\n");
            return ;
        }

        // Sweep through the tests
        log_array(impl->name);
        for (size_t n=config.N.start; n<=config.N.stop; n+=config.N.step) {
            LOG("N = %zu\n", n);
            log_struct(NULL);
            log_int("N", n);
            if (run_test(n, impl, 2.75, 0) < 0) {
                LOG("ERROR: Test failed.\n");
                return ;
            }
            log_struct_end();
        }
        log_array_end();

    }

    log_struct_end();
}

int data_setup()
{
    size_t n = config.N.stop;
    config.p = NULL;
    config.q = NULL;

    // allocate memory
    config.p = malloc(n*sizeof(double));
    if (!config.p) {
        return 0;
    }

    config.q = malloc((n+1)*sizeof(double));
    if (!config.q) {
        free(config.p);
        return 0;
    }

    // generate data
    double total = 0.0;
    srand(config.seed);

    for (size_t i=0; i<n; ++i) {
        double v = rand();
        config.p[i] = v;
        total += v;
    }
    for (size_t i=0; i<(n+1); ++i) {
        double v = rand();
        config.q[i] = v;
        total += v;
    }

    // normalize data
    for (size_t i=0; i<n; ++i)     config.p[i] /= total;
    for (size_t i=0; i<(n+1); ++i) config.q[i] /= total;

#ifdef DEBUG
    for (size_t i=0; i<n; ++i)     printf("%lf ", config.p[i]);
    for (size_t i=0; i<(n+1); ++i) printf("%lf ", config.q[i]);
#endif

    return 1;
}

void data_cleanup()
{
    if (config.p) free(config.p);
    if (config.q) free(config.q);
}

const char* usage_str =
"\nUSAGE:\n"
"./bst_driver [options] <input_sizes> <implementations>\n"
"  options:\n\n"
"  --logfile <path>\n"
"    writes log messages to <path> (default: test.log)\n\n"
"  --seed <seed>\n"
"    sets the random number seed to <seed> (default: 42 ;)\n\n"
" <input_sizes>: matlab-like range definition (start:step:stop)\n"
"\n"
" Example:"
" $ ./bst_driver 10:10:100 ref/bst_ref.c\n\n"
"   Executes the reference implementation with inputs sizes 10,20,..,100\n";

void print_usage_and_exit() {
    printf("%s", usage_str);
    exit(1);
}

int main(int argc, char *argv[])
{
    int ret;

    // set defaults
    config.logfile = "test.log";
    config.seed    = 42;

    int c;
    while (true) {
        static struct option long_options[] = {
            {"logfile", required_argument, 0, 'l'},
            {"seed",    required_argument, 0, 's'},
            {0,0,0,0}
        };

        int option_index = 0;
        c = getopt_long(argc, argv, "t", long_options, &option_index);
        if (c == -1) { // -1 indicates end of options reached
            break;
        }
        switch (c) {
        case 0:
            // the long option with name long_options[option_index].name is
            // found
            ERROR_MSG("getopt error on long option %s\n",
                   long_options[option_index].name);
            break;

        case 'l':
            config.logfile = optarg;
            break;
        case 's':
            config.seed = atoi(optarg);
            break;
        case '?':
            printf("getopt: error on character %c\n", optopt);
            break;
        default:
            printf("getopt: general error\n");
            abort();
        }
    }

    int argc_remain = argc - optind;
    int tmp_optind = optind;
    // check whether enough input arguments were provided.
    if( argc_remain < 2 ) {
        ERROR_MSG("Too few arguments.\n");
        print_usage_and_exit();
    }

    // read the test range from the cli
    if( read_urange( &config.N, argv[tmp_optind++] ) < 0 ) {
        ERROR_MSG("Error parsing argument: %s\n", argv[1] );
        print_usage_and_exit();
    }

    config.tests = &argv[tmp_optind];

    ret = data_setup();
    if (!ret) {
        ERROR_MSG("Error allocating data.\n");
        return -1;
    }

    ret = log_setup(config.logfile);
    if (!ret) {
        ERROR_MSG("Error opening log.\n");
        return -1;
    }

    // log_fmt("some crazy double", "%08.1lf", 2.0);

    // Initialize perfmon -------------------------------------------------
    char *events[] = {
        // on Haswell, can use at most 4, otherwise '0' results
        "PERF_COUNT_HW_CPU_CYCLES",
        "CACHE-REFERENCES",
        "CACHE-MISSES",
        "DTLB-LOAD-MISSES",
        NULL
    };

    ret = perf_init(events, &config.perf_data);
    if (ret < 0) {
        ERROR_MSG("Could not initialize perfmon.\n");
        perf_cleanup(config.perf_data);
        return 0;
    }

    run_configuration();

    log_cleanup();
    perf_cleanup(config.perf_data);
    data_cleanup();

    return 0;
}
