#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdarg.h>

#include "bst.h"
#include "perf/perfmon_wrapper.h"

#define def_str(s) _def_str(s)
#define _def_str(s) #s

struct {
    // Test Setup
    size_t from;     ///< minimum N
    size_t to;       ///< maximum N
    size_t step;     ///< stepsize for N
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
    }
};
#define impl_size (sizeof(implementations)/sizeof(bst_impl_t))

static int log_indent  = 0;
#define log_do_indent() \
    {for (int i=0; i<log_indent; ++i) fprintf(config.fd, "  ");}
int log_setup()
{
    config.fd = NULL;

    if (config.logfile) {
        config.fd = fopen(config.logfile, "w");
    }

    if (config.fd) {
        fprintf(config.fd, "{\n");
        log_indent++;
    } else {

    }

    return (config.fd != NULL);
}

void log_cleanup()
{
    if (config.fd) {
        fprintf(config.fd, "}\n");
        fclose(config.fd);
    }
}

void log_fmt(char *name, char *fmt, ...)
{
    if (config.fd) {
        log_do_indent();
        fprintf(config.fd, "'%s' : ", name);

        va_list argptr;
        va_start(argptr, fmt);
        vfprintf(config.fd, fmt, argptr);
        va_end(argptr);

        fprintf(config.fd, "\n");
    }
}

void log_int(char *name, int v)
{
    log_fmt(name, "%d", v);
}

void log_idouble(char *name, double v)
{
    log_fmt(name, "%.0lf", v);
}

void log_size(char *name, size_t v)
{
    log_fmt(name, "%zu", v);
}

void log_str(char *name, char *str)
{
    log_fmt(name, "'%s'", str);
}

void log_part(char *name, char delim)
{
    if (config.fd) {
        log_do_indent();
        if (name) {
            fprintf(config.fd, "'%s' : %c\n", name, delim);
        } else {
            fprintf(config.fd, "%c\n", delim);
        }
        log_indent++;
    }
}

void log_part_end(char delim)
{
    if (config.fd) {
        log_indent--;
        log_do_indent();
        fprintf(config.fd, "%c\n", delim);
    }
}

void log_struct(char *name)
{
    log_part(name, '{');
}

void log_struct_end()
{
    log_part_end('}');
}

void log_array(char *name)
{
    log_part(name, '[');
}

void log_array_end()
{
    log_part_end(']');
}

#define LOG(...) printf(__VA_ARGS__)

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
    log_size("from", config.from);
    log_size("to",   config.to);
    log_size("step", config.step);
    log_fmt("seed", "%u", config.seed);
    log_str("userflags", def_str(M_ENV_USERFLAGS));
    log_str("git-revision", def_str(M_ENV_GITREV));

    log_struct("runs");

    // debug output
    LOG("%d implementations available:\n", impl_size);
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
        for (size_t n=config.from; n<=config.to; n+=config.step) {
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
    size_t n = config.to;
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

int main(int argc, char *argv[])
{
    int ret;

    // size_t n = 5;
    // double p[] = {     .15, .1 , .05, .1 , .2};
    // double q[] = {.05, .1 , .05, .05, .05, .1};

    // config.p     = p;
    // config.q     = q;

    // config
    char *tests[] = {"ref/bst_ref.c", "ref/bst_ref.c_dummy", "ref/bst_ref.c", NULL};
    config.tests = tests;

    config.from    = 20;
    config.to      = 100;
    config.step    = 20;
    config.logfile = "test.log";
    config.seed    = 42;

    ret = data_setup();
    if (!ret) {
        printf("Error allocating data.\n");
        return -1;
    }

    ret = log_setup();
    if (!ret) {
        printf("Error opening log.\n");
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
        fprintf(stderr, "Could not initialize perfmon.\n");
        perf_cleanup(config.perf_data);
        return 0;
    }

    run_configuration();

    log_cleanup();
    perf_cleanup(config.perf_data);
    data_cleanup();

    return 0;
}
