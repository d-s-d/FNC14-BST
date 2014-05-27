#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdarg.h>
#include <getopt.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#include "urange.h"
#include "variants.h"
#include "bst.h"
#include "perf/rdtsc.h"
#include "perfplot.h"

#define def_str(s) _def_str(s) // for some reason, need (s)
#define _def_str(s) #s

#define TS_FMT "%Y-%m-%d_%H-%M-%S"
#define FNAME_NAMES_SIZES "names_sizes.txt"

int nb1;

struct {
    // Test Setup
    urange_t N;
    char   **tests;  ///< NULL-terminated list of implementations to test
    unsigned int seed; ///< seed for random number generator
    int      calibrate;

    urange_t nb1_range; ///< range for first blocking parameter

    // Validation
    bst_impl_t* validation;
    double* valid_values;

    FILE   *fd;      ///< handle for names_sizes.txt
    // Determined at runtime
    double *p;
    double *q;
    struct perf_data *perf_data;
} config;

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

#define CYCLES_REQUIRED 1e8
int run_test(size_t n, bst_impl_t *impl, double ref, int verify)
{
    assert(impl);

    void *bst_data = impl->alloc(n);
    if (!bst_data) {
        RUN_TEST_ERR("Could not allocate memory.\n");
        return -1;
    }

    double e;
    int num_runs = 1;
    if (config.calibrate) {
        // warm up cache
        impl->compute(bst_data, config.p, config.q, n);

        // calibrate
        tsc_counter start, end;
        while (num_runs < (1<<14)) {
            CPUID(); RDTSC(start);
            for (int i=0; i<num_runs; ++i) {
                impl->compute(bst_data, config.p, config.q, n);
            }
            RDTSC(end); CPUID();

            double cycles = (double)(COUNTER_DIFF(end, start));
            if (cycles >= CYCLES_REQUIRED) break;

            num_runs *=2;
        }

        
        measurement_start();
        for (int i=0; i<num_runs; ++i) {
            e = impl->compute(bst_data, config.p, config.q, n);
        }
        measurement_stop(1);
    } else {
        measurement_start();
        e = impl->compute(bst_data, config.p, config.q, n);
        measurement_stop(1);
    }

#ifdef DEBUG
    printf("Cost is: %lf. Root is %d.\n", e, impl->root(bst_data, 1, n));
    traverse(1, n, 0, bst_data, impl->root);
#endif

    int pass = 0;

    if (verify) {
        if (fabs(e - ref) > 1e-5) {
            pass = -2;
            RUN_TEST_ERR("Wrong result: %lf (diff: %lf)\n", e, e-ref);
        } else {
            printf("Verified result: %.3lf (diff: %.3lf)\n", e, e-ref);
        }
    }

    impl->free(bst_data);
    return pass;
}

bst_impl_t* get_implementation(char* name) {
    bst_impl_t *impl = NULL;
    for (size_t i=0; i<impl_size && !impl; ++i) {
        if (strncmp(implementations[i].name, name, 100) == 0) {
            impl = &implementations[i];
        }
    }
    return impl;
}

// escape strings from start to stop
size_t escape_character( char *restrict trgt_str, const char *restrict str,
        char chr, size_t start, size_t stop ) {
    size_t str_len = strlen(str);
    size_t i, j = 0;
    // count occurences
    for( i = 0; i < str_len+1; i++ ) {
        if( (str[i] == chr) && (i>=start) && (i<=stop) ) trgt_str[j++] = '\\';
        trgt_str[j++] = str[i];
    }
    return i;
}

void run_configuration()
{
    size_t ENV_USERFLAGS_LEN =
        sizeof(def_str(M_ENV_USERFLAGS))/sizeof(def_str(M_ENV_USERFLAGS)[0]);
    size_t ENV_GITREV_LEN =
        sizeof(def_str(M_ENV_GITREV))/sizeof(def_str(M_ENV_GITREV)[0]);
    char str_env_userflags[2*ENV_USERFLAGS_LEN];
    char str_env_gitrev[2*ENV_GITREV_LEN];
    escape_character( str_env_userflags, def_str(M_ENV_USERFLAGS), '"', 
        1, ENV_USERFLAGS_LEN-3 );
    escape_character( str_env_gitrev, def_str(M_ENV_GITREV), '"',
        1, ENV_GITREV_LEN-3 );
    // some output log
    long custom_counters[]
        = { 0x530110, ///< FP_COMP_OPS_EXE.X87
            0x538010, ///< FP_COMP_OPS_EXE.SSE_SCALAR_DOUBLE
            0x531010, ///< FP_COMP_OPS_EXE.SSE_FP_PACKED_DOUBLE
            0x530211, ///< SIMD_FP_256.PACKED_DOUBLE
            0x530110, // whatever
            0x530110,
            0x530110,
            0x530110 };

    printf("starting measurement\n");
    measurement_init( custom_counters, 0, 0 );
    // debug output
    LOG("%zd implementations available:\n", impl_size);
    for (size_t i=0; i<impl_size; ++i) {
        LOG("  %s\n", implementations[i].name);
    }

    // run over all implementations in config.tests
    for (char **impl_name=config.tests; *impl_name; ++impl_name) {

        // See if implementation available
        LOG("Evaluating implementation '%s'.\n", *impl_name);

        bst_impl_t* impl = get_implementation(*impl_name);

        if (!impl) {
            LOG("ERROR: Implementation not found.\n");
            return ;
        }

        int validate = config.validation!=NULL;
        // precalculate validation results
        if( validate ) {
            bst_impl_t* val_impl = config.validation;
            int i = 0;
            config.valid_values = (double*) malloc(
                        urange_get_size(&config.N) * sizeof(double) );
            if( config.valid_values ) {
                for (size_t n=config.N.start; n<=config.N.stop;
                     n+=config.N.step) {
                    void* bst_data = val_impl->alloc(n);
                    if( !bst_data ) {
                        LOG("ERROR: could not allocate memory for bst_data."
                            "\n");
                        config.validation = NULL;
                        break;
                    }
                    config.valid_values[i++] =
                            val_impl->compute(bst_data, config.p, config.q, n);
                    val_impl->free(bst_data);
#ifdef DEBUG
                    printf("Calculated value %lf for input size %zu.\n",
                           config.valid_values[i-1], n);
#endif
                }
            } else {
                LOG("ERROR: could not allocate memory for valid_values.\n");
                config.validation = NULL;
            }
        }

        size_t final_i_len = strlen(impl->name)+11;
        char final_impl_name[final_i_len];
        int test_ret = 0;
        for( nb1=config.nb1_range.start; nb1<=config.nb1_range.stop
                && !test_ret; nb1+=config.nb1_range.step ) {
            snprintf( final_impl_name, final_i_len, "%s:%d", impl->name,
                nb1 );
            // Sweep through the tests
            LOG("NB1=%d.\n", nb1);
            int i = 0;
            for (size_t n=config.N.start; n<=config.N.stop && !test_ret;
                 n+=config.N.step) {
                LOG(" N = %zu\n", n);
                fflush(stdout);
                test_ret = run_test(n, impl, validate ? config.valid_values[i] :
                                                        0.0, validate);
                fprintf(config.fd, "%s %zu\n", impl->name, n);
                fflush(config.fd);
                i++;
            }
        }
        if( test_ret )
            LOG("ERROR: Test failed.\n");

        if( config.valid_values )
            free(config.valid_values);
    }
    measurement_end();
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
"  --seed <seed>\n"
"    sets the random number seed to <seed> (default: 42 ;)\n\n"
"  --validate <implementation>\n"
"    validate against <implementation>.\n\n"
"  --nb1 <range>\n"
"    set range for the first blocking parameter (default 1:1:1)\n\n"
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
#define LOGFN_LEN 40
    int ret;
    config.seed    = 42;

    config.nb1_range.start = 1;
    config.nb1_range.stop  = 1;
    config.nb1_range.step  = 1;

    int c;
    while (true) {
        static struct option long_options[] = {
            {"seed",      required_argument, 0, 's'},
            {"validate",  required_argument, 0, 'v'},
            {"calibrate", required_argument, 0, 'c'},
            {"nb1",       required_argument, 0, 'n'},
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

        case 's':
            config.seed = atoi(optarg);
            break;
        case 'v':
            config.validation = get_implementation(optarg);
            break;
        case 'c':
            config.calibrate = atoi(optarg);
            break;
        case 'n':
            if( read_urange( &config.nb1_range, optarg ) < 0 ) {
                ERROR_MSG("Error parsing argument: %s \n", optarg);
                print_usage_and_exit();
            }
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
    
    if( !(config.fd = fopen(FNAME_NAMES_SIZES, "w")) ) {
        ERROR_MSG("Error opening "FNAME_NAMES_SIZES".\n");
        return -1;
    }

    run_configuration();

    // perf_cleanup(config.perf_data);
    data_cleanup();

    fclose( config.fd );
    return 0;
}
