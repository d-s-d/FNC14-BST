#include <stdlib.h>
#include <stdio.h>

#include "perfmon_wrapper.h"

// NOTE: COST_MEASURE needs to be defined before including this file if
//       counters should be used (e.g. as compiler flag).
#include "cost_model.h"

COST_VARIABLES_HERE

extern void other_computation(); // defined in other.c

#define N 300         // could put this into main, but valgrind doesn't
static double A[N*N]; // like big data on the stack.
static double B[N*N];
static double C[N*N];

int main(int argc, char *argv[])
{
    int ret;

    // For a list of valid performance events, see
    //   libpfm/examples/showevtinfo -L
    // or
    //   Intel 64 and IA-32 Architectures Software Developer’s Manual,
    //     Volume 3B: System Programming Guide, Part 2
    
    // Initialize perfmon -------------------------------------------------
    char *events[] = {
        // on Haswell, can use at most 4, otherwise '0' results
        "PERF_COUNT_HW_CPU_CYCLES",
        "CACHE-REFERENCES",
        "CACHE-MISSES",
        "DTLB-LOAD-MISSES",
        NULL
    };

    struct perf_data *data;
    
    ret = perf_init(events, &data);
    if (ret < 0) {
        fprintf(stderr, "Could not initialize perfmon.\n");
        perf_cleanup(data);
        return 0;
    }

    // Example 1: Measure performance for MMM in this file. ---------------

    perf_start(data);

    for (int i=0; i<N; ++i) {
        for (int j=0; j<N; ++j) {
            for (int k=0; k<N; ++k) {
                COST_INC_ADD(1); COST_INC_MUL(1);
                C[i*N+j] = A[i*N+k] * B[k*N+j] + C[i*N+j];
            }
        }
    }

    ret = perf_stop(data);
    if (ret < 0) {
        fprintf(stderr, "Error measuring performance.\n");
        perf_cleanup(data);
        return 0;
    }

    // short cycles check
    perf_update_values(data);
    printf("Cycles after example 1: %llu\n", data[0].value);

    // Example 2: Measure external computation ----------------------------

    perf_start(data);    // if desired, do perf_reset() first
    other_computation();
    perf_stop(data);

    // Report results -----------------------------------------------------

    perf_update_values(data);

    printf("Overall results:\n");
    printf("Performance Counters:\n");
    for (struct perf_data *iter = data; iter->name; ++iter) {
        printf("  %-50s : %"PRId64"\n", iter->name, iter->value);
    }

    printf("Cost Model:\n");
    printf("  Add: %"PRI_COST"\n", COST_ADD);
    printf("  Mul: %"PRI_COST"\n", COST_MUL);
    printf("  Div: %"PRI_COST"\n", COST_DIV);
    printf("  Pow: %"PRI_COST"\n", COST_POW);
    printf("  Abs: %"PRI_COST"\n", COST_ABS);

    printf("Derived Results\n");
    double flops      = COST_ADD + COST_MUL + COST_DIV;
    double cycles     = data[0].value; // index as in 'events'
    double cache_load = data[1].value;
    double cache_miss = data[2].value;
    printf("  Flops           : %.0lf\n", flops);
    printf("  Performance     : %.3lf\n", flops/cycles);
    printf("  Cache Miss Rate : %.3lf\n", cache_miss/cache_load);

    // Cleanup ------------------------------------------------------------
    perf_cleanup(data);

    return 0;
}
