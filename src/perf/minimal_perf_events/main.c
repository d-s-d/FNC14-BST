/**
 * This code shows how to run measurements using the linux 'perf events'
 * kernel interface. It is adopted from the Fastcode HW01.
 *
 * This code computes a simple MMM. If Sandy Bridge is used, D_MISSES or
 * D_FLOPS can be defined to get more raw PMU measurements. libpfm would
 * allow to use generic names rather than micro-architecture specific raw
 * events.
 *
 * 2014, Jeremia Baer <baerj@student.ethz.ch>
 */


#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <sys/ioctl.h>
#include <linux/perf_event.h>
#include <asm/unistd.h>
#include "measure.h"
#include "mmm.h"

#define NELEMS(x)  (sizeof(x) / sizeof(x[0]))

/**
 * for Intel Sandy Bridge
 *
 * Note: Raw events work on Sandy Bridge only. Haswell returns '0' results.
 * Note: If too many events are measured (M_FLOPS & M_MISSES),
 *       '0' values are returned for ALL events.
 */
event events [] = {
		{.fd=-1, .name="HW_CPU_CYCLES",                        .type=PERF_TYPE_HARDWARE, .config=PERF_COUNT_HW_CPU_CYCLES, .value=0},
#ifdef M_FLOPS
		{.fd=-1, .name="FP_COMP_OPS_EXE:X87",                  .type=PERF_TYPE_RAW, .config=0x530110, .value=0},
		{.fd=-1, .name="FP_COMP_OPS_EXE:SSE_FP_PACKED_DOUBLE", .type=PERF_TYPE_RAW, .config=0x531010, .value=0},
		{.fd=-1, .name="FP_COMP_OPS_EXE:SSE_SCALAR_DOUBLE",    .type=PERF_TYPE_RAW, .config=0x538010, .value=0},
		{.fd=-1, .name="SIMD_FP_256:PACKED_DOUBLE",            .type=PERF_TYPE_RAW, .config=0x530211, .value=0},
#endif // M_FLOPS
#ifdef M_MISSES
		{.fd=-1, .name="DTLB_LOAD_MISSES:MISS_CAUSES_A_WALK",  .type=PERF_TYPE_RAW, .config=0x530108, .value=0},
		{.fd=-1, .name="DTLB_STORE_MISSES:MISS_CAUSES_A_WALK", .type=PERF_TYPE_RAW, .config=0x530149, .value=0},
		{.fd=-1, .name="MEM_LOAD_RETIRED:L1_HIT",              .type=PERF_TYPE_RAW, .config=0x5301d1, .value=0},
		{.fd=-1, .name="MEM_LOAD_RETIRED:L2_HIT",              .type=PERF_TYPE_RAW, .config=0x5302d1, .value=0},
		{.fd=-1, .name="MEM_LOAD_RETIRED:L3_HIT",              .type=PERF_TYPE_RAW, .config=0x5304d1, .value=0},
		{.fd=-1, .name="MEM_LOAD_RETIRED:L3_MISS",             .type=PERF_TYPE_RAW, .config=0x5320d1, .value=0},
#endif // M_MISSES
};

int main(int argc, char **argv) {

	int i; double flops; double cycles;

	// if (argc!=4) {printf("usage: mmm <m> <k> <n>\n"); return -1;}
	// m = atoi(argv[1]);
	// k = atoi(argv[2]);
	// n = atoi(argv[3]);
	// printf("m = %d k = %d n = %d\n", m, k, n);
	m = n = k = 1000;

	A = (double *) malloc (m * k * sizeof(double));
	B = (double *) malloc (k * n * sizeof(double));
	C = (double *) malloc (m * n * sizeof(double));

	for(i = 0; i < m * k; i++) A[i]=i;
	for(i = 0; i < k * n; i++) B[i]=i;

	compute(); // warm cache
	compute(); // warm cache

	measurements_init (events, NELEMS(events));

	measurement_start (events);
	compute();
	measurement_stop (events);

	measurements_finish (events, NELEMS(events));

	cycles = (double) events[0].value;
	flops  = (double) (events[2].value + events[3].value);

	printf("Obtained performance: %lf flops/cycles\n", flops/cycles);

}
