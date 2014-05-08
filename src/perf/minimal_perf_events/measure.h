#ifndef MEASURE_H_
#define MEASURE_H_

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <sys/ioctl.h>
#include <linux/perf_event.h>
#include <asm/unistd.h>
#include <errno.h>

typedef struct event_t {
	int fd;
	__u64 type;
	__u64 config;
	__u64 value;
	const char * const name;
} event;

inline void print_perf_event_error(int fd);

inline void measurements_init (event * events, int count)
{
	int i;
	for (i = 0; i < count; ++i) {
		struct perf_event_attr pe;
		memset(&pe, 0, sizeof(struct perf_event_attr));
		pe.type = events[i].type;
		pe.size = sizeof(struct perf_event_attr);
		pe.config = events[i].config;
		pe.disabled = 1;
		pe.exclude_kernel = 1;
		pe.exclude_hv = 1;
		pe.read_format = PERF_FORMAT_GROUP | PERF_FORMAT_ID;
		events[i].fd = syscall(__NR_perf_event_open, &pe, 0, -1, events->fd, 0);
		if (events[i].fd == -1) {
			fprintf(stderr, "Error opening event %llx\n", pe.config);
			print_perf_event_error(events[i].fd);
			exit(EXIT_FAILURE);
		} else {
			fprintf(stdout, "Event 0x%06x initialization success\n", (unsigned int) events[i].config);
		}
	}
	fprintf(stdout, "Monitoring: %d events\n\n", count);
}

inline void measurement_start (event * events) {
	ioctl(events->fd, PERF_EVENT_IOC_RESET, PERF_IOC_FLAG_GROUP);
	ioctl(events->fd, PERF_EVENT_IOC_ENABLE, PERF_IOC_FLAG_GROUP);
}

inline void measurement_stop (event * events) {
	ioctl(events->fd, PERF_EVENT_IOC_DISABLE, PERF_IOC_FLAG_GROUP);
}

inline void measurements_finish (event * events, int count) {
	int i; int size = (1 + count * 2) * sizeof(__u64);
	__u64 * v = (__u64 *) malloc (size);
	read (events->fd, v, size);
	printf("Events report:\n");
	printf("-------------------\n");
	for (i = 0; i < count; ++i) {
		close(events[i].fd);
		events[i].value = v[1 + 2*i];
		printf("[ 0x%06x | %-40s ] %lld\n", (unsigned int) events[i].config, events[i].name, events[i].value);
	}
	printf("\n");
}

inline void print_perf_event_error(int fd) {
	switch (fd) {
	case EINVAL:
		fprintf(stderr, "\tExtra bits in config1\n");
		fprintf(stderr, "\tCache generalized event parameter is out of range\n");
		fprintf(stderr, "\tGeneralized event setting in kernel is -1\n");
		fprintf(stderr, "\tScheduling the events failed (conflict)\n");
		fprintf(stderr, "\tToo many events\n");
		fprintf(stderr, "\tInvalid flags setting\n");
		fprintf(stderr, "\tInvalid parameters in attr\n");
		fprintf(stderr, "\tfrequency setting higher than value set by sysctl\n");
		fprintf(stderr, "\tspecified CPU does not exist\n");
		fprintf(stderr, "\tnon-group leader marked as exclusive or pinned\n"); break;
	case EACCES: fprintf(stderr, "\tRequired root permissions:\n"); break;
	case ENODEV: fprintf(stderr, "\tHardware not supported\n"); break;
	case ENOENT:
		fprintf(stderr, "\tGeneralized event set to 0 in kernel\n");
		fprintf(stderr, "\tInvalid attr.type setting\n"); break;
	case ENOSYS: fprintf(stderr, "\tPERF_SOUREC_STACK_TRACE not supported\n"); break;
	case E2BIG: fprintf(stderr, "\tattr structure bigger than expected and non-zero\n"); break;
	case EOPNOTSUPP:
		fprintf(stderr, "\tPMU interrupt not available and requested sampling\n");
		fprintf(stderr, "\tRequest branch tracing and not available\n");
		fprintf(stderr, "\tRequest low-skid events and not available\n"); break;
	case ENOMEM: fprintf(stderr, "\tKernel failed while allocating memory\n"); break;
	default: fprintf(stderr, "\tUnknown reason: %d\n", fd); break;
	}
}


#endif /* MEASURE_H_ */
