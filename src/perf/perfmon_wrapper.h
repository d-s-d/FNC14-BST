/**
 * Helper functions to easily use performance counters on modern Intel
 * processors.
 *
 * 2014, Jeremia Baer <baerj@student.ethz.ch>
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>

#include <perfmon/pfmlib.h>
#include <perfmon/pfmlib_perf_event.h>

/**
 * Performance Measurment Configuration Data
 */
struct perf_data {
    char     *name;
    int      fd;
    uint64_t type;
    uint64_t config;
    uint64_t value;
};

/**
 * Initializes performance measurements using libpfm and perf_events.
 *
 * NOTE: Even if perf_init fails, you have to call perf_cleanup on the data.
 *
 * @param events: list of event names to measure
 * @param data:   address of pointer to store data. allocated by perf_init
 * @returns:      0 on success, -1 on failure
 */
static int perf_init(char **events, struct perf_data **data);

/**
 * Reset counter values.
 */
static inline void perf_reset(struct perf_data *data);

/**
 * Starts the performance counters of the given configuration.
 */
static inline void perf_start(struct perf_data *data);

/**
 * Stops performance measurements of the given configuration.
 */
static inline void perf_stop(struct perf_data *data);

/**
 * Updates the counter values of the given configuration.
 */
static int perf_update_values(struct perf_data *data);

/**
 * Frees data in struct perf_data.
 */
static void perf_cleanup(struct perf_data *data);

// ------------------------------------------------------------------------
// IMPLEMENTATION from here on
// ------------------------------------------------------------------------

#define PERF_ERR(...)                  \
    fprintf(stderr, "%s: ", __func__); \
    fprintf(stderr, __VA_ARGS__)

static inline void perf_fd_error(int fd) {
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

// if failed, still need to call cleanup!
static int perf_init(char **events, struct perf_data **data)
{
    int ret;      // generic return value
    *data = NULL; // in case of failure

    // check for valid argument
    if (!events) {
        PERF_ERR("No events given!\n");
        return -1;
    }

    // get number of events
    int event_count = 0;
    for (; events[event_count]; event_count++) { }
    if (event_count == 0) {
        PERF_ERR("No events given!\n");
        return -1;
    }

    // initialize libpfm
    ret = pfm_initialize();
    if (ret != PFM_SUCCESS) {
        PERF_ERR("Error initilizing pfm library: %s\n", pfm_strerror(ret));
        return -1;
    }

    // get event's perf configuration. perf_data NULL-name-terminated
    struct perf_data *perf_data = malloc(
            (event_count+1)*sizeof(struct perf_data));
    *data = perf_data;
    if (!perf_data) {
        PERF_ERR("Could not allocate memory.\n");
        return -1;
    }
    memset(perf_data, 0, (event_count+1)*sizeof(struct perf_data));

    for (int i=0; i<event_count; ++i) {
        char *fqstr = NULL;
        struct perf_event_attr perf_attr;
        memset(&perf_attr, 0, sizeof(perf_attr));

        pfm_perf_encode_arg_t e;
        memset(&e, 0, sizeof(e));
        e.fstr = &fqstr;
        e.attr = &perf_attr;

        ret = pfm_get_os_event_encoding(events[i],
                PFM_PLM3, PFM_OS_PERF_EVENT, &e);
        if (ret != PFM_SUCCESS) {
            PERF_ERR("Cannot encode event %s: %s\n",
                    events[i], pfm_strerror(ret));
            return -1;
        }

        perf_data[i].name   = fqstr;
        perf_data[i].fd     = -1;
        perf_data[i].type   = perf_attr.type;
        perf_data[i].config = perf_attr.config;
        perf_data[i].value  = 0;
    }

    // printf("events converted. they are:\n");
    // for (int i=0; i<event_count; ++i) {
    //  printf("%d: %s\n", i, perf_data[i].name);
    // }
    
    // now that we have the necessary perf configurations, let's setup
    // the measurements.
    for (int i=0; i<event_count; ++i) {
        struct perf_event_attr pe;
        memset(&pe, 0, sizeof(struct perf_event_attr));

        pe.type           = perf_data[i].type;
        pe.size           = sizeof(struct perf_event_attr);
        pe.config         = perf_data[i].config;
        pe.disabled       = 1;
        pe.exclude_kernel = 1;
        pe.exclude_hv     = 1;
        pe.read_format    = PERF_FORMAT_GROUP | PERF_FORMAT_ID;

        perf_data[i].fd = syscall(__NR_perf_event_open, &pe, 0, -1,
                perf_data->fd, 0);
        if (perf_data[i].fd == -1) {
            PERF_ERR("Error opening event %"PRIx64" (%s)\n",
                    pe.config, perf_data[i].name);
            perf_fd_error(perf_data[i].fd);
            return -1;
        }

    }

    perf_reset(perf_data);

    return 0;
}

static inline void perf_reset(struct perf_data *data)
{
    ioctl(data->fd, PERF_EVENT_IOC_RESET,  PERF_IOC_FLAG_GROUP);
}

static inline void perf_start(struct perf_data *data)
{
    ioctl(data->fd, PERF_EVENT_IOC_ENABLE, PERF_IOC_FLAG_GROUP);
}

static inline void perf_stop(struct perf_data *data)
{
    ioctl(data->fd, PERF_EVENT_IOC_DISABLE,PERF_IOC_FLAG_GROUP);
}

static int perf_update_values(struct perf_data *data)
{
    int event_count = 0;
    for(; data[event_count].name; ++event_count);

    int values_size = (1 + 2*event_count)*sizeof(uint64_t);
    uint64_t *values = malloc(values_size);
    if (!values) {
        PERF_ERR("Could not allocate memory\n");
        return -1;
    }

    // XXX: ret == sizeof(values) ??
    read(data->fd, values, values_size);
    for (int i=0; i<event_count; ++i) {
        data[i].value = values[1+2*i];
    }

    free(values);
    return 0;
}

static void perf_cleanup(struct perf_data *data)
{
    if (!data) return;

    // free all names
    for (struct perf_data *iter = data; iter->name; ++iter) {
        close(iter->fd);
        free(iter->name);
    }

    // free struct
    free(data);

    pfm_terminate();
}
