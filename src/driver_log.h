#ifndef DRIVER_LOG_H
#define  DRIVER_LOG_H

/**
 * Setup data logging.
 * @param filename: name of logfile
 * @returns:        true if logging can be used, false otherwise
 */
int log_setup(char *filename);

/**
 * Close logging file etc.
 */
void log_cleanup();

/**
 * Use this to log an int.
 * @param name: can be NULL (no tag, e.g. in arrays)
 */
void log_int(char *name, int v);

/**
 * Use this to log an integral double.
 * @param name: can be NULL (no tag, e.g. in arrays)
 */
void log_idouble(char *name, double v);

/**
 * Use this to log a size_t.
 * @param name: can be NULL (no tag, e.g. in arrays)
 */
void log_size(char *name, size_t v);

/**
 * Use this to log a string.
 * @param name: can be NULL (no tag, e.g. in arrays)
 */
void log_str(char *name, char *str);

/**
 * Use this if you need a custom format.
 *
 * E.g, you want to log a double as percentage, i.e. precision 2.
 * You can use
 *   log_fmt("some-percentage", "%.2lf", v)
 *
 * @param name: can be NULL (no tag, e.g. in arrays)
 */
void log_fmt(char *name, char *fmt, ...);

/**
 * Use this to declare a new struct block.
 * @param name: can be NULL (no tag, e.g. in arrays)
 */
void log_struct(char *name);

/**
 * Use this to close the struct block.
 */
void log_struct_end();

/**
 * Use this to declare a new array block.
 * @param name: can be NULL (no tag, e.g. in arrays)
 */
void log_array(char *name);

/**
 * Use this to close the array block.
 */
void log_array_end();

#endif
