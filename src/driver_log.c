#include <stdio.h>
#include <stdarg.h>

#include "driver_log.h"

static FILE *fd = NULL;

static int log_indent  = 0;
static int log_level_first[100];

static void log_go_in() {
    log_indent++;
    log_level_first[log_indent] = 1;
}

static void log_go_out() {
    fprintf(fd, "\n");
    log_indent--;
}

static void log_do_indent() {
    for (int i=0; i<log_indent; ++i) fprintf(fd, "  ");
}

static void log_do_comma() {
    if (!log_level_first[log_indent]) {
        fprintf(fd, ",\n");
    } else {
        fprintf(fd, "\n");
        log_level_first[log_indent] = 0;
    }
}

static void log_new_entry() {
    log_do_comma();
    log_do_indent();
}

void log_fmt(char *name, char *fmt, ...)
{
    if (fd) {
        log_new_entry();
        if (name) {
            fprintf(fd, "\"%s\" : ", name);
        }

        va_list argptr;
        va_start(argptr, fmt);
        vfprintf(fd, fmt, argptr);
        va_end(argptr);
    }
}

void log_int(char *name, int v)
{
    log_fmt(name, "%d", v);
}

void log_size_t(char* name, size_t v)
{
    log_fmt(name, "%zu", v);
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
    log_fmt(name, "\"%s\"", str);
}

static void log_part(char *name, char delim)
{
    if (fd) {
        log_new_entry();
        if (name) {
            fprintf(fd, "\"%s\" : %c", name, delim);
        } else {
            fprintf(fd, "%c", delim);
        }
        log_go_in();
    }
}

static void log_part_end(char delim)
{
    if (fd) {
        log_go_out();
        log_do_indent();
        fprintf(fd, "%c", delim);
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

int log_setup(char *filename)
{
    if (filename) {
        fd = fopen(filename, "w");
    }

    if (fd) {
        fprintf(fd, "{");
        log_go_in();
    } else {

    }

    return (fd != NULL);
}

void log_cleanup()
{
    if (fd) {
        log_go_out();
        fprintf(fd, "}\n");
        fclose(fd);
    }

    if (log_indent != 0) {
        fprintf(stderr, "LOGERR: Indent level not 0 at end. Called array/struct_end?\n");
    }
}
