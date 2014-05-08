#include <stdio.h>
#include <stdarg.h>

#include "driver_log.h"

static FILE *fd = NULL;

static int log_indent  = 0;
#define log_do_indent() \
    {for (int i=0; i<log_indent; ++i) fprintf(fd, "  ");}

int log_setup(char *filename)
{
    if (filename) {
        fd = fopen(filename, "w");
    }

    if (fd) {
        fprintf(fd, "{\n");
        log_indent++;
    } else {

    }

    return (fd != NULL);
}

void log_cleanup()
{
    if (fd) {
        fprintf(fd, "}\n");
        fclose(fd);
    }
}

void log_fmt(char *name, char *fmt, ...)
{
    if (fd) {
        log_do_indent();
        if (name) {
            fprintf(fd, "'%s' : ", name);
        }

        va_list argptr;
        va_start(argptr, fmt);
        vfprintf(fd, fmt, argptr);
        va_end(argptr);

        fprintf(fd, "\n");
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

static void log_part(char *name, char delim)
{
    if (fd) {
        log_do_indent();
        if (name) {
            fprintf(fd, "'%s' : %c\n", name, delim);
        } else {
            fprintf(fd, "%c\n", delim);
        }
        log_indent++;
    }
}

static void log_part_end(char delim)
{
    if (fd) {
        log_indent--;
        log_do_indent();
        fprintf(fd, "%c\n", delim);
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
