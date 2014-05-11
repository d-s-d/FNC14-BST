#ifndef URANGE_H
#define URANGE_H

#include<stddef.h>

typedef struct {
    size_t start;
    size_t stop;
    size_t step;
} urange_t;

size_t urange_get_max_val( urange_t* urange );
void urange_print( urange_t* urange );
int read_urange( urange_t* urange, char* str );
size_t urange_get_size( urange_t* urange );

#endif // URANGE_H
