#include"urange.h"
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>

size_t urange_get_max_val( urange_t* urange ) {
    assert( urange != NULL );
    size_t diff = urange->stop - urange->start;
    size_t q = diff/urange->step;
    return urange->start + q*(urange->step);
}

void urange_print( urange_t* urange ) {
    printf("%.4lu:%.4lu:%.4lu", urange->start, urange->step, urange->stop);
}

int read_urange( urange_t* urange, char* str ) {
    char** endptr = &str;
    urange->step = 1;
    urange->start = strtoul( str, endptr, 10 );
    urange->stop = urange->start;
    if( **endptr == ':' ) {
        urange->stop = strtod( (*endptr)+1, endptr );
        if( **endptr == ':' ) {
            urange->step = urange->stop;
            urange->stop = strtod( (*endptr)+1, endptr );
        }
    }
    if ( **endptr != '\0' )
        return -1;
    return 0;
}

size_t urange_get_size( urange_t* urange ) {
    assert( urange != NULL );
    if( urange->start > urange->stop ) return 0;
    return ((urange->stop - urange->start)/urange->step)+1;
}
