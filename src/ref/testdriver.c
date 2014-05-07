#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>

#include"bst_test.h"

typedef struct {
  bool (*f)();
  char* desc;
} testfunc_t;

// include list of files

// #include "matrix_test.h"

/*
 * Returning an empty string means success. Any non-empty string means failure.
 * In case of a failure, the test function is supposed to return a string
 * describing the reason of failure.
 *
 */
bool dummy_success( ) {
    return true;
}

bool dummy_fail() {
    printf("This function failed for some reason.\n");
    return false;
}

size_t register_test( testfunc_t** tests, size_t n, bool (*f)(),
                      char* desc ) {
    size_t new_n = n+1;
    *tests = (testfunc_t*) realloc( *tests, (new_n) * sizeof( testfunc_t ) );
    (*tests)[new_n-1] = (testfunc_t) { f, desc };
    return new_n;
}

// register the function that you want to test here!

size_t register_tests( testfunc_t** tests, size_t n ) {
    n = register_test( tests, n, dummy_success, "dummy_success" );
    n = register_test( tests, n, dummy_fail, "dummy_fail" );
    n = register_test( tests, n, simple_bst_test,
                       "simple_bst_test" );
}

void run_tests( testfunc_t* tests, size_t n ) {
    size_t i;
    bool res;
    for( i = 0; i < n; i++ ) {
        printf("\n*** Running test: %s ***\n", tests[i].desc);
        res = tests[i].f();
        if( res ) printf("%s: SUCCEEDED\n", tests[i].desc);
        else      printf("%s: FAILED\n", tests[i].desc);
    }
}

int main( int argc, char* argv[] ) {
    testfunc_t* tests = NULL;
    size_t test_count = 0;

    test_count = register_tests( &tests, test_count );
    run_tests( tests, test_count );

    free(tests);
    return 0;
}
