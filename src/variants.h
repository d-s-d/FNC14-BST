#ifndef VARIANTS_H
#define VARIANTS_H

void* bst_alloc_001_transposed( size_t n );
double bst_compute_001_transposed( void*_bst_obj, double* p, double* q,
                                   size_t n );

size_t bst_get_root_001_transposed( void* _bst_obj, size_t i, size_t j );
void bst_free_001_transposed( void* _mem );

#endif // VARIANTS_H
