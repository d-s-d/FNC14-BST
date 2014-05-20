#include<bst.h>
#include<math.h>
#include<string.h>


/*
 * Disclaimer: The indices differ from the high-level pseudo-code description
 * found in the book [1].
 *
 * Here, e[i,j] is the expected value for the tree containing keys i to j-1,
 * since the first key is p_0 and not p_1 as in the book. w[.,.] and root[.,.]
 * are defined similarly.
 *
 */

/*
 * unroll 16 instead of 8 in bst_114
 */

#define STRIDE (n+1)
#define IDX(i,j) ((i)*STRIDE + j)

typedef struct {
    double* e;
    double* w;
    int* r;
    size_t n;
} segments_t;

void* bst_alloc_115_block_unroll_16( size_t n ) {
    segments_t* mem = (segments_t*) malloc( sizeof(segments_t) );
    size_t sz = (n+1)*(n+1);
    // XXX: for testing: calloc
    mem->e = (double*) calloc(1,  sz * sizeof(double) );
    mem->w = (double*) calloc(1,  sz * sizeof(double) );
    mem->r = (int*)    calloc(1,  sz * sizeof(int) );
    memset( mem->r, -1, sz * sizeof(int) );
    return mem;
}

#define NB 16

#define FULL_BLOCK(i, j, r, start, stop)                               \
    {                                                                  \
        double fb_t_min = e[IDX((i),(j))]; \
        int fb_r_min; \
        for ((r)=(start); (r)<(stop); ++(r)) {                         \
            t = e[IDX((i),(r))] + e[IDX((r)+1,(j))] + w[IDX((i),(j))]; \
            if (t < fb_t_min) {                                 \
                fb_t_min = t; \
                fb_r_min = r; \
            }                                                          \
        }                                                              \
        e[IDX((i),(j))] = fb_t_min;                                   \
        root[IDX((i),(j))] = fb_r_min;                              \
    }

double bst_compute_115_block_unroll_16( void*_bst_obj, double* p, double* q, size_t nn ) {
    segments_t* mem = (segments_t*) _bst_obj;
    int n;
    int i, l, r, j;
    double t, t_min, w_cur;
    int r_min;
    double* e = mem->e, *w = mem->w;
    int* root = mem->r;
    // initialization
    mem->n = nn;
    n = nn; // subtractions with n potentially negative. say hello to all the bugs
    for( i = 0; i < n+1; i++ ) {
        e[IDX(i,i)] = q[i];
        w[IDX(i,i)] = q[i];
    }

    int ib;

    // compute bottom right triangle NBxNB (i.e. bottom NB rows)
    for (i = n-1; (i >= 0) && (i > (n-NB)); --i) {
        for (j = i+1; j < n+1; ++j) {
            e[IDX(i,j)] = INFINITY;
            w[IDX(i,j)] = w[IDX(i,j-1)] + p[j-1] + q[j];
            FULL_BLOCK(i, j, r, i, j);
            /*
            for (r=i; r<j; ++r) {
                t = e[IDX(i,r)] + e[IDX(r+1,j)] + w[IDX(i,j)];
                if (t < e[IDX(i,j)]) {
                    e[IDX(i,j)] = t;
                    root[IDX(i,j)] = r;
                }
            }
            */
        }
    }

    // compute the remaining rows
    // printf("Rest of rows: i=%d\n", i);
    for (; i >= 0; --i) {
        // First, the starting NB values in this row are computed.
        // This corresponds to completing the NBxNB triangle right down from (i,i)
        for (j = i+1; j < (i+NB); ++j) {
            e[IDX(i,j)] = INFINITY;
            w[IDX(i,j)] = w[IDX(i,j-1)] + p[j-1] + q[j];
            FULL_BLOCK(i, j, r, i, j);
            /*
            for (r=i; r<j; ++r) {
                t = e[IDX(i,r)] + e[IDX(r+1,j)] + w[IDX(i,j)];
                if (t < e[IDX(i,j)]) {
                    e[IDX(i,j)] = t;
                    root[IDX(i,j)] = r;
                }
            }
            */
        }

        // Now we compute the rest of the row, but do updates in chunk of NB's.
        // Since we now have the first NB values in this row, we can compute
        // the first NB iterations of the r-loop for all the remaining values
        // in this row.
        // (this needs to be seperated from the loop afterwards since we also do
        //  initialization to INFINITY))
        for (; j < (n+1); ++j) {
            e[IDX(i,j)] = INFINITY;
            w[IDX(i,j)] = w[IDX(i,j-1)] + p[j-1] + q[j];
            FULL_BLOCK(i, j, r, i, i+NB);
            /*
            for (r=i; r < (i+NB); ++r) {
                t = e[IDX(i,r)] + e[IDX(r+1,j)] + w[IDX(i,j)];
                if (t < e[IDX(i,j)]) {
                    e[IDX(i,j)] = t;
                    root[IDX(i,j)] = r;
                }
            }
            */
        }

        // We now continue to update the values in this row in chunks of NB
        // as long as possible.
        for (ib = i+NB; (ib+NB) < (n+1); ib += NB) {
            //printf("got in here for i=%d\n", i);

            // Again we start by finishing computing the next NB values of
            // row 'i'. The last values needed for that are from the triangle
            // in down right from (ib,ib).
            for (j = (ib+1); j < (ib+NB); ++j) {
                FULL_BLOCK(i, j, r, ib, j);
                /*
                for (r=ib; r<j; ++r) {
                    t = e[IDX(i,r)] + e[IDX(r+1,j)] + w[IDX(i,j)];
                    if (t < e[IDX(i,j)]) {
                        e[IDX(i,j)] = t;
                        root[IDX(i,j)] = r;
                    }
                }
                */
            }

            // Now, having NB new values in row 'i', we compute the next NB
            // r-iterations for the remaining values in row 'i'.
            // THIS IS THE MAIN LOOP

            for (; j < (n+1); ++j) {
#if (NB != 16)
#error NB must be 16
#endif

#if 1
                double t01, t02, t03, t04, t05, t06, t07, t08;
                double t09, t10, t11, t12, t13, t14, t15, t16;
                double t0102, t0304, t0506, t0708;
                double t0910, t1112, t1314, t1516;
                double t01020304, t05060708;
                double t09101112, t13141516;
                double t0102030405060708, t0910111213141516;
                double t01020304050607080910111213141516;
                int r0102, r0304, r0506, r0708;
                int r0910, r1112, r1314, r1516;
                int r01020304, r05060708, r09101112, r13141516;
                int r0102030405060708, r0910111213141516;
                int r01020304050607080910111213141516;

                double w_cur = w[IDX(i,j)];
                int idx_i_ib = IDX(i,ib);

                t01 = e[idx_i_ib+ 0] + e[IDX(ib+ 1,j)] + w_cur;
                t02 = e[idx_i_ib+ 1] + e[IDX(ib+ 2,j)] + w_cur;
                t03 = e[idx_i_ib+ 2] + e[IDX(ib+ 3,j)] + w_cur;
                t04 = e[idx_i_ib+ 3] + e[IDX(ib+ 4,j)] + w_cur;
                t05 = e[idx_i_ib+ 4] + e[IDX(ib+ 5,j)] + w_cur;
                t06 = e[idx_i_ib+ 5] + e[IDX(ib+ 6,j)] + w_cur;
                t07 = e[idx_i_ib+ 6] + e[IDX(ib+ 7,j)] + w_cur;
                t08 = e[idx_i_ib+ 7] + e[IDX(ib+ 8,j)] + w_cur;
                t09 = e[idx_i_ib+ 8] + e[IDX(ib+ 9,j)] + w_cur;
                t10 = e[idx_i_ib+ 9] + e[IDX(ib+10,j)] + w_cur;
                t11 = e[idx_i_ib+10] + e[IDX(ib+11,j)] + w_cur;
                t12 = e[idx_i_ib+11] + e[IDX(ib+12,j)] + w_cur;
                t13 = e[idx_i_ib+12] + e[IDX(ib+13,j)] + w_cur;
                t14 = e[idx_i_ib+13] + e[IDX(ib+14,j)] + w_cur;
                t15 = e[idx_i_ib+14] + e[IDX(ib+15,j)] + w_cur;
                t16 = e[idx_i_ib+15] + e[IDX(ib+16,j)] + w_cur;

                /*
                printf("i=%d, j=%d, ib=%d\n", i, j, ib);
                printf("T1=%.3lf, "
                       "T2=%.3lf, "
                       "T3=%.3lf, "
                       "T4=%.3lf, "
                       "T5=%.3lf, "
                       "T6=%.3lf, "
                       "T7=%.3lf, "
                       "T8=%.3lf, \n",
                       t1, t2, t3, t4, t5, t6, t7);
                       */

                /*
                if (t01 < t02) { t0102 = t01; r0102 = ib+ 0; }
                else           { t0102 = t02; r0102 = ib+ 1; }
                if (t03 < t04) { t0304 = t03; r0304 = ib+ 2; }
                else           { t0304 = t04; r0304 = ib+ 3; }
                if (t05 < t06) { t0506 = t05; r0506 = ib+ 4; }
                else           { t0506 = t06; r0506 = ib+ 5; }
                if (t07 < t08) { t0708 = t07; r0708 = ib+ 6; }
                else           { t0708 = t08; r0708 = ib+ 7; }
                if (t09 < t10) { t0910 = t09; r0910 = ib+ 8; }
                else           { t0910 = t10; r0910 = ib+ 9; }
                if (t11 < t12) { t1112 = t11; r1112 = ib+10; }
                else           { t1112 = t12; r1112 = ib+11; }
                if (t13 < t14) { t1314 = t13; r1314 = ib+12; }
                else           { t1314 = t14; r1314 = ib+13; }
                if (t15 < t16) { t1516 = t15; r1516 = ib+14; }
                else           { t1516 = t16; r1516 = ib+15; }

                if (t12 < t34) { t1234 = t12; r1234 = r12; }
                else           { t1234 = t34; r1234 = r34; }
                if (t56 < t78) { t5678 = t56; r5678 = r56; }
                else           { t5678 = t78; r5678 = r78; }
                if (t1234 < t5678) { t = t1234; r = r1234; }
                else               { t = t5678; r = r5678; }
                */

#define UNROLL_MIN(a,b) \
                { if (t##a < t##b) { t##a##b = t##a; r##a##b = r##a; }   \
                  else             { t##a##b = t##b; r##a##b = r##b; } }
                const int r01 =  0;
                const int r02 =  1;
                const int r03 =  2;
                const int r04 =  3;
                const int r05 =  4;
                const int r06 =  5;
                const int r07 =  6;
                const int r08 =  7;
                const int r09 =  8;
                const int r10 =  9;
                const int r11 = 10;
                const int r12 = 11;
                const int r13 = 12;
                const int r14 = 13;
                const int r15 = 14;
                const int r16 = 15;

                UNROLL_MIN(01, 02);
                UNROLL_MIN(03, 04);
                UNROLL_MIN(05, 06);
                UNROLL_MIN(07, 08);
                UNROLL_MIN(09, 10);
                UNROLL_MIN(11, 12);
                UNROLL_MIN(13, 14);
                UNROLL_MIN(15, 16);

                UNROLL_MIN(0102, 0304);
                UNROLL_MIN(0506, 0708);
                UNROLL_MIN(0910, 1112);
                UNROLL_MIN(1314, 1516);

                UNROLL_MIN(01020304, 05060708);
                UNROLL_MIN(09101112, 13141516);

                UNROLL_MIN(0102030405060708, 0910111213141516);

                t = t01020304050607080910111213141516;
                r = r01020304050607080910111213141516;

                if (t < e[IDX(i,j)]) {
                    e[IDX(i,j)] = t;
                    root[IDX(i,j)] = r;
                }
#else
                //FULL_BLOCK(i, j, r, ib, ib+NB);
                printf("i=%d, j=%d, ib=%d\n", i, j, ib);
                for (r=ib; r<(ib+NB); ++r) {
                    t = e[IDX(i,r)] + e[IDX(r+1,j)] + w[IDX(i,j)];
                    printf("T%d=%.3lf, ", r-ib+1, t);
                    if (t < e[IDX(i,j)]) {
                        e[IDX(i,j)] = t;
                        root[IDX(i,j)] = r;
                    }
                }
                printf("\n");
#endif
            }
        }

        // There are less than NB elements remaining in row 'i'. The values
        // missing come from the triangle down right of (ib,ib)
        for (j = (ib+1); j < (n+1); ++j) {
            FULL_BLOCK(i, j, r, ib, j);
            /*
            for (r=ib; r<j; ++r) {
                t = e[IDX(i,r)] + e[IDX(r+1,j)] + w[IDX(i,j)];
                if (t < e[IDX(i,j)]) {
                    e[IDX(i,j)] = t;
                    root[IDX(i,j)] = r;
                }
            }
            */
        }
    }

    return e[IDX(0,n)];
}

size_t bst_get_root_115_block_unroll_16( void* _bst_obj, size_t i, size_t j )
{
    // [i,j], in table: [i-1, j]+1
    segments_t *mem = _bst_obj;
    size_t n = mem->n;
    int *root = mem->r;
    return (size_t) root[(i-1)*(n+1)+j]+1;
}

void bst_free_115_block_unroll_16( void* _mem ) {
    segments_t* mem = (segments_t*) _mem;

    /*
    size_t n = mem->n;
    for (size_t i=0; i<=n; ++i) {
        for (size_t j=0; j<=n; ++j) {
            printf(" %.3lf", mem->e[i*(n+1)+j]);
        }
        printf("\n");
    }
    */

    free( mem->e );
    free( mem->w );
    free( mem->r );
    free( mem );
}

size_t bst_flops_115_block_unroll_16( size_t n ) {
    size_t n3 = n*n*n;
    size_t n2 = n*n;
    return (n3 + 5*n2)/2 + 2*n;
}
