// Provides other_computation to test whether the cost model works fine with
// multiple files involved.

#include "cost_model.h"

#define N 2000000
static double A[N];
static double B[N];
static double C[N];

void other_computation()
{
    for (int i=0; i<N; ++i) {
        COST_INC_ADD(1);
        COST_INC_MUL(2);
        COST_INC_DIV(1);

        C[N] = (3.0*A[N]) / B[N] + 2.0*C[N];
    }
}
