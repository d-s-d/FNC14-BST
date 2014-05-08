#include "mmm.h"

void compute() {
  int i,j,h;
  for(i = 0; i < m; ++i) {
    for(j = 0; j < n; ++j) {
      for(h = 0; h < k; ++h) {
        C[i*n+j] += A[i*k+h] * B[h*n+j];
      }
    }
  }
}
