// copyright: see the header of "qsp.c"
#include "stdlib2.h"
#include "iqsort.h"

double *A;

// index qsort
void iqsort(UINT64 *index, double *score, UINT64 n) {
  UINT64 i;

  A = (double *)malloc(sizeof(double) * n);
  memcpy(A, score, sizeof(double) * n);

  for (i = 0; i < n; i++) {
    index[i] = i;
  }

  qsort(index, n, sizeof(UINT64), compare);

  free(A);
}

// for index qsort
int compare(const void *a, const void *b) {
  if (A[*(int *)a] > A[*(int *)b]) {
    return -1;
  } else if (A[*(int *)a] == A[*(int *)b]) {
    return 0;
  }
  return 1;
}
