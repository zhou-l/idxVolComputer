// copyright: see the header of "qsp.c"
#include "../StdDefines.h"
// index qsort
extern void iqsort(UINT64 *index, double *score, UINT64 n);

// comparison function for index qsort
extern int compare(const void *a , const void *b);
