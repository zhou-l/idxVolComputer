// copyright: see the header of "qsp.c"
// compute qsp for every data point
#include "../StdDefines.h"
extern void qsp(double *X, UINT64 n, UINT64 d, UINT64 n_sample, UINT64 seed, double *score);

// random sampling
extern void sampling(UINT64 n, UINT64 min, UINT64 max, UINT64 seed, UINT64 *id_sample);
extern UINT64 checkArray(UINT64 id_current, UINT64 n_sample, UINT64 *id_sample);

// normalization of data (divide by SDs for each dimension)
extern void normalize(double *X, UINT64 n, UINT64 d);
