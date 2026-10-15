#ifndef DREAMPLACE_UTILITY_OMP_H
#define DREAMPLACE_UTILITY_OMP_H
#ifdef _OPENMP
#include <omp.h>
#else
// The same CPU kernels run serially when OpenMP is not requested.
inline int omp_get_thread_num() { return 0; }
inline int omp_get_max_threads() { return 1; }
#endif
#endif
