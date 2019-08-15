#ifndef __ALIGNED_VECTOR_H__
#define __ALIGNED_VECTOR_H__
#include <vector>
#include "allocator.h"

//some useful definitions
template <typename T>
using aligned_vector = std::vector< T, aligned_allocator<T, 64> >;

using aligned_double_vector = aligned_vector<double>;

#endif
