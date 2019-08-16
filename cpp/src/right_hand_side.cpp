#include "right_hand_side.h"

#include <math.h>

real RightHandSide::l2residual(const aligned_vector<real> &c) {
  aligned_vector<real> dcdt;
  (*this)(c, dcdt, 0);
  real res = 0.;
#pragma omp parallel for simd reduction(+:res)
  for (int i = 0; i < dcdt.size(); ++i) {
    res += dcdt[i] * dcdt[i];
  }
  return std::sqrt(res);
}
