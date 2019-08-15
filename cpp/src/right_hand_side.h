#ifndef __RIGHT_HAND_SIDE_H__
#define __RIGHT_HAND_SIDE_H__

#include <vector>
#include <string>
#include "chparams.h"

class RightHandSide {
public:
  virtual void rhs(const aligned_vector<real> &c, aligned_vector<real> &dcdt,
                   const real t) = 0;
  virtual void write_state( const aligned_vector<real> &x , const int idx , const int nx , const int ny , std::string& outdir) = 0;
  virtual void setInitialConditions(aligned_vector<real> &x) = 0;
  void operator()(const aligned_vector<real> &c, aligned_vector<real> &dcdt, const real t)
  {
    rhs(c,dcdt,t);
  }
  real l2residual(const aligned_vector<real> &c);
  
};

#endif
