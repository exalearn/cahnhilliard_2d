#ifndef __RIGHT_HAND_SIDE_H__
#define __RIGHT_HAND_SIDE_H__

#include <vector>
#include <string>
#include "allocator.h"

class RightHandSide {
public:
  virtual void rhs(const aligned_vector<double> &c, aligned_vector<double> &dcdt,
                   const double t) = 0;
  virtual void write_state( const aligned_vector<double> &x , const int idx , const int nx , const int ny , std::string& outdir) = 0;
  virtual void setInitialConditions(aligned_vector<double> &x) = 0;
  void operator()(const aligned_vector<double> &c, aligned_vector<double> &dcdt, const double t)
  {
    rhs(c,dcdt,t);
  }
  double l2residual(const aligned_vector<double> &c);
  
};

#endif
