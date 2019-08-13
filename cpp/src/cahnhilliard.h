#ifndef __CAHNHILLIARD_H__
#define __CAHNHILLIARD_H__

#include <vector>
#include <algorithm>
#include <random>
#include "chparams.h"
#include "right_hand_side.h"
#include "timer.h"

class CahnHilliard2DRHS : public RightHandSide {

 public:

  CahnHilliard2DRHS(CHparamsScalar& chp , SimInfo& info);
  CahnHilliard2DRHS(CHparamsVector& chp , SimInfo& info);
  ~CahnHilliard2DRHS();
  void rhs(const aligned_vector<real> &c, aligned_vector<real> &dcdt, const real t) override;
  void write_state(const aligned_vector<real> &x , const int idx , const int nx , const int ny , std::string& outdir) override;
  void setInitialConditions(aligned_vector<real> &x);
  
 private:

  CHparamsVector chpV_;
  SimInfo& info_;
  void (*ch_rhs_) (const aligned_vector<real>&, aligned_vector<real>&, real, CHparamsVector&, SimInfo&);

  std::default_random_engine generator_;
  std::normal_distribution<real> noise_dist_;
    
};

#endif
