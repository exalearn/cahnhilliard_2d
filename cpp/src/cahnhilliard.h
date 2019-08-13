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
  void rhs(const aligned_vector<double> &c, aligned_vector<double> &dcdt, const double t) override;
  void write_state(const aligned_vector<double> &x , const int idx , const int nx , const int ny , std::string& outdir) override;
  void setInitialConditions(aligned_vector<double> &x);
  
 private:

  CHparamsVector chpV_;
  SimInfo& info_;
  void (*ch_rhs_) (const aligned_vector<double>&, aligned_vector<double>&, double, CHparamsVector&, SimInfo&);

  std::default_random_engine generator_;
  std::normal_distribution<double> noise_dist_;
    
};

#endif
