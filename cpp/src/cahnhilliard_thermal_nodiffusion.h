#ifndef __CAHNHILLIARD_THERMAL_NODIFFUSION_H__
#define __CAHNHILLIARD_THERMAL_NODIFFUSION_H__

#include <vector>
#include <algorithm>
#include <random>
#include "timer.h"
#include "chparams.h"
#include "right_hand_side.h"

class CahnHilliard2DRHS_thermal_nodiffusion : public RightHandSide {

 public:

  CahnHilliard2DRHS_thermal_nodiffusion(CHparamsScalar& chp , SimInfo& info);
  CahnHilliard2DRHS_thermal_nodiffusion(CHparamsVector& chp , SimInfo& info);
  ~CahnHilliard2DRHS_thermal_nodiffusion();
  void rhs(const aligned_vector<real> &c, aligned_vector<real> &dcdt, const real t) override;
  void write_state( const aligned_vector<real> &x , const int idx , const int nx , const int ny , std::string& outdir) override;
  void setInitialConditions(aligned_vector<real> &x);
  void printTimers() const;
  
 private:

  CHparamsVector chpV_;
  SimInfo& info_;
  void (*ch_rhs_) (const aligned_vector<real>&, aligned_vector<real>&, real, CHparamsVector&, SimInfo&);

  std::default_random_engine generator_;
  std::normal_distribution<real> noise_dist_;
  
  timer t_params, t_nonlocal;
};





#endif
