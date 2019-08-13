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
  void rhs(const aligned_vector<double> &c, aligned_vector<double> &dcdt, const double t) override;
  void write_state( const aligned_vector<double> &x , const int idx , const int nx , const int ny , std::string& outdir) override;
  void setInitialConditions(aligned_vector<double> &x);
  void printTimers() const;
  
 private:

  CHparamsVector chpV_;
  SimInfo& info_;
  void (*ch_rhs_) (const aligned_vector<double>&, aligned_vector<double>&, double, CHparamsVector&, SimInfo&);

  std::default_random_engine generator_;
  std::normal_distribution<double> noise_dist_;
  
  timer t_params, t_nonlocal;
};





#endif
