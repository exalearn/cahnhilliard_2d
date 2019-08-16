#ifndef __CAHNHILLIARD_THERMAL_H__
#define __CAHNHILLIARD_THERMAL_H__

#include <vector>
#include <algorithm>
#include <random>
#include "chparams.h"
#include "right_hand_side.h"

class CahnHilliard2DRHS_thermal : public RightHandSide {

 public:

  CahnHilliard2DRHS_thermal(CHparamsScalar& chp , SimInfo& info);
  CahnHilliard2DRHS_thermal(CHparamsVector& chp , SimInfo& info);
  ~CahnHilliard2DRHS_thermal();
  void rhs(const aligned_vector<real> &c, aligned_vector<real> &dcdt, const real t) override;
  void write_state( const aligned_vector<real> &x , const int idx , const int nx , const int ny , std::string& outdir ) const override;
  void setInitialConditions(aligned_vector<real> &x) const;

  //  struct PetscContext {
  //  PetscContext(CahnHilliard2DRHS_thermal &outer) : instance_(outer) {}
    
  //  CahnHilliard2DRHS_thermal &instance_;
  //};

private:

  CHparamsVector& chpV_;
  SimInfo& info_;
  void (*ch_rhs_) (const aligned_vector<real>&, aligned_vector<real>&, real, CHparamsVector&, SimInfo&);

  std::default_random_engine generator_;
  std::normal_distribution<real> noise_dist_;

};



#endif
