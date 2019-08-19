#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <omp.h>
#include <boost/numeric/odeint.hpp>
#include "cahnhilliard_thermal_nodiffusion.h"
#include "utils_ch.h"
  /*
  Cahn-Hilliard:
  
  dc/dt = laplacian( u*c^3 - b*c ) - eps_2*biharm(c) - sigma*(c - m) + sigma_noise * N(0,1^2)
  
  expanding out RHS into individual differentials:
  D*laplacian( u*c^3 - b*c) - D*eps_2*biharm(c)
  assuming constant eps_2.

  need a d^4 and a d^2 operator.
  */

CahnHilliard2DRHS_thermal_nodiffusion::CahnHilliard2DRHS_thermal_nodiffusion(CHparamsScalar& chp , SimInfo& info)
  : noise_dist_(0.0,1.0), chpV_(*(new CHparamsVector())), info_(info)
  {    
    chpV_.eps_2    = aligned_vector<real>( info_.nx*info_.ny , chp.eps_2     );
    chpV_.b        = aligned_vector<real>( info_.nx*info_.ny , chp.b         );
    chpV_.u        = aligned_vector<real>( info_.nx*info_.ny , chp.u         );
    chpV_.sigma    = aligned_vector<real>( info_.nx*info_.ny , chp.sigma     );
    chpV_.m        = aligned_vector<real>( info_.nx*info_.ny , chp.m  );
    chpV_.DT       = aligned_vector<real>( info_.nx*info_.ny , chp.DT  );
    chpV_.f_T      = aligned_vector<real>( info_.nx*info_.ny , chp.f_T  );
    chpV_.sigma_noise    = chp.sigma_noise;
    chpV_.T_const        = aligned_vector<real>( info_.nx*info_.ny , chp.T_const  );

    if ( info.bc.compare("dirichlet") == 0) {
      ch_rhs_ = &compute_ch_nonlocal_stationary_boundaries;
      if(info.verbosity>=1)
        std::cout << "Initialized Cahn-Hilliard equation: scalar parameters, dirichlet BCs, thermal coefficient dependence, no thermal diffusion" << std::endl;
    }
    else if ( info.bc.compare("neumann") == 0) {
      ch_rhs_ = &compute_ch_nonlocal_neumannBC;
      if(info.verbosity>=1)
        std::cout << "Initialized Cahn-Hilliard equation: scalar parameters, neumann BCs, thermal coefficient dependence, no thermal diffusion" << std::endl;
    }
    else {
      ch_rhs_ = &compute_ch_nonlocal;
      if(info.verbosity>=1)
        std::cout << "Initialized Cahn-Hilliard equation: scalar parameters, periodic BCs, thermal coefficient dependence, no thermal diffusion" << std::endl;
    }
  }

CahnHilliard2DRHS_thermal_nodiffusion::CahnHilliard2DRHS_thermal_nodiffusion(CHparamsVector& chp , SimInfo& info)
  : noise_dist_(0.0,1.0) , chpV_(chp) , info_(info)
  {

    if ( info.bc.compare("dirichlet") == 0) {
      ch_rhs_ = &compute_ch_nonlocal_stationary_boundaries;
      if(info.verbosity>=1)
        std::cout << "Initialized Cahn-Hilliard equation: spatial-field parameters, dirichlet BCs, thermal coefficient dependence, no thermal diffusion" << std::endl;
    }
    else if ( info.bc.compare("neumann") == 0) {
      ch_rhs_ = &compute_ch_nonlocal_neumannBC;
      if(info.verbosity>=1)
        std::cout << "Initialized Cahn-Hilliard equation: spatial-field parameters, neumann BCs, thermal coefficient dependence, no thermal diffusion" << std::endl;
    }
    else {
      ch_rhs_ = &compute_ch_nonlocal;
      if(info.verbosity>=1)
        std::cout << "Initialized Cahn-Hilliard equation: spatial-field parameters, periodic BCs, thermal coefficient dependence, no thermal diffusion" << std::endl;
    }
    
  }

CahnHilliard2DRHS_thermal_nodiffusion::~CahnHilliard2DRHS_thermal_nodiffusion() { };

void CahnHilliard2DRHS_thermal_nodiffusion::rhs(const aligned_vector<real> &c, aligned_vector<real> &dcdt, const real t)
  {
    
    dcdt.resize(info_.nx*info_.ny);
    
    // evaluate CH parameter dependencies on temperature
    //chpV_ = compute_chparams_using_temperature( chpV_ , info_ , chpV_.T_const );
    compute_eps2_and_sigma_from_polymer_params( chpV_ , info_ , chpV_.T_const );
    
    // evaluate deterministic nonlocal dynamics
    compute_ch_nonlocal(c, dcdt, t, chpV_, info_);
  }


void CahnHilliard2DRHS_thermal_nodiffusion::setInitialConditions(aligned_vector<real> &x) const
  { 
    x.resize(info_.nx * info_.ny);

    std::default_random_engine generator;
    std::uniform_real_distribution<real> distribution(-1.0,1.0);

    for (int i = 0; i < info_.ny; ++i) {
      for (int j = 0; j < info_.nx; ++j) {
        x[info_.idx2du(i,j)]   = distribution(generator) * 0.005;
      }
    }

    // Set BCs if needed
    if ( info_.bc.compare("dirichlet") == 0) {
      x = apply_dirichlet_bc( x , info_ );
    }
    else if ( info_.bc.compare("neumann") == 0 ) {
      x = apply_neumann_bc( x , info_ );
    }
  }


void CahnHilliard2DRHS_thermal_nodiffusion::write_state(const aligned_vector<real> &x , const int idx , const int nx , const int ny , std::string& outdir) const
{
  if ( outdir.back() != '/' )
    outdir += '/';
  std::ofstream outC;
  outC.open( outdir + "C_" + std::to_string(idx) + ".out" );
  outC.precision(16);
  
  for (int i = 0; i < ny; ++i){
    for (int j = 0; j < nx; ++j){
      outC << x[i * ny + j] << " ";
    }
  }

  outC.close();
}
