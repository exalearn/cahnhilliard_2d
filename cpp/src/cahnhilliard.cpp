#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <omp.h>
#include <boost/numeric/odeint.hpp>
#include <string>
#include "cahnhilliard.h"
#include "utils_ch.h"


  /*
  Cahn-Hilliard:
  
  dc/dt = laplacian( u*c^3 - b*c ) - eps_2*biharm(c) - sigma*(c - m) + sigma_noise * N(0,1^2)
  
  expanding out RHS into individual differentials:
  D*laplacian( u*c^3 - b*c) - D*eps_2*biharm(c)
  assuming constant eps_2.

  need a d^4 and a d^2 operator.
  */

CahnHilliard2DRHS::CahnHilliard2DRHS(CHparamsScalar& chp , SimInfo& info)
  : noise_dist_(0.0,1.0), chpV_(*(new CHparamsVector())), info_(info)
  {    
    chpV_.eps_2    = aligned_vector<real>( info_.nx*info_.ny , chp.eps_2     );
    chpV_.b        = aligned_vector<real>( info_.nx*info_.ny , chp.b         );
    chpV_.u        = aligned_vector<real>( info_.nx*info_.ny , chp.u         );
    chpV_.sigma    = aligned_vector<real>( info_.nx*info_.ny , chp.sigma     );
    chpV_.m        = aligned_vector<real>( info_.nx*info_.ny , chp.m  );
    chpV_.sigma_noise    = chp.sigma_noise;

    if ( info.bc.compare("dirichlet") == 0 ) {
      ch_rhs_ = &compute_ch_nonlocal_stationary_boundaries;
      if(info.verbosity>=1)
        std::cout << "Initialized Cahn-Hilliard equation: scalar parameters, dirichlet BCs, no thermal dependence" << std::endl;
    }
    else if ( info.bc.compare("neumann") == 0 ) {
      ch_rhs_ = &compute_ch_nonlocal_neumannBC;
      if(info.verbosity>=1)
        std::cout << "Initialized Cahn-Hilliard equation: scalar parameters, neumann BCs, no thermal dependence" << std::endl;
    }
    else if ( info.bc.compare("mixed_neumann_bottom_dirichlet") == 0 ) {
      ch_rhs_ = &compute_ch_nonlocal_mixedBC_neumann_with_bottom_dirichlet;
      if(info.verbosity>=1)
        std::cout << "Initialized Cahn-Hilliard equation: scalar parameters, mixed BC (neumann + bottom dirichlet), no thermal dependence" << std::endl;
    }
    else if ( info.bc.compare("mixed_neumann_top_dirichlet") == 0 ) {
      ch_rhs_ = &compute_ch_nonlocal_mixedBC_neumann_with_top_dirichlet;
      if(info.verbosity>=1)
        std::cout << "Initialized Cahn-Hilliard equation: scalar parameters, mixed BC (neumann + top dirichlet), no thermal dependence" << std::endl;
    }
    else {
      ch_rhs_ = &compute_ch_nonlocal;
      if(info.verbosity>=1)
        std::cout << "Initialized Cahn-Hilliard equation: scalar parameters, periodic BCs, no thermal dependence" << std::endl;
    }
    
  }


CahnHilliard2DRHS::CahnHilliard2DRHS(CHparamsVector& chp , SimInfo& info)
  : noise_dist_(0.0,1.0) , chpV_(chp) , info_(info)
  {
    if ( info.bc.compare("dirichlet") == 0) {
      ch_rhs_ = &compute_ch_nonlocal_stationary_boundaries;
      if(info.verbosity>=1)
        std::cout << "Initialized Cahn-Hilliard equation with spatial-field parameters, dirichlet BCs, no thermal dependence" << std::endl;
    }
    else if ( info.bc.compare("neumann") == 0) {
      ch_rhs_ = &compute_ch_nonlocal_neumannBC;
      if(info.verbosity>=1)
        std::cout << "Initialized Cahn-Hilliard equation with spatial-field parameters, neumann BCs, no thermal dependence" << std::endl;
    }
    else if ( info.bc.compare("mixed_neumann_bottom_dirichlet") == 0 ) {
      ch_rhs_ = &compute_ch_nonlocal_mixedBC_neumann_with_bottom_dirichlet;
      if(info.verbosity>=1)
        std::cout << "Initialized Cahn-Hilliard equation: scalar parameters, mixed BC (neumann + bottom dirichlet), no thermal dependence" << std::endl;
    }
    else if ( info.bc.compare("mixed_neumann_top_dirichlet") == 0 ) {
      ch_rhs_ = &compute_ch_nonlocal_mixedBC_neumann_with_top_dirichlet;
      if(info.verbosity>=1)
        std::cout << "Initialized Cahn-Hilliard equation: scalar parameters, mixed BC (neumann + top dirichlet), no thermal dependence" << std::endl;
    }
    else {
      ch_rhs_ = &compute_ch_nonlocal;
      if(info.verbosity>=1)
        std::cout << "Initialized Cahn-Hilliard equation: scalar parameters, periodic BCs, no thermal dependence" << std::endl;
    }
  }


CahnHilliard2DRHS::~CahnHilliard2DRHS() { };


void CahnHilliard2DRHS::rhs(const aligned_vector<real> &c, aligned_vector<real> &dcdt, const real t)
  {
    dcdt.resize(info_.nx * info_.ny);
    (*ch_rhs_)(c, dcdt, t, chpV_, info_);
  }


void CahnHilliard2DRHS::setInitialConditions(aligned_vector<real> &x) const
  { 
    //resize if necessary
    x.resize(info_.nx * info_.ny);
    
    std::default_random_engine generator;
    std::uniform_real_distribution<real> distribution(-1.0,1.0);

    // real initial_value = -1.0;
    for (int i = 0; i < info_.ny; ++i) {
      for (int j = 0; j < info_.nx; ++j) {
        x[info_.idx2d(i,j)] = distribution(generator) * 0.005;
      }
    }
    // Set BCs if needed
    if ( info_.bc.compare("dirichlet") == 0) {
      x = apply_dirichlet_bc( x , info_ );
    }
    else if ( info_.bc.compare("mixed_neumann_bottom_dirichlet") == 0 ) {
      x = apply_mixed_bc_neumann_with_bottom_dirichlet( x , info_ );
    }
    else if ( info_.bc.compare("mixed_neumann_top_dirichlet") == 0 ) {
      x = apply_mixed_bc_neumann_with_top_dirichlet( x , info_ );
    }
    else if ( info_.bc.compare("neumann") == 0 ) {
      x = apply_neumann_bc( x , info_ );
    }
  }


void CahnHilliard2DRHS::write_state(const aligned_vector<real> &x , const int idx , const int nx , const int ny , std::string& outdir) const
{
  if ( outdir.back() != '/' )
    outdir += '/';
  std::ofstream out;
  out.open( outdir + "C_" + std::to_string(idx) + ".out" ); 
  out.precision(16);
  
  for (int i = 0; i < nx; ++i){
    for (int j = 0; j < ny; ++j){
      out << x[i * ny + j] << " ";
    }
  }

  out.close();
};
