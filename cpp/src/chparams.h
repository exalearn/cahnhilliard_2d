#ifndef __CHPARAMS_H__
#define __CHPARAMS_H__

#include <vector>
#include <algorithm>
#include <string>
#include "allocator.h"

// ************************************************************
// These classes form the main user-interface for the CH solver
// Use them to set all simulation options/parameters
// ************************************************************
 
//switch for doing double or single precision
typedef double real;

static bool abs_compare(int a, int b)
{
  return (std::abs(a) < std::abs(b));
}

class SimInfo
{

  // User-interface for simulation info, grid properties, boundary condition, etc.
  
 public:

  SimInfo() { };
  ~SimInfo() { };
  
  real t0, tf;                // Initial and final simulation times
  int iter = 0;                 // Current simulation iteration
  real dx, dy;                // Spatial discretization
  int nx, ny;                   // Number of grid points in (x,y)
  aligned_vector<real> x;        // Optional, used to specify current state IFF t0 != 0 (MAKE THIS MORE ROBUST)
  std::string bc = "periodic";  // Boundary condition type: "periodic", "neumann", "dirichlet", "mixed_neumann_bottom_dirichlet", "mixed_neumann_top_dirichlet"
  std::string rhs_type = "ch_non_thermal"; // RHS type: "ch_non_thermal", "ch_thermal_no_diffusion", "ch_thermal_with_diffusion"
  real BC_dirichlet_ch;       // Used to specify BC value for dirichlet BC
  std::string outdir = "./";    // Filepath to the output directory

  //safe version (used for indices with offsets):
  inline int idx2d(int i, int j){
    i = (i + ny) % ny;
    j = (j + nx) % nx;
    
    return j * ny + i;
  }
  
  //unsafe version:
  inline int idx2du(int i, int j){
    return j * ny + i;
  }
  
 //private:
  
};

class CHparamsScalar
{

  // User-interface for CH parameter specification where all parameters are scalar-valued
  
 public:

  CHparamsScalar();
  ~CHparamsScalar() { };

  // CH parameters: dc/dt = -eps_2 * \nabla^4( c ) + \nabla^2( u*c^3 - b*c ) - sigma*( c - m ) + sigma_noise*eta 
  real eps_2;
  real b;
  real u;
  real sigma;
  real m;
  real sigma_noise;
  // Thermal dynamics: dT/dt = DT * \nabla^2( T ) + f_T
  real DT;
  real f_T;
  real eps2_min, eps2_max, sigma_min, sigma_max, T_min, T_max; // Limiters on eps_2, sigma, and T
  real T_const; 
  // Polymer parameters
  real L_kuhn;            // Kuhn-statistical length
  real N;                 // Total polymer chain length
  real L_omega;           // Length of physical domain (in 2D, this is \sqrt( area ) )
  real X_min, X_max;      // Min/max values of Flory-Huggins

  inline real compute_stability_limit(real& dx , real& dy){
    real dmin = std::min( dx , dy );
    return 0.5 * dmin * dmin * dmin * dmin / eps_2;
  }
  
  inline real convert_temperature_to_flory_huggins( 
                 const real& T ,
					       const real& T_min ,
					       const real& T_max ,
					       const real& X_min ,
					       const real& X_max ){
                   const real dX_dTinv   = ( X_max  - X_min ) / ( 1.0 / T_min - 1.0 / T_max );  
                   const real dTinv      = 1.0 / T - 1.0 / T_max;
                   const real X          = dX_dTinv * dTinv + X_min;

                   return X;
					       }
                 
  inline real compute_eps2_from_polymer_params(     
                 const real& T ,
					       const real& m ,
					       const real& L_kuhn ,
					       const real& N ){
                   const real m_scaled   = 0.5 * ( 1.0 - m );
                   const real Eps_2      = L_kuhn * L_kuhn / ( 3.0 * m_scaled * (1.0 - m_scaled) * (1.0 - m_scaled) * L_kuhn * L_kuhn * T * N * N );

                   return Eps_2;
					       }
                 
  inline real compute_sigma_from_polymer_params(    
                 const real& T ,
					       const real& m ,
					       const real& L_kuhn ,
					       const real& L_omega ,
					       const real& N ){
                   const real m_scaled   = 0.5 * ( 1.0 - m );
                   const real Sigma      = 36.0 * L_omega * L_omega / ( m_scaled * m_scaled * (1.0 - m_scaled) * (1.0 - m_scaled) * L_kuhn * L_kuhn * T * N * N );

                   return Sigma;
					       }
                 
  void compute_and_set_eps2_and_sigma_from_polymer_params( const real& T );

};

class CHparamsVector
{

  // User-interface for CH parameter specification where all parameters are vector-valued fields
  
 public:

  CHparamsVector() { };
  CHparamsVector(int nx , int ny);
  ~CHparamsVector() { };

  // CH parameters: dc/dt = -eps_2 * \nabla^4( c ) + \nabla^2( u*c^3 - b*c ) - sigma*( c - m ) + sigma_noise*eta
  aligned_vector<real> eps_2;
  aligned_vector<real> b;
  aligned_vector<real> u;
  aligned_vector<real> sigma;
  aligned_vector<real> m;
  real sigma_noise;
  // Thermal dynamics: dT/dt = DT * \nabla^2( T ) + f_T
  aligned_vector<real> DT;
  aligned_vector<real> f_T;
  aligned_vector<real> T_const;
  real eps2_min, eps2_max, sigma_min, sigma_max, T_min, T_max; // Limiters on eps_2, sigma, and T
  // Polymer parameters
  real L_kuhn;        // Kuhn-statistical length
  real N;             // Total polymer chain length
  real L_omega;       // Length of physical domain (in 2D, this is \sqrt( area ) )
  real X_min, X_max;  // Min/max values of Flory-Huggins

  inline real compute_stability_limit(real& dx , real& dy){
    real dmin  = std::min( dx , dy );
    int idx_gmax = std::distance( eps_2.begin() , std::max_element( eps_2.begin() , eps_2.end() , abs_compare ) );
    real gmax  = eps_2[ idx_gmax ];
    return 0.5 * dmin * dmin * dmin * dmin / gmax;
  }
  
  inline real convert_temperature_to_flory_huggins( 
                 const real& T ,
					       const real& T_min ,
					       const real& T_max ,
					       const real& X_min ,
					       const real& X_max ){
                   const real dX_dTinv   = ( X_max  - X_min ) / ( 1.0 / T_min - 1.0 / T_max );  
                   const real dTinv      = 1.0 / T - 1.0 / T_max;
                   const real X          = dX_dTinv * dTinv + X_min;

                   return X;
	}
                 
  inline real compute_eps2_from_polymer_params(     
                 const real& T ,
					       const real& m ,
					       const real& L_kuhn ,
					       const real& N ){
                   const real m_scaled   = 0.5 * ( 1.0 - m );
                   const real Eps_2      = L_kuhn * L_kuhn / ( 3.0 * m_scaled * (1.0 - m_scaled) * T * L_omega * L_omega );

                   return Eps_2;
					       }
                 
  inline real compute_sigma_from_polymer_params(    
                 const real& T ,
					       const real& m ,
					       const real& L_kuhn ,
					       const real& L_omega ,
					       const real& N ){
                   const real m_scaled   = 0.5 * ( 1.0 - m );
                   const real Sigma      = 36.0 * L_omega * L_omega / ( m_scaled * m_scaled * (1.0 - m_scaled) * (1.0 - m_scaled) * L_kuhn * L_kuhn * T * N * N );

                   return Sigma;
					       }
                 
  void compute_and_set_eps2_and_sigma_from_polymer_params( const real T ,
							   SimInfo& info );
  
  
};

#endif
