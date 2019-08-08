#ifndef __CHPARAMS_H__
#define __CHPARAMS_H__

#include <vector>
#include <algorithm>
#include <string>

// ************************************************************
// These classes form the main user-interface for the CH solver
// Use them to set all simulation options/parameters
// ************************************************************
 
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
  
  double t0, tf;                // Initial and final simulation times
  int iter = 0;                 // Current simulation iteration
  double dx, dy;                // Spatial discretization
  int nx, ny;                   // Number of grid points in (x,y)
  std::vector<double> x;        // Optional, used to specify current state IFF t0 != 0 (MAKE THIS MORE ROBUST)
  std::string bc = "periodic";  // Boundary condition type: "periodic", "neumann", "dirichlet", "mixed_neumann_bottom_dirichlet", "mixed_neumann_top_dirichlet"
  std::string rhs_type = "ch_non_thermal"; // RHS type: "ch_non_thermal", "ch_thermal_no_diffusion", "ch_thermal_with_diffusion"
  double BC_dirichlet_ch;       // Used to specify BC value for dirichlet BC
  std::string outdir = "./";    // Filepath to the output directory

  inline int idx2d(int i, int j){
    i = (i + ny) % ny;
    j = (j + nx) % nx;
    
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
  double eps_2;
  double b;
  double u;
  double sigma;
  double m;
  double sigma_noise;
  // Thermal dynamics: dT/dt = DT * \nabla^2( T ) + f_T
  double DT;
  double f_T;
  double eps2_min, eps2_max, sigma_min, sigma_max, T_min, T_max; // Limiters on eps_2, sigma, and T
  double T_const; 
  // Polymer parameters
  double L_kuhn;            // Kuhn-statistical length
  double N;                 // Total polymer chain length
  double L_omega;           // Length of physical domain (in 2D, this is \sqrt( area ) )
  double X_min, X_max;      // Min/max values of Flory-Huggins

  inline double compute_stability_limit(double& dx , double& dy){
    double dmin = std::min( dx , dy );
    return 0.5 * dmin * dmin * dmin * dmin / eps_2;
  }
  
  inline double convert_temperature_to_flory_huggins( 
                 const double& T ,
					       const double& T_min ,
					       const double& T_max ,
					       const double& X_min ,
					       const double& X_max ){
                   const double dX_dTinv   = ( X_max  - X_min ) / ( 1.0 / T_min - 1.0 / T_max );  
                   const double dTinv      = 1.0 / T - 1.0 / T_max;
                   const double X          = dX_dTinv * dTinv + X_min;

                   return X;
					       }
                 
  inline double compute_eps2_from_polymer_params(     
                 const double& T ,
					       const double& m ,
					       const double& L_kuhn ,
					       const double& N ){
                   const double m_scaled   = 0.5 * ( 1.0 - m );
                   const double Eps_2      = L_kuhn * L_kuhn / ( 3.0 * m_scaled * (1.0 - m_scaled) * (1.0 - m_scaled) * L_kuhn * L_kuhn * T * N * N );

                   return Eps_2;
					       }
                 
  inline double compute_sigma_from_polymer_params(    
                 const double& T ,
					       const double& m ,
					       const double& L_kuhn ,
					       const double& L_omega ,
					       const double& N ){
                   const double m_scaled   = 0.5 * ( 1.0 - m );
                   const double Sigma      = 36.0 * L_omega * L_omega / ( m_scaled * m_scaled * (1.0 - m_scaled) * (1.0 - m_scaled) * L_kuhn * L_kuhn * T * N * N );

                   return Sigma;
					       }
                 
  void compute_and_set_eps2_and_sigma_from_polymer_params( const double& T );

};

class CHparamsVector
{

  // User-interface for CH parameter specification where all parameters are vector-valued fields
  
 public:

  CHparamsVector() { };
  CHparamsVector(int nx , int ny);
  ~CHparamsVector() { };

  // CH parameters: dc/dt = -eps_2 * \nabla^4( c ) + \nabla^2( u*c^3 - b*c ) - sigma*( c - m ) + sigma_noise*eta
  std::vector<double> eps_2;
  std::vector<double> b;
  std::vector<double> u;
  std::vector<double> sigma;
  std::vector<double> m;
  double sigma_noise;
  // Thermal dynamics: dT/dt = DT * \nabla^2( T ) + f_T
  std::vector<double> DT;
  std::vector<double> f_T;
  std::vector<double> T_const;
  double eps2_min, eps2_max, sigma_min, sigma_max, T_min, T_max; // Limiters on eps_2, sigma, and T
  // Polymer parameters
  double L_kuhn;        // Kuhn-statistical length
  double N;             // Total polymer chain length
  double L_omega;       // Length of physical domain (in 2D, this is \sqrt( area ) )
  double X_min, X_max;  // Min/max values of Flory-Huggins

  inline double compute_stability_limit(double& dx , double& dy){
    double dmin  = std::min( dx , dy );
    int idx_gmax = std::distance( eps_2.begin() , std::max_element( eps_2.begin() , eps_2.end() , abs_compare ) );
    double gmax  = eps_2[ idx_gmax ];
    return 0.5 * dmin * dmin * dmin * dmin / gmax;
  }
  
  inline double convert_temperature_to_flory_huggins( 
                 const double& T ,
					       const double& T_min ,
					       const double& T_max ,
					       const double& X_min ,
					       const double& X_max ){
                   const double dX_dTinv   = ( X_max  - X_min ) / ( 1.0 / T_min - 1.0 / T_max );  
                   const double dTinv      = 1.0 / T - 1.0 / T_max;
                   const double X          = dX_dTinv * dTinv + X_min;

                   return X;
	}
                 
  inline double compute_eps2_from_polymer_params(     
                 const double& T ,
					       const double& m ,
					       const double& L_kuhn ,
					       const double& N ){
                   const double m_scaled   = 0.5 * ( 1.0 - m );
                   const double Eps_2      = L_kuhn * L_kuhn / ( 3.0 * m_scaled * (1.0 - m_scaled) * T * L_omega * L_omega );

                   return Eps_2;
					       }
                 
  inline double compute_sigma_from_polymer_params(    
                 const double& T ,
					       const double& m ,
					       const double& L_kuhn ,
					       const double& L_omega ,
					       const double& N ){
                   const double m_scaled   = 0.5 * ( 1.0 - m );
                   const double Sigma      = 36.0 * L_omega * L_omega / ( m_scaled * m_scaled * (1.0 - m_scaled) * (1.0 - m_scaled) * L_kuhn * L_kuhn * T * N * N );

                   return Sigma;
					       }
                 
  void compute_and_set_eps2_and_sigma_from_polymer_params( const double T ,
							   SimInfo& info );
  
  
};

#endif
