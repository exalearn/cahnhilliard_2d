#include <vector>
#include <algorithm>
#include <math.h>
#include "chparams.h"


//real CHparamsScalar::compute_stability_limit(real dx , real dy) {
//  real dmin = std::min( dx , dy );
//  return 0.5 * dmin * dmin * dmin * dmin / eps_2;
//};

//real CHparamsScalar::convert_temperature_to_flory_huggins( const real T ,
//                                                             const real T_min ,
//                                                             const real T_max ,
//                                                             const real X_min ,
//                                                             const real X_max ) {
//
//  const real dX_dTinv   = ( X_max  - X_min ) / ( 1.0 / T_min - 1.0 / T_max );  
//  const real dTinv      = 1.0 / T - 1.0 / T_max;
//  const real X          = dX_dTinv * dTinv + X_min;
//
//  return X;
//
//};

//real CHparamsScalar::compute_eps2_from_polymer_params( const real X ,
//                                                         const real m ,
//                                                         const real L_kuhn ,
//                                                         const real N ) {
//  
//  const real m_scaled   = 0.5 * ( 1.0 - m );
//  const real Eps_2      = L_kuhn * L_kuhn / ( 3.0 * m_scaled * (1.0 - m_scaled) * (1.0 - m_scaled) * L_kuhn * L_kuhn * X * N * N );
//
//  return Eps_2;
//};
//
//real CHparamsScalar::compute_sigma_from_polymer_params( const real X ,
//                                                          const real m ,
//                                                          const real L_kuhn ,
//                                                          const real L_omega ,
//                                                          const real N  ) {
//  
//  const real m_scaled   = 0.5 * ( 1.0 - m );
//  const real Sigma      = 36.0 * L_omega * L_omega / ( m_scaled * m_scaled * (1.0 - m_scaled) * (1.0 - m_scaled) * L_kuhn * L_kuhn * X * N * N );
//
//  return Sigma;
//};

void CHparamsScalar::compute_and_set_eps2_and_sigma_from_polymer_params( const real& T ) {

  const real X = convert_temperature_to_flory_huggins( T , T_min , T_max , X_min , X_max );
  eps_2          = compute_eps2_from_polymer_params( X , m , L_kuhn , N );
  sigma          = compute_sigma_from_polymer_params( X , m , L_kuhn , L_omega , N );

};

CHparamsScalar::CHparamsScalar( ) {
  
  // CH parameter defaults
  b            = 1.0;
  u            = 1.0;
  m            = 0.0;
  eps2_min     = 0.0;
  eps2_max     = 1.0;
  sigma_min    = 0.0;
  sigma_max    = pow( 1.0 , 10 );
  sigma_noise  = 0.0;

  // Thermal dynamics defaults
  DT        = 1.0;
  f_T       = 0.0;
  T_min     = 0.1;
  T_max     = 1.0;
  T_const   = 0.5 * ( T_min + T_max );  

  // Polymer physics defaults
  X_min           = 0.055;
  X_max           = 0.5;
  N               = 0.5 * ( 200.0 + 2000.0 );
  real L_repeat = 0.5 * ( 20.0 + 80.0 );   // nanometers
  int n_repeat    = 15;
  L_omega         = n_repeat * L_repeat;
  L_kuhn          = 0.5 * ( 0.5 + 3.0 );     // nanometers

  // Set CH parameters based on polymer defaults
  real T_mid = 0.5 * ( T_min + T_max );
  real X     = convert_temperature_to_flory_huggins(          T_mid , T_min , T_max , X_min , X_max );
  eps_2        = compute_eps2_from_polymer_params(      T_mid , m     , L_kuhn , N );
  sigma        = compute_sigma_from_polymer_params(     T_mid , m     , L_kuhn , L_omega , N );

};

CHparamsVector::CHparamsVector( const int& nx , const int& ny ) {
  
  // CH parameter defaults
  eps_2.resize(   nx * ny );
  sigma.resize(   nx * ny );
  b.resize(       nx * ny ); std::fill( b.begin() , b.end() , 1 );
  u.resize(       nx * ny ); std::fill( u.begin() , u.end() , 1 );
  m.resize(       nx * ny ); std::fill( m.begin() , m.end() , 0 );
  eps2_min     = 0.0;
  eps2_max     = 1.0;
  sigma_min    = 0.0;
  sigma_max    = pow( 1.0 , 10 );
  sigma_noise  = 0.0;

  // Thermal dynamics defaults
  DT.resize(      nx * ny ); std::fill( eps_2.begin() , eps_2.end() , 1 );
  f_T.resize(     nx * ny ); std::fill( eps_2.begin() , eps_2.end() , 0 );
  T_min     = 0.1;
  T_max     = 1.0;
  T_const.resize( nx * ny ); std::fill( eps_2.begin() , eps_2.end() , 0.5 * ( T_min + T_max ) );  

  // Polymer physics defaults
  X_min           = 0.055;
  X_max           = 0.5;
  N               = 0.5 * ( 200.0 + 2000.0 );
  real L_repeat = 0.5 * ( 20.0 + 80.0 );   // nanometers
  int n_repeat    = 15;
  L_omega         = n_repeat * L_repeat;
  L_kuhn          = 0.5 * ( 0.5 + 3.0 );     // nanometers

  // Set CH parameters based on polymer defaults
  real T_mid = 0.5 * ( T_min + T_max );
  real X     = convert_temperature_to_flory_huggins(          T_mid , T_min , T_max , X_min , X_max );
  real eps2_default  = compute_eps2_from_polymer_params(      T_mid , m[0]  , L_kuhn , N );
  real sigma_default = compute_sigma_from_polymer_params(     T_mid , m[0]  , L_kuhn , L_omega , N );
  std::fill( eps_2.begin() , eps_2.end() , eps2_default  );
  std::fill( sigma.begin() , sigma.end() , sigma_default );

};

//real CHparamsVector::compute_stability_limit(real dx , real dy) {
//  real dmin  = std::min( dx , dy );
//  int idx_gmax = std::distance( eps_2.begin() , std::max_element( eps_2.begin() , eps_2.end() , abs_compare ) );
//  real gmax  = eps_2[ idx_gmax ];
//  return 0.5 * dmin * dmin * dmin * dmin / gmax;
//};

//real CHparamsVector::convert_temperature_to_flory_huggins( const real T ,
//                                                             const real T_min ,
//                                                             const real T_max ,
//                                                             const real X_min ,
//                                                             const real X_max ) {
//
//  const real dX_dTinv   = ( X_max  - X_min ) / ( 1.0 / T_min - 1.0 / T_max );  
//  const real dTinv      = 1.0 / T - 1.0 / T_max;
//  const real X          = dX_dTinv * dTinv + X_min;
//
//  return X;
//
//};

//real CHparamsVector::compute_eps2_from_polymer_params( const real X ,
//                                                         const real m ,
//                                                         const real L_kuhn ,
//                                                         const real N ) {
//  
//  const real m_scaled   = 0.5 * ( 1.0 - m );
//  const real Eps_2      = L_kuhn * L_kuhn / ( 3.0 * m_scaled * (1.0 - m_scaled) * X * L_omega * L_omega );
//
//  return Eps_2;
//};

//real CHparamsVector::compute_sigma_from_polymer_params( const real X ,
//                                                          const real m ,
//                                                          const real L_kuhn ,
//                                                          const real L_omega ,
//                                                          const real N  ) {
//  
//  const real m_scaled   = 0.5 * ( 1.0 - m );
//  const real Sigma      = 36.0 * L_omega * L_omega / ( m_scaled * m_scaled * (1.0 - m_scaled) * (1.0 - m_scaled) * L_kuhn * L_kuhn * X * N * N );
//
//  return Sigma;
//};

void CHparamsVector::compute_and_set_eps2_and_sigma_from_polymer_params( const real T ,
                                                                         SimInfo& info ) {
  #pragma omp parallel for
  for (int i = 0; i < info.nx*info.ny; i++ ) {
    const real X = convert_temperature_to_flory_huggins( T , T_min , T_max , X_min , X_max );
    eps_2[i]       = compute_eps2_from_polymer_params( X , m[i] , L_kuhn , N );
    sigma[i]       = compute_sigma_from_polymer_params( X , m[i] , L_kuhn , L_omega , N );
  }

};

