#ifndef __UTILS_CH_H__
#define __UTILS_CH_H__

#include <vector>
#include "chparams.h"

void compute_ch_nonlocal(const std::vector<double> &c,
			 std::vector<double> &dcdt,
			 const double t,
			 CHparamsVector& chpV,
			 SimInfo& info);

void compute_ch_nonlocal_stationary_boundaries(const std::vector<double> &c,
					       std::vector<double> &dcdt,
					       const double t,
					       CHparamsVector& chpV,
					       SimInfo& info);

void compute_ch_nonlocal_neumannBC(const std::vector<double> &c,
				   std::vector<double> &dcdt,
				   const double t,
				   CHparamsVector& chpV,
				   SimInfo& info);

void compute_ch_nonlocal_mixedBC_neumann_with_bottom_dirichlet(const std::vector<double> &c,
							       std::vector<double> &dcdt,
							       const double t,
							       CHparamsVector& chpV,
							       SimInfo& info);

void compute_ch_nonlocal_mixedBC_neumann_with_top_dirichlet(const std::vector<double> &c,
							    std::vector<double> &dcdt,
							    const double t,
							    CHparamsVector& chpV,
							    SimInfo& info);

std::vector<double>& set_boundary_values_to_zero( std::vector<double> &dcdt ,
						  SimInfo& info );

std::vector<double>& apply_dirichlet_bc( std::vector<double>& c ,
					 SimInfo& info );

std::vector<double>& apply_neumann_bc( std::vector<double>& c ,
				       SimInfo& info );

std::vector<double>& apply_mixed_bc_neumann_with_bottom_dirichlet( std::vector<double>& c ,
								   SimInfo& info );

std::vector<double>& apply_mixed_bc_neumann_with_top_dirichlet( std::vector<double>& c ,
								SimInfo& info );

std::vector<double>& freeze_corners( std::vector<double>& dcdt ,
				     SimInfo& info );

inline double laplace_component(const int& i ,
                         const std::vector<double>& c ,
                         const std::vector<double>& u ,
                         const std::vector<double>& b ){
                           return u[i] * (c[i] * c[i] * c[i]) - b[i] * c[i];
                         }

inline CHparamsVector compute_chparams_using_temperature( CHparamsVector& chpV0 ,
						   SimInfo& info,
						   std::vector<double>& T ){
                 const double dX_dTinv   = ( chpV.X_max  - chpV.X_min ) / ( 1.0 / chpV.T_min - 1.0 / chpV.T_max );  
                 const double dTinv      = 1.0 / T - 1.0 / chpV.T_max;
                 const double X          = dX_dTinv * dTinv + chpV.X_min;

                 return X;
						   }

double convert_temperature_to_flory_huggins( CHparamsVector& chpV ,
					     SimInfo& info,
					     const double T );

CHparamsVector compute_eps2_and_sigma_from_polymer_params( CHparamsVector& chpV0 ,
							   SimInfo& info,
							   std::vector<double>& T );



#endif
