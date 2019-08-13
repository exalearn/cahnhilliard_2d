#ifndef __UTILS_CH_H__
#define __UTILS_CH_H__

#include <vector>
#include "chparams.h"
#include <math.h>
#include "profile.h"

void compute_ch_nonlocal(const aligned_vector<double> &c,
			 aligned_vector<double> &dcdt,
			 const double t,
			 CHparamsVector& chpV,
			 SimInfo& info);

void compute_ch_nonlocal_stationary_boundaries(const aligned_vector<double> &c,
					       aligned_vector<double> &dcdt,
					       const double t,
					       CHparamsVector& chpV,
					       SimInfo& info);

void compute_ch_nonlocal_neumannBC(const aligned_vector<double> &c,
				   aligned_vector<double> &dcdt,
				   const double t,
				   CHparamsVector& chpV,
				   SimInfo& info);

void compute_ch_nonlocal_mixedBC_neumann_with_bottom_dirichlet(const aligned_vector<double> &c,
							       aligned_vector<double> &dcdt,
							       const double t,
							       CHparamsVector& chpV,
							       SimInfo& info);

void compute_ch_nonlocal_mixedBC_neumann_with_top_dirichlet(const aligned_vector<double> &c,
							    aligned_vector<double> &dcdt,
							    const double t,
							    CHparamsVector& chpV,
							    SimInfo& info);

aligned_vector<double>& set_boundary_values_to_zero( aligned_vector<double> &dcdt ,
						  SimInfo& info );

aligned_vector<double>& apply_dirichlet_bc( aligned_vector<double>& c ,
					 SimInfo& info );

aligned_vector<double>& apply_neumann_bc( aligned_vector<double>& c ,
				       SimInfo& info );

aligned_vector<double>& apply_mixed_bc_neumann_with_bottom_dirichlet( aligned_vector<double>& c ,
								   SimInfo& info );

aligned_vector<double>& apply_mixed_bc_neumann_with_top_dirichlet( aligned_vector<double>& c ,
								SimInfo& info );

aligned_vector<double>& freeze_corners( aligned_vector<double>& dcdt ,
				     SimInfo& info );

//inline double laplace_component(const int& i ,
//                         const aligned_vector<double>& c ,
//                         const aligned_vector<double>& u ,
//                         const aligned_vector<double>& b ){
//                           return u[i] * (c[i] * c[i] * c[i]) - b[i] * c[i];
//                         }
                         
inline double laplace_component(const double& c,
                         const double& u,
                         const double& b ){
                           return u * (c * c * c) - b * c;
                         }

CHparamsVector compute_chparams_using_temperature( CHparamsVector& chpV0 ,
						   SimInfo& info,
						   aligned_vector<double>& T );

inline double convert_temperature_to_flory_huggins( const CHparamsVector& chpV ,
					     const SimInfo& info,
					     const double T ){
                 const double dX_dTinv   = ( chpV.X_max  - chpV.X_min ) / ( 1.0 / chpV.T_min - 1.0 / chpV.T_max );  
                 const double dTinv      = 1.0 / T - 1.0 / chpV.T_max;
                 const double X          = dX_dTinv * dTinv + chpV.X_min;
                 
                 return X;
					     }

//inline version
void compute_eps2_and_sigma_from_polymer_params( CHparamsVector& chpV ,
							   SimInfo& info,
							   const aligned_vector<double>& T );


#endif
