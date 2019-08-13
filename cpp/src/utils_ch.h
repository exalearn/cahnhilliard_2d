#ifndef __UTILS_CH_H__
#define __UTILS_CH_H__

#include <vector>
#include "chparams.h"
#include <math.h>
#include "profile.h"

void compute_ch_nonlocal(const aligned_vector<real> &c,
			 aligned_vector<real> &dcdt,
			 const real t,
			 CHparamsVector& chpV,
			 SimInfo& info);

void compute_ch_nonlocal_stationary_boundaries(const aligned_vector<real> &c,
					       aligned_vector<real> &dcdt,
					       const real t,
					       CHparamsVector& chpV,
					       SimInfo& info);

void compute_ch_nonlocal_neumannBC(const aligned_vector<real> &c,
				   aligned_vector<real> &dcdt,
				   const real t,
				   CHparamsVector& chpV,
				   SimInfo& info);

void compute_ch_nonlocal_mixedBC_neumann_with_bottom_dirichlet(const aligned_vector<real> &c,
							       aligned_vector<real> &dcdt,
							       const real t,
							       CHparamsVector& chpV,
							       SimInfo& info);

void compute_ch_nonlocal_mixedBC_neumann_with_top_dirichlet(const aligned_vector<real> &c,
							    aligned_vector<real> &dcdt,
							    const real t,
							    CHparamsVector& chpV,
							    SimInfo& info);

aligned_vector<real>& set_boundary_values_to_zero( aligned_vector<real> &dcdt ,
						  SimInfo& info );

aligned_vector<real>& apply_dirichlet_bc( aligned_vector<real>& c ,
					 SimInfo& info );

aligned_vector<real>& apply_neumann_bc( aligned_vector<real>& c ,
				       SimInfo& info );

aligned_vector<real>& apply_mixed_bc_neumann_with_bottom_dirichlet( aligned_vector<real>& c ,
								   SimInfo& info );

aligned_vector<real>& apply_mixed_bc_neumann_with_top_dirichlet( aligned_vector<real>& c ,
								SimInfo& info );

aligned_vector<real>& freeze_corners( aligned_vector<real>& dcdt ,
				     SimInfo& info );

//inline real laplace_component(const int& i ,
//                         const aligned_vector<real>& c ,
//                         const aligned_vector<real>& u ,
//                         const aligned_vector<real>& b ){
//                           return u[i] * (c[i] * c[i] * c[i]) - b[i] * c[i];
//                         }
                         
inline real laplace_component(const real& c,
                         const real& u,
                         const real& b ){
                           return u * (c * c * c) - b * c;
                         }

void compute_chparams_using_temperature( CHparamsVector& chpV0 ,
						   SimInfo& info,
						   aligned_vector<real>& T );

inline real convert_temperature_to_flory_huggins( const CHparamsVector& chpV ,
					     const SimInfo& info,
					     const real T ){
                 const real dX_dTinv   = ( chpV.X_max  - chpV.X_min ) / ( 1.0 / chpV.T_min - 1.0 / chpV.T_max );  
                 const real dTinv      = 1.0 / T - 1.0 / chpV.T_max;
                 const real X          = dX_dTinv * dTinv + chpV.X_min;
                 
                 return X;
					     }

//inline version
void compute_eps2_and_sigma_from_polymer_params( CHparamsVector& chpV ,
							   SimInfo& info,
							   const aligned_vector<real>& T );


#endif
