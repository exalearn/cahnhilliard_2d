#include "utils_ch.h"

void compute_ch_nonlocal(const std::vector<double> &c,
			 std::vector<double> &dcdt,
			 const double t,
			 CHparamsVector& chpV,
			 SimInfo& info) {

  // Computes deterministic nonlocal CH dynamics
  // dc/dt = laplacian( u*c^3 - b*c ) - eps_2*biharm(c) - sigma*(c - m)

  // evaluate the second order term, 5 point central stencil
  const double dxn = 1.0 / (info.dx * info.dx);
  const double dyn = 1.0 / (info.dx * info.dx);
  
  //cache blocking
  //const int bsy = 64;
  //const int bsx = 4;
  //int nby = static_cast<int>(ceil(info.ny / bsy));
  //int nbx = static_cast<int>(ceil(info.nx / bsx));
  
  //#pragma omp parallel for collapse(2) firstprivate(dxn, dyn)
  //for (int jb = 0; jb < info.nx; jb += bsx) {
  //  for (int ib = 0; ib < info.ny; ib += bsy) {
  //    
  //    for (int jj = 0; jj < bsx; ++jj) {
  //      for (int ii = 0; ii < bsy; ++ii) {
  //    
  //        //compute global indices
  //        int j = jj + jb;
  //        int i = ii + ib;
  //        
  //        if( (j < info.ny) && (i < info.nx) ){
  
  const double * __restrict__ c_data = c.data();
  double * __restrict__  dcdt_data = dcdt.data();
  double* __restrict__ u_data = chpV.u.data();
  double* __restrict__ b_data = chpV.b.data();
  double* __restrict__ eps_2_data = chpV.eps_2.data();
  double* __restrict__ sigma_data = chpV.sigma.data();
  double* __restrict__ m_data = chpV.m.data();
  
  #pragma omp parallel for firstprivate(dxn, dyn)
  for (int j = 0; j < info.nx; ++j) {
    for (int i = 0; i < info.ny; ++i) {
    
      //get NN components
      double c_i   = c_data[info.idx2du(i    , j    ) ];
      double c_im1 = c_data[info.idx2d (i - 1, j    ) ];
      double c_ip1 = c_data[info.idx2d (i + 1, j    ) ];
      double c_jm1 = c_data[info.idx2d (i    , j - 1) ];
      double c_jp1 = c_data[info.idx2d (i    , j + 1) ];
    
      // evaluate the second order term, 5 point central stencil
      double l_i   = laplace_component( c_i   , u_data[info.idx2du(i    , j    )] , b_data[info.idx2du(i    , j    )] );
      double l_im1 = laplace_component( c_im1 , u_data[info.idx2d (i - 1, j    )] , b_data[info.idx2d (i - 1, j    )] );
      double l_ip1 = laplace_component( c_ip1 , u_data[info.idx2d (i + 1, j    )] , b_data[info.idx2d (i + 1, j    )] );
      double l_jm1 = laplace_component( c_jm1 , u_data[info.idx2d (i    , j - 1)] , b_data[info.idx2d (i    , j - 1)] );
      double l_jp1 = laplace_component( c_jp1 , u_data[info.idx2d (i    , j + 1)] , b_data[info.idx2d (i    , j + 1)] );
    
      //compute derivative
      double dxx = dxn * ( l_jp1 + l_jm1 - 2.0 * l_i );
      double dyy = dyn * ( l_ip1 + l_im1 - 2.0 * l_i );
      double dcdt_temp = dxx + dyy;
    
    
      // evaluate the 4th order term, 9 point central stencil
      //nearest neighbor
      //c_i   = c[info.idx2d(i, j)];
      //c_im1 = c[info.idx2d(i - 1, j)];
      //c_ip1 = c[info.idx2d(i + 1, j)];
      //c_jm1 = c[info.idx2d(i, j - 1)];
      //c_jp1 = c[info.idx2d(i, j + 1)];
      //next to nearest neighbor
      double c_im2 = c_data[info.idx2d(i - 2, j)];
      double c_ip2 = c_data[info.idx2d(i + 2, j)];
      double c_jm2 = c_data[info.idx2d(i, j - 2)];
      double c_jp2 = c_data[info.idx2d(i, j + 2)];
      double c_ul  = c_data[info.idx2d(i-1 , j-1)];
      double c_ur  = c_data[info.idx2d(i-1 , j+1)];
      double c_bl  = c_data[info.idx2d(i+1 , j-1)];
      double c_br  = c_data[info.idx2d(i+1 , j+1)];

      // y-direction u_yyyy
      double dyyyy = dyn * dyn * (c_ip2 - 4.0*c_ip1 + 6.0*c_i - 4.0*c_im1 + c_im2);
      // x-direction u_xxxx
      double dxxxx = dxn * dxn * (c_jp2 - 4.0*c_jp1 + 6.0*c_i - 4.0*c_jm1 + c_jm2);
      // mixed term 2*u_xxyy
      double dxxyy = dxn * dyn * 2.0 * (4.0*c_i - 2.0*(c_im1 + c_ip1 + c_jm1 + c_jp1) + c_ul + c_ur + c_bl + c_br );
      //add together
      dcdt_temp += -eps_2_data[info.idx2du(i,j)] * ( dxxxx + dyyyy + dxxyy );

    
      // evaluate linear term
      dcdt_temp  += -sigma_data[info.idx2du(i,j)] * ( c_i - m_data[info.idx2du(i,j)] );
      
      
      //store result
      dcdt_data[info.idx2du(i, j)] = dcdt_temp;
    }
  }
        //}
      //}
    //}

  //// evaluate the 4th order term, 9 point central stencil
  //# pragma omp parallel for
  //for (int j = 0; j < info.nx; ++j) {
  //  for (int i = 0; i < info.ny; ++i) {
  //    
  //    const double c_i   = c[info.idx2d(i, j)];
  //    const double c_im1 = c[info.idx2d(i - 1, j)];
  //    const double c_ip1 = c[info.idx2d(i + 1, j)];
  //    const double c_im2 = c[info.idx2d(i - 2, j)];
  //    const double c_ip2 = c[info.idx2d(i + 2, j)];
  //    const double c_jm1 = c[info.idx2d(i, j - 1)];
  //    const double c_jp1 = c[info.idx2d(i, j + 1)];
  //    const double c_jm2 = c[info.idx2d(i, j - 2)];
  //    const double c_jp2 = c[info.idx2d(i, j + 2)];
  //    const double c_ul  = c[info.idx2d(i-1 , j-1)];
  //    const double c_ur  = c[info.idx2d(i-1 , j+1)];
  //    const double c_bl  = c[info.idx2d(i+1 , j-1)];
  //    const double c_br  = c[info.idx2d(i+1 , j+1)];
  //
  //    // y-direction u_yyyy
  //    double dyyyy = 1.0 / (info.dy * info.dy * info.dy * info.dy) * 
  //      (c_ip2 - 4.0*c_ip1 + 6.0*c_i - 4.0*c_im1 + c_im2);
  //
  //    // x-direction u_xxxx
  //    double dxxxx = 1.0 / (info.dx * info.dx * info.dx * info.dx) * 
  //      (c_jp2 - 4.0*c_jp1 + 6.0*c_i - 4.0*c_jm1 + c_jm2);
  //
  //    // mixed term 2*u_xxyy
  //    double dxxyy = 1.0 / (info.dx * info.dx * info.dy * info.dy) * 
  //      2 * (4*c_i - 2*(c_im1 + c_ip1 + c_jm1 + c_jp1) + c_ul + c_ur + c_bl + c_br );
  //
  //    dcdt[info.idx2d(i,j)] += -chpV.eps_2[info.idx2d(i,j)] * ( dxxxx + dyyyy + dxxyy );
  //    
  //  }
  //}

  //// evaluate linear term
  //# pragma omp parallel for
  //for (int j = 0; j < info.nx; ++j){
  //  for (int i = 0; i < info.ny; ++i){
  //      
  //    const double c_i        = c[info.idx2d(i, j)];
  //    dcdt[info.idx2d(i,j)]  += -chpV.sigma[info.idx2d(i,j)] * ( c_i - chpV.m[info.idx2d(i,j)] );
  //
  //  }
  //}
}

std::vector<double>& apply_dirichlet_bc( std::vector<double>& c ,
					 SimInfo& info ) {

  # pragma omp parallel for simd
  for (int i = 0; i < info.nx; ++i) {
    
    c[info.idx2du(0, i)]         = info.BC_dirichlet_ch;
    c[info.idx2du(1, i)]         = info.BC_dirichlet_ch;
    c[info.idx2du(info.ny-1, i)] = info.BC_dirichlet_ch;
    c[info.idx2du(info.ny-2, i)] = info.BC_dirichlet_ch;

  }

  # pragma omp parallel for simd
  for (int i = 0; i < info.ny; ++i) {

    c[info.idx2du(i, 0)]         = info.BC_dirichlet_ch;
    c[info.idx2du(i, 1)]         = info.BC_dirichlet_ch;
    c[info.idx2du(i, info.nx-1)] = info.BC_dirichlet_ch;
    c[info.idx2du(i, info.nx-2)] = info.BC_dirichlet_ch;
    
  }

  return c;

}

std::vector<double>& apply_neumann_bc( std::vector<double>& c ,
				       SimInfo& info ) {

  # pragma omp parallel for simd
  for (int i = 0; i < info.nx; ++i) {
    
    c[info.idx2du(0, i)]         = c[info.idx2du(4, i)];
    c[info.idx2du(1, i)]         = c[info.idx2du(3, i)];
    c[info.idx2du(info.ny-1, i)] = c[info.idx2du(info.ny-5, i)];
    c[info.idx2du(info.ny-2, i)] = c[info.idx2du(info.ny-4, i)];

  }
  
  # pragma omp parallel for simd
  for (int i = 0; i < info.ny; ++i) {

    c[info.idx2du(i, 0)]         = c[info.idx2du(i, 4)];
    c[info.idx2du(i, 1)]         = c[info.idx2du(i, 3)];
    c[info.idx2du(i, info.nx-1)] = c[info.idx2du(i, info.nx-5)];
    c[info.idx2du(i, info.nx-2)] = c[info.idx2du(i, info.nx-4)];

  }
  
  return c;

}

std::vector<double>& apply_mixed_bc_neumann_with_top_dirichlet( std::vector<double>& c ,
                                                                SimInfo& info ) {

  c = apply_neumann_bc( c , info );
  
  # pragma omp parallel for simd
  for (int i = 0; i < info.nx; ++i) {

    c[info.idx2du(info.ny-2, i)] = info.BC_dirichlet_ch;
    c[info.idx2du(info.ny-1, i)] = info.BC_dirichlet_ch;

  }

  return c;
  
}

std::vector<double>& apply_mixed_bc_neumann_with_bottom_dirichlet( std::vector<double>& c ,
								   SimInfo& info ) {

  c = apply_neumann_bc( c , info );
  
  # pragma omp parallel for simd
  for (int i = 0; i < info.nx; ++i) {
    
    c[info.idx2du(0, i)] = info.BC_dirichlet_ch;
    c[info.idx2du(1, i)] = info.BC_dirichlet_ch;

  }

  return c;
  
}

std::vector<double>& set_boundary_values_to_zero( std::vector<double> &dcdt ,
						  SimInfo& info ) {

  # pragma omp parallel for simd
  for (int i = 0; i < info.nx; ++i) {
    
    dcdt[info.idx2du(0, i)]         = 0;
    dcdt[info.idx2du(1, i)]         = 0;
    dcdt[info.idx2du(info.ny-1, i)] = 0;
    dcdt[info.idx2du(info.ny-2, i)] = 0;

  }

  # pragma omp parallel for simd
  for (int i = 0; i < info.ny; ++i) {

    dcdt[info.idx2du(i, 0)]         = 0;
    dcdt[info.idx2du(i, 1)]         = 0;
    dcdt[info.idx2du(i, info.nx-1)] = 0;
    dcdt[info.idx2du(i, info.nx-2)] = 0;

  }

  return dcdt;
  
}

void compute_ch_nonlocal_stationary_boundaries(const std::vector<double> &c,
					       std::vector<double> &dcdt,
					       const double t,
					       CHparamsVector& chpV,
					       SimInfo& info) {

  compute_ch_nonlocal(c, dcdt, t, chpV, info);
  dcdt = set_boundary_values_to_zero( dcdt , info );

}

std::vector<double>& freeze_corners( std::vector<double>& dcdt , SimInfo& info ) {

  dcdt[info.idx2du(0,0)] = 0;         dcdt[info.idx2du(0,1)] = 0;         dcdt[info.idx2du(0,info.nx-2)] = 0;         dcdt[info.idx2du(0,info.nx-1)] = 0;
  dcdt[info.idx2du(1,0)] = 0;         dcdt[info.idx2du(1,1)] = 0;         dcdt[info.idx2du(1,info.nx-2)] = 0;         dcdt[info.idx2du(1,info.nx-1)] = 0;
  dcdt[info.idx2du(info.ny-2,0)] = 0; dcdt[info.idx2du(info.ny-2,1)] = 0; dcdt[info.idx2du(info.ny-2,info.nx-2)] = 0; dcdt[info.idx2du(info.ny-2,info.nx-1)] = 0;
  dcdt[info.idx2du(info.ny-1,0)] = 0; dcdt[info.idx2du(info.ny-1,1)] = 0; dcdt[info.idx2du(info.ny-1,info.nx-2)] = 0; dcdt[info.idx2du(info.ny-1,info.nx-1)] = 0;

  return dcdt;
  
}

void compute_ch_nonlocal_neumannBC(const std::vector<double> &c,
                                   std::vector<double> &dcdt,
                                   const double t,
                                   CHparamsVector& chpV,
                                   SimInfo& info) {

  compute_ch_nonlocal(c, dcdt, t, chpV, info);
  dcdt = apply_neumann_bc( dcdt , info );
  dcdt = freeze_corners( dcdt , info ); // 4 corners don't affect dynamics
  
}

void compute_ch_nonlocal_mixedBC_neumann_with_bottom_dirichlet(const std::vector<double> &c,
                                                               std::vector<double> &dcdt,
                                                               const double t,
                                                               CHparamsVector& chpV,
                                                               SimInfo& info) {

  compute_ch_nonlocal(c, dcdt, t, chpV, info);
  dcdt = apply_mixed_bc_neumann_with_bottom_dirichlet( dcdt , info );
  dcdt = freeze_corners( dcdt , info );
  
}

void compute_ch_nonlocal_mixedBC_neumann_with_top_dirichlet(const std::vector<double> &c,
                                                            std::vector<double> &dcdt,
                                                            const double t,
                                                            CHparamsVector& chpV,
                                                            SimInfo& info) {

  compute_ch_nonlocal(c, dcdt, t, chpV, info);
  dcdt = apply_mixed_bc_neumann_with_top_dirichlet( dcdt , info );
  dcdt = freeze_corners( dcdt , info );
  
}

//double laplace_component(int i ,
//                         const std::vector<double>& c ,
//                         const std::vector<double>& u ,
//                         const std::vector<double>& b ) {
//
//  return u[i] * (c[i] * c[i] * c[i]) - b[i] * c[i];
//}

CHparamsVector compute_chparams_using_temperature( CHparamsVector& chpV0,
                                                   SimInfo& info,
                                                   std::vector<double>& T ) {

  CHparamsVector chpV = chpV0;
  double deps2_dT     = ( chpV.eps2_max  - chpV.eps2_min )  / ( chpV.T_max - chpV.T_min );
  double dsigma_dT    = ( chpV.sigma_max - chpV.sigma_min ) / ( chpV.T_max - chpV.T_min );
  
  double * __restrict__ T_data = T.data();
  double * __restrict__ eps_2_data = chpV.eps_2.data();
  double * __restrict__ sigma_data = chpV.sigma.data();
  
  # pragma omp parallel for
  for (int j = 0; j < info.nx; ++j) {
    for (int i = 0; i < info.ny; ++i) {

      const double dT         = T_data[info.idx2du(i, j)] - chpV.T_min;
      const double eps2_fit   = deps2_dT  * dT + chpV.eps2_min;
      const double sigma_fit  = dsigma_dT * dT + chpV.sigma_min;
      eps_2_data[info.idx2du(i, j)] = std::min( std::max( eps2_fit  , chpV.eps2_min )  , chpV.eps2_max );
      sigma_data[info.idx2du(i, j)] = std::min( std::max( sigma_fit , chpV.sigma_min ) , chpV.sigma_max );

    }
  }
  
  return chpV;

}


CHparamsVector compute_eps2_and_sigma_from_polymer_params( CHparamsVector& chpV0,
                                                           SimInfo& info,
                                                           std::vector<double>& T ) {

  CHparamsVector chpV = chpV0;
  
  const double * __restrict__ T_data = T.data();
  double * __restrict__ eps_2_data = chpV.eps_2.data();
  double * __restrict__ sigma_data = chpV.sigma.data();
  double * __restrict__ m_data = chpV.m.data();
  
  # pragma omp parallel for
  for (int j = 0; j < info.nx; ++j) {
    for (int i = 0; i < info.ny; ++i) {
      
      const int idx_ij        = info.idx2du(i, j);
      const double X          = convert_temperature_to_flory_huggins( chpV , info , T_data[idx_ij] );
      const double m_scaled   = 0.5 * ( 1.0 - m_data[idx_ij] );
      const double eps_2      = chpV.L_kuhn * chpV.L_kuhn / ( 3.0 * m_scaled * (1.0 - m_scaled) * X * chpV.L_omega * chpV.L_omega );
      const double sigma      = 36.0 * chpV.L_omega * chpV.L_omega / ( m_scaled * m_scaled * (1.0 - m_scaled) * (1.0 - m_scaled) * chpV.L_kuhn * chpV.L_kuhn * X * chpV.N * chpV.N );
      
      eps_2_data[idx_ij] = std::min( std::max( eps_2 , chpV.eps2_min )  , chpV.eps2_max );
      sigma_data[idx_ij] = std::min( std::max( sigma , chpV.sigma_min ) , chpV.sigma_max );

    }
  }
  
  return chpV;

}

//double convert_temperature_to_flory_huggins( CHparamsVector& chpV,
//                                             SimInfo& info,
//                                             const double T ) {
//
//  const double dX_dTinv   = ( chpV.X_max  - chpV.X_min ) / ( 1.0 / chpV.T_min - 1.0 / chpV.T_max );  
//  const double dTinv      = 1.0 / T - 1.0 / chpV.T_max;
//  const double X          = dX_dTinv * dTinv + chpV.X_min;
//
//  return X;
//
//}
