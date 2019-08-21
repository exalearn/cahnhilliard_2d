#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <omp.h>
#include "run_ch_solver.h"
#include "profile.h"

int main(int argc, char* argv[])
{

  // ***************************************************************
  // Example driver program for 2D modified Cahn-Hilliard
  // No thermal dependence
  // No polymer property dependence
  // No thermal dynamics
  // Scalar coefficients
  // ***************************************************************

  //int tlevel_provided;
  //MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &tlevel_provided);
  /*int ierr = MPI_Init(&argc, &argv);
  if ( ierr != 0 )
  {
    std::cout << "\n";
    std::cout << "ch2d - Fatal error!\n";
    std::cout << "  MPI_Init returned ierr = " << ierr << "\n";
    exit ( 1 );
    }*/
  //init markers
  LIKWID_MARKER_INIT;
  
#pragma omp parallel
  {
    //init thread markers
    LIKWID_MARKER_THREADINIT;
  }
  
  
  LIKWID_MARKER_START("ch2d");
  
  // ********* POLYMER PARAMETERS *********
  real Xmin     = 0.055;
  real Xmax     = 0.5;
  real N         = 1100;
  real L_repeat = (1e-9) * 50.;
  real n_repeat = 15.;
  real L_omega  = n_repeat * L_repeat;
  real L_kuhn   = (1e-9) * 1.75;
    // **************************************
  
  // *********  Inputs  ***********
  //SimInfo
  SimInfo info;
  info.nx               = 128;
  info.ny               = 128;
  info.dx               = 1.0 / info.nx;
  info.dy               = 1.0 / info.ny;
  info.t0               = 0.0;
  info.bc               = "neumann";
  info.rhs_type         = "ch_thermal_no_diffusion";
  
  //CHParams
  CHparamsVector chparams(info.nx, info.ny);
  chparams.sigma_noise  = 0.0;
  chparams.eps2_min     = 0.0;
  chparams.eps2_max     = 1.0;
  chparams.sigma_min    = 0.0;
  chparams.sigma_max    = 1.0e10;
  chparams.T_min        = 0.1;
  chparams.T_max        = 1.0;
  chparams.L_kuhn       = L_kuhn;
  chparams.N            = N;
  chparams.L_omega      = L_omega;
  chparams.X_min        = Xmin;
  chparams.X_max        = Xmax;
  
  int n_tsteps        = 100;
  real n_dt         = 2000.0;
  // ******************************
  
  //fill vectors
  for(unsigned int i=0; i<info.nx*info.ny; i++){
    chparams.b[i] = 1.;
    chparams.u[i] = 1.;
    chparams.m[i] = 0.;
    chparams.T_const[i] = chparams.T_max;
  }
  chparams.compute_and_set_eps2_and_sigma_from_polymer_params( chparams.T_max , info );
  
  real eps_2_max = chparams.eps_2[0];
  real u_max = chparams.u[0];
  real b_max = chparams.b[0];
  real sigma_max = chparams.sigma[0];
  for(unsigned int i=1; i<info.nx*info.ny; i++){
    eps_2_max = std::max(eps_2_max, chparams.eps_2[i]);
    u_max = std::max(u_max, chparams.u[i]);
    b_max = std::max(b_max, chparams.b[i]);
    sigma_max = std::max(sigma_max, chparams.sigma[i]);
  }
  
  // Compute linear timescales
  real dt_biharm  = (info.dx * info.dx * info.dx * info.dx) / eps_2_max;
  real dt_diff    = info.dx * info.dx / std::max(u_max, b_max);
  real dt_lin     = 1.0 / sigma_max;
  real dt_stiff   = std::min(std::min(dt_biharm, dt_diff), dt_lin);

  // Setup temporal checkpointing
  real tf         = n_dt * dt_stiff;
  real dt_check   = tf / n_tsteps;
  
  //temperatures
  aligned_vector<double> T(n_tsteps);
  for(int t=0; t<n_tsteps/4; t++) T[t] = chparams.T_max;
  for(int t=0; t<n_tsteps/4; t++) T[t+n_tsteps/4] = (chparams.T_max*(float(n_tsteps/4.)-1.-t)+chparams.T_min*t)/(float(n_tsteps)/4.-1.);
  for(int t=n_tsteps/4; t<n_tsteps; t++) T[t] = chparams.T_min;

  std::cout << "Biharmonic timescale dt_biharm = " << dt_biharm << std::endl;
  std::cout << "Diffusion timescale dt_diff = " << dt_diff/dt_biharm << " dt_biharm" << std::endl;
  std::cout << "Linear timescale dt_lin = " << dt_lin/dt_biharm << " dt_biharm" << std::endl;

  //set up timer
  timer run_solver("run_ch_solver");
  
  // Run solver
  for (int i=0; i<n_tsteps; i++) {
    info.t0 = i * dt_stiff;
    info.tf = (i+1) * dt_stiff;
    std::cout << "t0 = " << info.t0/dt_lin << " dt_lin , tf = " << info.tf/dt_lin << " dt_lin" << std::endl;
    for(unsigned int k=0; k<info.nx*info.ny; k++){
      chparams.T_const[k] = T[i];
    }
    run_solver.start();
    run_ch_solver(chparams , info);
    run_solver.stop();
  }
  run_solver.print();

  LIKWID_MARKER_STOP("ch2d");

  //close likwid markers
  LIKWID_MARKER_CLOSE;
  //finalize MPI
  //MPI_Finalize();
  
}
