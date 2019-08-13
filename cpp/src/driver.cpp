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
  
  CHparamsScalar chparams;
  SimInfo info;
  
  // *********  Inputs  ***********
  info.nx               = 128;
  info.ny               = 128;
  info.dx               = 1.0 / info.nx;
  info.dy               = 1.0 / info.ny;
  info.t0               = 0.0;
  info.bc               = "periodic";
  info.rhs_type         = "ch_non_thermal";

  double eps_2          = pow( 0.01 ,2 );
  chparams.eps_2        = eps_2;
  chparams.b            = eps_2 / info.dx / info.dx;
  chparams.u            = eps_2 / info.dx / info.dx;
  chparams.sigma        = eps_2 / info.dx / info.dx / info.dx / info.dx / 200.0;
  chparams.m            = 0.0;
  chparams.sigma_noise  = 0.0;
  
  int n_tsteps        = 25;
  double n_dt         = 300.0;
  // ******************************

  // Compute linear timescales
  double dt_biharm  = (info.dx * info.dx * info.dx * info.dx) / chparams.eps_2;
  double dt_diff    = info.dx * info.dx / chparams.u;
  double dt_lin     = 1.0 / chparams.sigma;

  // Setup temporal checkpointing
  double tf         = n_dt * dt_biharm;
  double dt_check   = tf / n_tsteps;


  std::cout << "Biharmonic timescale dt_biharm = " << dt_biharm << std::endl;
  std::cout << "Diffusion timescale dt_diff = " << dt_diff/dt_biharm << " dt_biharm" << std::endl;
  std::cout << "Linear timescale dt_lin = " << dt_lin/dt_biharm << " dt_biharm" << std::endl;

  //set up timer
  timer run_solver("run_ch_solver");
  
  // Run solver
  for (int i=0; i<n_tsteps; i++) {
    info.t0 = i * dt_check;
    info.tf = (i+1) * dt_check;
    std::cout << "t0 = " << info.t0/dt_biharm << " dt_biharm , tf = " << info.tf/dt_biharm << " dt_biharm" << std::endl;
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
