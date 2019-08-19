#include <boost/numeric/odeint.hpp>
#include "stochastic_euler.hpp"
#include "cahnhilliard.h"
#include "cahnhilliard_thermal.h"
#include "cahnhilliard_thermal_nodiffusion.h"
#include "run_ch_solver.h"

real run_ch_solver_non_thermal( CHparamsVector& chparams , SimInfo& info )
{

  // Instantiate rhs
  CahnHilliard2DRHS rhs = CahnHilliard2DRHS( chparams , info );
  
  //aligned_vector<real> x;
  if (info.t0 == 0) {
    if(info.verbosity>=2)
      std::cout << "applying initial conditions: " << info.bc << std::endl;
    rhs.setInitialConditions(info.x);
    int iter = 0;
    
    if (chparams.sigma_noise < 1e-2){
      if(info.verbosity>=2)
        std::cout << "Solving deterministic (noise-free) CH" << std::endl;
    }
    else{
      if(info.verbosity>=2)
        std::cout << "Solving stochastic CH" << std::endl;
    }
  }
  else {
    //x        = info.x;
    int iter = info.iter;
  }
  
  // define adaptive stepper
  typedef boost::numeric::odeint::runge_kutta_cash_karp54<aligned_vector<real>> error_stepper_type;

  // define runge kutta
  typedef boost::numeric::odeint::controlled_runge_kutta<error_stepper_type> controlled_stepper_type;

  controlled_stepper_type controlled_stepper;

  const real stability_limit = chparams.compute_stability_limit(info.dx , info.dy); // just an estimate
  const real res0            = rhs.l2residual(info.x);
  
  //std::cout << "residual at initial condition: " << res0 << std::endl;
  if (info.iter == 0)
    rhs.write_state(info.x,0,info.nx,info.ny,info.outdir);

  if (chparams.sigma_noise < 1e-2) {
    integrate_adaptive(controlled_stepper, rhs, info.x, info.t0, info.tf, static_cast<real>(stability_limit/2.));
    //boost::numeric::odeint::integrate_const(controlled_stepper, rhs, x, info.t0, info.tf, stability_limit/2.);
  }
  else {
    boost::mt19937 rng;
    boost::numeric::odeint::integrate_const( stochastic_euler() ,
		     std::make_pair( rhs , ornstein_stoch( rng , chparams.sigma_noise ) ),
		     info.x , info.t0 , info.tf , static_cast<real>(stability_limit/40.) );
  }
  info.iter += 1;
  
  //residual after solve
  const real res1 =  rhs.l2residual(info.x);
  
  if(info.verbosity>=1)
    std::cout << "iter: " << info.iter << " , t0 = " << info.t0 << " , tf = " << info.tf << ", relative residual: " <<  res1 / res0 << std::endl;
  rhs.write_state(info.x,info.iter,info.nx,info.ny,info.outdir);
  //info.x = x;
  return res1;
};


real run_ch_solver_thermal_no_diffusion( CHparamsVector& chparams , SimInfo& info )
{

  // Instantiate rhs
  CahnHilliard2DRHS_thermal_nodiffusion rhs = CahnHilliard2DRHS_thermal_nodiffusion( chparams , info );
  
  //aligned_vector<real> x;
  if (info.t0 == 0) {
    if(info.verbosity>=2)
      std::cout << "applying initial conditions: " << info.bc << std::endl;
    rhs.setInitialConditions(info.x);
    int iter = 0;
    
    if (chparams.sigma_noise < 1e-2){
      if(info.verbosity>=2)
        std::cout << "Solving deterministic (noise-free) CH" << std::endl;
    }
    else{
      if(info.verbosity>=2)
        std::cout << "Solving stochastic CH" << std::endl;
    }
  }
  else {
    //x        = info.x;
    int iter = info.iter;
  }
  
  // define adaptive stepper
  typedef boost::numeric::odeint::runge_kutta_cash_karp54<aligned_vector<real>> error_stepper_type;

  // define runge kutta
  typedef boost::numeric::odeint::controlled_runge_kutta<error_stepper_type> controlled_stepper_type;

  controlled_stepper_type controlled_stepper;

  const real stability_limit = chparams.compute_stability_limit(info.dx , info.dy); // just an estimate
  const real res0            = rhs.l2residual(info.x);

  //std::cout << "residual at initial condition: " << res0 << std::endl;
  if (info.iter == 0)
    rhs.write_state(info.x,0,info.nx,info.ny,info.outdir);

  if (chparams.sigma_noise < 1e-2) {
    integrate_adaptive(controlled_stepper, rhs, info.x, info.t0, info.tf, static_cast<real>(stability_limit/2.));
    //boost::numeric::odeint::integrate_const(controlled_stepper, rhs, x, info.t0, info.tf, stability_limit/2.);
  }
  else {
    boost::mt19937 rng;
    boost::numeric::odeint::integrate_const( stochastic_euler() ,
		     std::make_pair( rhs , ornstein_stoch( rng , chparams.sigma_noise ) ),
		     info.x , info.t0 , info.tf , static_cast<real>(stability_limit/40.) );
  }
  info.iter += 1;
  
  //residual after solve
  const real res1 =  rhs.l2residual(info.x);
  
  if(info.verbosity>=1)
    std::cout << "iter: " << info.iter << " , t0 = " << info.t0 << " , tf = " << info.tf << ", relative residual: " << res1 / res0 << std::endl;
  rhs.write_state(info.x,info.iter,info.nx,info.ny,info.outdir);
  //info.x = x;
  return res1;
};


real run_ch_solver_thermal_with_diffusion( CHparamsVector& chparams , SimInfo& info )
{

  // Instantiate rhs
  CahnHilliard2DRHS_thermal rhs = CahnHilliard2DRHS_thermal( chparams , info );
  
  //aligned_vector<real> x;
  if (info.t0 == 0) {
    if(info.verbosity>=2)
      std::cout << "applying initial conditions: " << info.bc << std::endl;
    rhs.setInitialConditions(info.x);
    int iter = 0;
    
    if (chparams.sigma_noise < 1e-2){
      if(info.verbosity>=2)
        std::cout << "Solving deterministic (noise-free) CH" << std::endl;
    }
    else{
      if(info.verbosity>=2)
        std::cout << "Solving stochastic CH" << std::endl;
    }
  }
  else {
    //x        = info.x;
    int iter = info.iter;
  }
  
  // define adaptive stepper
  typedef boost::numeric::odeint::runge_kutta_cash_karp54<aligned_vector<real>> error_stepper_type;

  // define runge kutta
  typedef boost::numeric::odeint::controlled_runge_kutta<error_stepper_type> controlled_stepper_type;

  controlled_stepper_type controlled_stepper;

  const real stability_limit = chparams.compute_stability_limit(info.dx , info.dy); // just an estimate
  const real res0            = rhs.l2residual(info.x);

  //std::cout << "residual at initial condition: " << res0 << std::endl;
  if (info.iter == 0)
    rhs.write_state(info.x,0,info.nx,info.ny,info.outdir);

  if (chparams.sigma_noise < 1e-2) {
    integrate_adaptive(controlled_stepper, rhs, info.x, info.t0, info.tf, static_cast<real>(stability_limit/2.) );
    //boost::numeric::odeint::integrate_const(controlled_stepper, rhs, x, info.t0, info.tf, stability_limit/2.);
  }
  else {
    boost::mt19937 rng;
    boost::numeric::odeint::integrate_const( stochastic_euler() ,
		     std::make_pair( rhs , ornstein_stoch( rng , chparams.sigma_noise ) ),
		     info.x , info.t0 , info.tf , static_cast<real>(stability_limit/40.) );
  }
  info.iter += 1;
  
  //residual after solve
  const real res1 =  rhs.l2residual(info.x);
  
  if(info.verbosity>=1)
    std::cout << "iter: " << info.iter << " , t0 = " << info.t0 << " , tf = " << info.tf << ", relative residual: " << res1 / res0 << std::endl;
  rhs.write_state(info.x,info.iter,info.nx,info.ny,info.outdir);
  
  //info.x = x;
  return res1;
};


real run_ch_solver_non_thermal( CHparamsScalar& chparams , SimInfo& info )
{

  // Instantiate rhs
  CahnHilliard2DRHS rhs = CahnHilliard2DRHS( chparams , info );
  
  //aligned_vector<real> x;
  if (info.t0 == 0) {
    if(info.verbosity>=2)
      std::cout << "applying initial conditions: " << info.bc << std::endl;
    rhs.setInitialConditions(info.x);
    int iter = 0;
    
    if (chparams.sigma_noise < 1e-2){
      if(info.verbosity>=2)
        std::cout << "Solving deterministic (noise-free) CH" << std::endl;
    }
    else{
      if(info.verbosity>=2)
        std::cout << "Solving stochastic CH" << std::endl;
    }
  }
  else {
    //x        = info.x;
    int iter = info.iter;
  }
  
  // define adaptive stepper
  typedef boost::numeric::odeint::runge_kutta_cash_karp54<aligned_vector<real>> error_stepper_type;

  // define runge kutta
  typedef boost::numeric::odeint::controlled_runge_kutta<error_stepper_type> controlled_stepper_type;

  controlled_stepper_type controlled_stepper;

  const real stability_limit = chparams.compute_stability_limit(info.dx , info.dy); // just an estimate
  const real res0            = rhs.l2residual(info.x);

  //std::cout << "residual at initial condition: " << res0 << std::endl;
  if (info.iter == 0)
    rhs.write_state(info.x,0,info.nx,info.ny,info.outdir);

  if (chparams.sigma_noise < 1e-2) {
    integrate_adaptive(controlled_stepper, rhs, info.x, info.t0, info.tf, static_cast<real>(stability_limit/2.));
    //boost::numeric::odeint::integrate_const(controlled_stepper, rhs, x, info.t0, info.tf, stability_limit/2.);
  }
  else {
    boost::mt19937 rng;
    boost::numeric::odeint::integrate_const( stochastic_euler() ,
		     std::make_pair( rhs , ornstein_stoch( rng , chparams.sigma_noise ) ),
		     info.x , info.t0 , info.tf , static_cast<real>(stability_limit/40.) );
  }
  info.iter += 1;
  
  //residual after solve
  const real res1 =  rhs.l2residual(info.x);
  
  if(info.verbosity>=1)
    std::cout << "iter: " << info.iter << " , t0 = " << info.t0 << " , tf = " << info.tf << ", relative residual: " << res1 / res0 << std::endl;
  rhs.write_state(info.x,info.iter,info.nx,info.ny,info.outdir);
  //info.x = x;
  
  //residual after solve
  return res1;
};


real run_ch_solver_thermal_no_diffusion( CHparamsScalar& chparams , SimInfo& info )
{

  // Instantiate rhs
  CahnHilliard2DRHS_thermal_nodiffusion rhs = CahnHilliard2DRHS_thermal_nodiffusion( chparams , info );
  
  //aligned_vector<real> x;
  if (info.t0 == 0) {
    if(info.verbosity>=2)
      std::cout << "applying initial conditions: " << info.bc << std::endl;
    rhs.setInitialConditions(info.x);
    int iter = 0;
    
    if (chparams.sigma_noise < 1e-2){
      if(info.verbosity>=2)
        std::cout << "Solving deterministic (noise-free) CH" << std::endl;
    }
    else{
      if(info.verbosity>=2)
        std::cout << "Solving stochastic CH" << std::endl;
    }
  }
  else {
    //x        = info.x;
    int iter = info.iter;
  }
  
  // define adaptive stepper
  typedef boost::numeric::odeint::runge_kutta_cash_karp54<aligned_vector<real>> error_stepper_type;

  // define runge kutta
  typedef boost::numeric::odeint::controlled_runge_kutta<error_stepper_type> controlled_stepper_type;

  controlled_stepper_type controlled_stepper;

  const real stability_limit = chparams.compute_stability_limit(info.dx , info.dy); // just an estimate
  const real res0            = rhs.l2residual(info.x);

  //std::cout << "residual at initial condition: " << res0 << std::endl;
  if (info.iter == 0)
    rhs.write_state(info.x,0,info.nx,info.ny,info.outdir);

  if (chparams.sigma_noise < 1e-2) {
    integrate_adaptive(controlled_stepper, rhs, info.x, info.t0, info.tf, static_cast<real>(stability_limit/2.));
    //boost::numeric::odeint::integrate_const(controlled_stepper, rhs, x, info.t0, info.tf, stability_limit/2.);
  }
  else {
    boost::mt19937 rng;
    boost::numeric::odeint::integrate_const( stochastic_euler() ,
		     std::make_pair( rhs , ornstein_stoch( rng , chparams.sigma_noise ) ),
		     info.x , info.t0 , info.tf , static_cast<real>(stability_limit/40.) );
  }
  info.iter += 1;
  
  //residual after solve
  const real res1 =  rhs.l2residual(info.x);
  
  if(info.verbosity>=1)
    std::cout << "iter: " << info.iter << " , t0 = " << info.t0 << " , tf = " << info.tf << ", relative residual: " << res1 / res0 << std::endl;
  rhs.write_state(info.x,info.iter,info.nx,info.ny,info.outdir);
  //info.x = x;
  return res1;
};


real run_ch_solver_thermal_with_diffusion( CHparamsScalar& chparams , SimInfo& info )
{

  // Instantiate rhs
  CahnHilliard2DRHS_thermal rhs = CahnHilliard2DRHS_thermal( chparams , info );
  
  //aligned_vector<real> x;
  if (info.t0 == 0) {
    if(info.verbosity>=2)
      std::cout << "applying initial conditions: " << info.bc << std::endl;
    rhs.setInitialConditions(info.x);
    int iter = 0;
    
    if (chparams.sigma_noise < 1e-2){
      if(info.verbosity>=2)
        std::cout << "Solving deterministic (noise-free) CH" << std::endl;
    }
    else{
      if(info.verbosity>=2)
        std::cout << "Solving stochastic CH" << std::endl;
    }
  }
  else {
    //x        = info.x;
    int iter = info.iter;
  }
  
  // define adaptive stepper
  typedef boost::numeric::odeint::runge_kutta_cash_karp54<aligned_vector<real>> error_stepper_type;

  // define runge kutta
  typedef boost::numeric::odeint::controlled_runge_kutta<error_stepper_type> controlled_stepper_type;

  controlled_stepper_type controlled_stepper;

  const real stability_limit = chparams.compute_stability_limit(info.dx , info.dy); // just an estimate
  const real res0            = rhs.l2residual(info.x);

  //std::cout << "residual at initial condition: " << res0 << std::endl;
  if (info.iter == 0)
    rhs.write_state(info.x,0,info.nx,info.ny,info.outdir);

  if (chparams.sigma_noise < 1e-2) {
    integrate_adaptive(controlled_stepper, rhs, info.x, info.t0, info.tf, static_cast<real>(stability_limit/2.));
    //boost::numeric::odeint::integrate_const(controlled_stepper, rhs, x, info.t0, info.tf, stability_limit/2.);
  }
  else {
    boost::mt19937 rng;
    boost::numeric::odeint::integrate_const( stochastic_euler() ,
		     std::make_pair( rhs , ornstein_stoch( rng , chparams.sigma_noise ) ),
		     info.x , info.t0 , info.tf , static_cast<real>(stability_limit/40.) );
  }
  info.iter += 1;
  
  //residual after solve
  const real res1 =  rhs.l2residual(info.x);
  
  if(info.verbosity>=1)
    std::cout << "iter: " << info.iter << " , t0 = " << info.t0 << " , tf = " << info.tf << ", relative residual: " << res1 / res0 << std::endl;
  rhs.write_state(info.x,info.iter,info.nx,info.ny,info.outdir);
  //info.x = x;
  
  return res1;
};

real run_ch_solver( CHparamsVector& chparams , SimInfo& info )
{
  
  if      (info.rhs_type.compare("ch_non_thermal") == 0) {
    return run_ch_solver_non_thermal( chparams , info );
  }
  else if (info.rhs_type.compare("ch_thermal_no_diffusion") == 0) {
    return run_ch_solver_thermal_no_diffusion( chparams , info );
  }
  else if (info.rhs_type.compare("ch_thermal_with_diffusion") == 0) {
    return run_ch_solver_thermal_with_diffusion( chparams , info );
  }

};

real run_ch_solver( CHparamsScalar& chparams , SimInfo& info )
{
  
  if      (info.rhs_type.compare("ch_non_thermal") == 0) {
    return run_ch_solver_non_thermal( chparams , info );
  }
  else if (info.rhs_type.compare("ch_thermal_no_diffusion") == 0) {
    return run_ch_solver_thermal_no_diffusion( chparams , info );
  }
  else if (info.rhs_type.compare("ch_thermal_with_diffusion") == 0) {
    return run_ch_solver_thermal_with_diffusion( chparams , info );
  }

};
