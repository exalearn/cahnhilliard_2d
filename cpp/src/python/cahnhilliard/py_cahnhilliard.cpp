#include <pybind11/stl_bind.h>
#include <pybind11/iostream.h>
#include "chparams.h"
#include "run_ch_solver.h"

namespace py = pybind11;

PYBIND11_MODULE(cahnhilliard, m) {
  
  //global functions
  m.def("run_ch_solver", [](CHparamsVector& chparamv, SimInfo& info){
    py::scoped_ostream_redirect stream(
            std::cout,                               // std::ostream&
            py::module::import("sys").attr("stdout") // Python output
        );
        run_ch_solver(chparamv, info);
  }, py::arg("chparamv"), py::arg("info"));
  
  m.def("run_ch_solver", [](CHparamsScalar& chparams, SimInfo& info){
    py::scoped_ostream_redirect stream(
            std::cout,                               // std::ostream&
            py::module::import("sys").attr("stdout") // Python output
        );
        run_ch_solver(chparams, info);
  }, py::arg("chparams"), py::arg("info"));
  //m.def("run_ch_solver", (void (*)(CHparamsScalar&, SimInfo&)) &run_ch_solver, py::arg("chparams"), py::arg("info"));
  
  //siminfo
  py::class_<SimInfo> sim_info(m, "SimInfo");
  
  //default constructor
  sim_info.def(py::init<>());
  
  //do all the scalar stuff first
  sim_info.def_property("nx", [](const SimInfo& si) { return si.nx; }, 
                                [](SimInfo& si, const int& val) -> void { si.nx = val; });
  
  sim_info.def_property("ny", [](const SimInfo& si) { return si.ny; }, 
                                [](SimInfo& si, const int& val) -> void { si.ny = val; });
  
  sim_info.def_property("t0", [](const SimInfo& si) { return si.t0; }, 
                              [](SimInfo& si, const real& val) -> void { si.t0 = val; });
  
  sim_info.def_property("tf", [](const SimInfo& si) { return si.tf; }, 
                              [](SimInfo& si, const real& val) -> void { si.tf = val; });
  
  sim_info.def_property("iter", [](const SimInfo& si) { return si.iter; }, 
                                [](SimInfo& si, const int& val) -> void { si.iter = val; });
                               
  sim_info.def_property("dx", [](const SimInfo& si) { return si.dx; }, 
                              [](SimInfo& si, const real& val) -> void { si.dx = val; });
  
  sim_info.def_property("dy", [](const SimInfo& si) { return si.dy; }, 
                              [](SimInfo& si, const real& val) -> void { si.dy = val; });
                              
  sim_info.def_property("bc", [](const SimInfo& si) { return si.bc; }, 
                              [](SimInfo& si, const std::string& val) -> void { si.bc = val; });
                              
  sim_info.def_property("rhs_type", [](const SimInfo& si) { return si.rhs_type; }, 
                                    [](SimInfo& si, const std::string& val) -> void { si.rhs_type = val; });
  
  sim_info.def_property("outdir", [](const SimInfo& si) { return si.outdir; }, 
                                  [](SimInfo& si, const std::string& val) -> void { si.outdir = val; });
  
  sim_info.def_property("BC_dirichlet_ch", [](const SimInfo& si) { return si.BC_dirichlet_ch; }, 
                                           [](SimInfo& si, const real& val) -> void { si.BC_dirichlet_ch = val; });
                                           
  sim_info.def_property("x", [](const SimInfo& si) { return si.x; }, 
                             [](SimInfo& si, const aligned_vector<real>& val) -> void { si.x = val; });


  //chparamsvec
  py::class_<CHparamsVector> chparamv(m, "CHparamsVector");
  
  //default constructor
  chparamv.def(py::init<>());
  chparamv.def(py::init<const int&, const int&>());
  
  //do all the scalar stuff first
  chparamv.def_property("sigma_noise", [](const CHparamsVector& chv) { return chv.sigma_noise; }, 
                              [](CHparamsVector& chv, const real& val) -> void { chv.sigma_noise = val; });
                              
  chparamv.def_property("L_kuhn", [](const CHparamsVector& chv) { return chv.L_kuhn; }, 
                                  [](CHparamsVector& chv, const real& val) -> void { chv.L_kuhn = val; });
  
  chparamv.def_property("N", [](const CHparamsVector& chv) { return chv.N; }, 
                             [](CHparamsVector& chv, const real& val) -> void { chv.N = val; });
  
  chparamv.def_property("L_omega", [](const CHparamsVector& chv) { return chv.L_omega; }, 
                                   [](CHparamsVector& chv, const real& val) -> void { chv.L_omega = val; });
  
  chparamv.def_property("X_min", [](const CHparamsVector& chv) { return chv.X_min; }, 
                                 [](CHparamsVector& chv, const real& val) -> void { chv.X_min = val; });
                                 
  chparamv.def_property("X_max", [](const CHparamsVector& chv) { return chv.X_max; }, 
                                 [](CHparamsVector& chv, const real& val) -> void { chv.X_max = val; });
                                 
  chparamv.def_property("eps2_min", [](const CHparamsVector& chv) { return chv.eps2_min; }, 
                                    [](CHparamsVector& chv, const real& val) -> void { chv.eps2_min = val; });
                                    
  chparamv.def_property("eps2_max", [](const CHparamsVector& chv) { return chv.eps2_max; }, 
                                    [](CHparamsVector& chv, const real& val) -> void { chv.eps2_max = val; });
                                    
  chparamv.def_property("sigma_min", [](const CHparamsVector& chv) { return chv.sigma_min; }, 
                                     [](CHparamsVector& chv, const real& val) -> void { chv.sigma_min = val; });
                                     
  chparamv.def_property("sigma_max", [](const CHparamsVector& chv) { return chv.sigma_max; }, 
                                     [](CHparamsVector& chv, const real& val) -> void { chv.sigma_max = val; });
  
  chparamv.def_property("T_min", [](const CHparamsVector& chv) { return chv.T_min; }, 
                                 [](CHparamsVector& chv, const real& val) -> void { chv.T_min = val; });
  
  chparamv.def_property("T_max", [](const CHparamsVector& chv) { return chv.T_max; }, 
                                 [](CHparamsVector& chv, const real& val) -> void { chv.T_max = val; });
  
  //vector properties
  chparamv.def_property("eps_2", [](const CHparamsVector& chv) { return chv.eps_2; }, 
                                 [](CHparamsVector& chv, const aligned_vector<real>& val) -> void { chv.eps_2 = val; });
  
  chparamv.def_property("b", [](const CHparamsVector& chv) { return chv.b; }, 
                             [](CHparamsVector& chv, const aligned_vector<real>& val) -> void { chv.b = val; });
                             
  chparamv.def_property("u", [](const CHparamsVector& chv) { return chv.u; }, 
                             [](CHparamsVector& chv, const aligned_vector<real>& val) -> void { chv.u = val; });
  
  chparamv.def_property("sigma", [](const CHparamsVector& chv) { return chv.sigma; }, 
                                 [](CHparamsVector& chv, const aligned_vector<real>& val) -> void { chv.sigma = val; });
                             
  chparamv.def_property("m", [](const CHparamsVector& chv) { return chv.m; }, 
                             [](CHparamsVector& chv, const aligned_vector<real>& val) -> void { chv.m = val; });
  
  chparamv.def_property("DT", [](const CHparamsVector& chv) { return chv.DT; }, 
                              [](CHparamsVector& chv, const aligned_vector<real>& val) -> void { chv.DT = val; });
  
  chparamv.def_property("f_T", [](const CHparamsVector& chv) { return chv.f_T; }, 
                               [](CHparamsVector& chv, const aligned_vector<real>& val) -> void { chv.f_T = val; });
  
  chparamv.def_property("T_const", [](const CHparamsVector& chv) { return chv.T_const; }, 
                                   [](CHparamsVector& chv, const aligned_vector<real>& val) -> void { chv.T_const = val; });
                                   
  chparamv.def("compute_and_set_eps2_and_sigma_from_polymer_params",
                &CHparamsVector::compute_and_set_eps2_and_sigma_from_polymer_params, py::arg("T"), py::arg("info") );
}