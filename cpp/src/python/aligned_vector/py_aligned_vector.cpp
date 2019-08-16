#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include <sstream>
#include "aligned_vector.h"

namespace py = pybind11;

PYBIND11_MODULE(aligned_vector, m) {
  py::class_<aligned_double_vector> aligned_vector(m, "aligned_double_vector", py::buffer_protocol());
  
  aligned_vector.def_buffer([](aligned_double_vector& m) -> py::buffer_info {
    return py::buffer_info(
      m.data(),                                /* Pointer to buffer */
      sizeof(double),                          /* Size of one scalar */
      py::format_descriptor<double>::format(), /* Python struct-style format descriptor */
      1,                                       /* Number of dimensions */
      { m.size() },                            /* Buffer dimensions */
      { m.get_allocator().get_alignment() } /* Strides (in bytes) for each index */
    );
  });


  aligned_vector.def(
    "__mult__",
    [](const aligned_double_vector& lhs, const aligned_double_vector& rhs) {
      auto n   = std::min(lhs.size(), rhs.size());
      auto tmp = lhs;
      for(uint64_t i = 0; i < n; ++i){
        tmp.at(i) *= rhs.at(i);
      }
      return tmp;
    }, py::is_operator());
    
    
  aligned_vector.def(
    "__add__",
    [](const aligned_double_vector& lhs, const aligned_double_vector& rhs) {
      auto n   = std::min(lhs.size(), rhs.size());
      auto tmp = lhs;
      for(uint64_t i = 0; i < n; ++i){
        tmp.at(i) += rhs.at(i);
      }
      return tmp;
    }, py::is_operator());
      
      
  aligned_vector.def("__len__", [](const aligned_double_vector& v) { return v.size(); });
  
  
  aligned_vector.def("__iter__", [](aligned_double_vector& v) {
    return py::make_iterator(v.begin(), v.end());
  }, py::keep_alive<0, 1>());
  
  
  auto init_aligned_vec = [](uint64_t size, double def) {
    return aligned_double_vector(size, def);
  };


  aligned_vector.def(py::init(init_aligned_vec), "Initialization",
    py::return_value_policy::automatic_reference, py::arg("size") = 0,
    py::arg("fill") = 0.0);
  
  
  auto init_aligned_vec_from_nparray = [](py::array_t<double, py::array::c_style | py::array::forcecast>& rhs) {
    py::buffer_info rhs_buf = rhs.request();
    //query array
    if (rhs_buf.ndim != 1)
      throw std::runtime_error("Number of dimensions must be one");
    //allocate output
    aligned_double_vector result(rhs_buf.shape[0]);
    double* ptr = static_cast<double*>(rhs_buf.ptr);
    for(uint64_t i = 0; i < result.size(); ++i){
      result[i] = ptr[i];
    }
    return result;
  };
  
  
  aligned_vector.def(py::init(init_aligned_vec_from_nparray), "Initialization",
    py::return_value_policy::automatic_reference, py::arg("array"));
    
    
  auto init_aligned_vec_from_list = [](py::list& rhs) {
    //allocate output
    aligned_double_vector result(rhs.size());
    for(uint64_t i = 0; i < result.size(); ++i){
      result[i] = rhs[i].cast<double>();
    }
    return result;
  };
  
  
  aligned_vector.def(py::init(init_aligned_vec_from_list), "Initialization",
    py::return_value_policy::automatic_reference, py::arg("list"));


  auto aligned_vec_info = [&](aligned_double_vector& self) {
    std::stringstream ss;
    ss << "aligned_vector_double( Alignment = " << std::alignment_of<aligned_double_vector>::value;
    ss << ", size = " << self.size();
    ss << ", data = [";
    for(uint64_t i = 0; i < self.size(); ++i){
      ss << self.at(i);
      if(i + 1 < self.size())
        ss << " ";
    }
    ss << "] )";
    return ss;
  };


  auto aligned_vec_generate = [&](aligned_double_vector& self, uint64_t n, double incr, double init) {
    self.resize(n);
    double val = init;
    std::generate(self.begin(), self.end(), [&val, &incr]() {
      double tmp = val;
      val += incr;
      return tmp;
    });
  };
  
  
  auto aligned_vec_str = [&](aligned_double_vector& self) {
    std::stringstream ss;
    if(self.size()==0) ss << "[]";
    if(self.size()==1) ss << "[" << self.at(0) << "]";
    else{
      ss << "[";
      for(uint64_t i = 0; i < self.size()-1; ++i){
        ss << self.at(i) << " ";
      }
      ss << self.at(self.size()-1) << "]";
    }
    return ss.str();
  };
  
  
  aligned_vector.def("__str__", aligned_vec_str, "Print String");
  aligned_vector.def("__repr__", aligned_vec_info, "Print Representation");
  aligned_vector.def("generate", aligned_vec_generate, "Generate data");
}
