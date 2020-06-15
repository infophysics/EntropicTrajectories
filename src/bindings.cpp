//  Python bindings
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include "vector.h"




namespace py = pybind11;


PYBIND11_MODULE(etraj, m) {

  py::class_<ET::Vector<double>>(m, "Vector")
    .def(py::init<>())
    //.def(py::init<unsigned int>())
    //.def(py::init<std::vector<double>>())
    ;
}
