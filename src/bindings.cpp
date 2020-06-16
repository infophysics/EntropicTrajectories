//  Python bindings
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <lapacke.h>
#include <cblas.h>
#include "vector.h"
#include "matrix.h"




namespace py = pybind11;


PYBIND11_MODULE(etraj, m) {

  py::class_<ET::Vector<double>>(m, "Vector")
    .def(py::init<>())
    .def(py::init<unsigned int>())
    .def(py::init<std::vector<double>>())
    ;

  py::class_<ET::Matrix<double>>(m, "Matrix")
    .def(py::init<>())
    .def(py::init<ET::Matrix<double>>())
    .def(py::init<unsigned int>())
    .def(py::init<std::string, unsigned int>())
    .def(py::init<unsigned int, unsigned int>())
    .def(py::init<std::string, unsigned int, unsigned int>())
    .def(py::init<unsigned int, unsigned int, const double&>())
    .def(py::init<std::string, unsigned int, unsigned int, const double&>())
    .def(py::init<unsigned int, std::vector<double>>())
    .def(py::init<std::string, unsigned int, std::vector<double>>())
    .def(py::init<unsigned int, unsigned int, std::vector<double>>())
    .def(py::init<std::string, unsigned int, unsigned int, std::vector<double>>())
    ;
}
