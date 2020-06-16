//  Python bindings
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/operators.h>
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
		//	getters
		.def("get_rows", &ET::Matrix<double>::get_rows)
		.def("get_cols", &ET::Matrix<double>::get_cols)
		.def("get_name", &ET::Matrix<double>::get_name)
		.def("get_mat", &ET::Matrix<double>::get_mat)
		.def("get_row", &ET::Matrix<double>::get_row)
		.def("get_col", &ET::Matrix<double>::get_col)
		//	setters
		.def("set_name", &ET::Matrix<double>::set_name)
		.def("set_row", &ET::Matrix<double>::set_row)
		.def("set_col", &ET::Matrix<double>::set_col)
    //  operator overloads
		.def(py::self == py::self)
    .def(py::self + py::self)
    .def(py::self += py::self)
    .def(py::self - py::self)
    .def(py::self -= py::self)
    .def(py::self * py::self)
    .def(py::self *= py::self)
    .def(py::self + double())
    //.def(double() + py::self)
    .def(py::self += double())
    .def(py::self - double())
    //.def(double() - py::self)
    .def(py::self -= double())
    .def(py::self * double())
    //.def(double() * py::self)
    .def(py::self *= double())
    .def(py::self / double())
    //.def(double() / py::self)
    .def(py::self /= double())

    .def("__call__", [](const ET::Matrix<double> &self, int i, int j) {
				return self(i,j);
		}, py::is_operator())

		.def("__call__", [](const ET::Matrix<double> &self, int i) {
				return self(i);
		}, py::is_operator())

		//	print functionality
		.def("__repr__", [](const ET::Matrix<double> &m) {
				return m.summary();
		})

		//	various functions
		.def("print", &ET::Matrix<double>::print)
    ;
}
