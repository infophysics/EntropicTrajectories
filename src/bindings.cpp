//  Python bindings
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>
#include <lapacke.h>
#include <cblas.h>
#include "vector.h"
#include "matrix.h"
#include "grid.h"




namespace py = pybind11;


PYBIND11_MODULE(etraj, m) {

  py::class_<ET::Vector<double>>(m, "Vector")
    .def(py::init<>())
    .def(py::init<unsigned int>())
    .def(py::init<std::vector<double>>())
    ;

  py::class_<ET::Matrix<double>>(m, "Matrix", py::buffer_protocol())
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
		.def(py::init<std::vector<std::vector<double>>>())
		.def(py::init<std::string, std::vector<std::vector<double>>>())
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
		.def(py::self != py::self)
    .def(py::self + py::self)
    .def(py::self += py::self)
    .def(py::self - py::self)
    .def(py::self -= py::self)
    .def(py::self * py::self)
    .def(py::self *= py::self)
    .def(py::self + double())
    .def(py::self += double())
    .def(py::self - double())
    .def(py::self -= double())
    .def(py::self * double())
    .def(py::self *= double())
    .def(py::self / double())
    .def(py::self /= double())
		.def(-py::self)
		.def(double() + py::self)
		.def(double() - py::self)
		.def(double() * py::self)
		.def(double() / py::self)

		//	Matrix array access.  This first method allows the user
		//	to write m[i,j] to get and set elements.
    .def("__getitem__", [](const ET::Matrix<double> &self,
					std::tuple<unsigned int, unsigned int> ij)
		{
				unsigned int i, j;
				std::tie(i, j) = ij;
				if (i < 0 || i >= self.get_rows())
					throw py::index_error("Index " + std::to_string(i) +
																" out of bounds for array with "
																+ std::to_string(self.get_rows()) + " rows!");
				if (j < 0 || j >= self.get_cols())
					throw py::index_error("Index " + std::to_string(j) +
															  " out of bounds for array with "
															  + std::to_string(self.get_cols()) + " columns!");
				return self(i,j);
		}, py::is_operator())

		.def("__setitem__", [](ET::Matrix<double> &self,
					std::tuple<unsigned int, unsigned int> ij, const double& val)
		{
				unsigned int i, j;
				std::tie(i, j) = ij;
				if (i < 0 || i >= self.get_rows())
					throw py::index_error("Index " + std::to_string(i) +
																" out of bounds for array with "
																+ std::to_string(self.get_rows()) + " rows!");
				if (j < 0 || j >= self.get_cols())
					throw py::index_error("Index " + std::to_string(j) +
															  " out of bounds for array with "
															  + std::to_string(self.get_cols()) + " columns!");
				self(i,j) = val;
		}, py::is_operator())

		.def("__getitem__", [](const ET::Matrix<double> &self, int i) {
			if (i < 0 || i >= self.get_rows())
				throw py::index_error("Index " + std::to_string(i) +
															" out of bounds for array with "
															+ std::to_string(self.get_rows()) + " rows!");
			return self.get_row(i);
		}, py::is_operator())

		//	print functionality
		.def("__repr__", [](const ET::Matrix<double> &m) {
				return m.summary();
		})

		//	various functions
		.def("print", &ET::Matrix<double>::print)

		.def_buffer([](ET::Matrix<double> &m) -> py::buffer_info {
        return py::buffer_info(
            & *m.get_mat().begin(),                  /* Pointer to buffer */
            sizeof(double),                          /* Size of one scalar */
            py::format_descriptor<double>::format(), /* Python struct-style format descriptor */
            2,                                       /* Number of dimensions */
            { m.get_rows(), m.get_cols() },          /* Buffer dimensions */
            { sizeof(double) * m.get_cols(),         /* Strides (in bytes) for each index */
              sizeof(double) }
					);
		   })
    ;

		py::class_<ET::Grid<double>>(m, "Grid")
			.def(py::init<>())
			.def(py::init<unsigned int>())
			.def(py::init<std::string, unsigned int>())
			.def(py::init<unsigned int, unsigned int>())
			.def(py::init<std::string, unsigned int, unsigned int>())
			//	getters
			.def("get_dim", &ET::Grid<double>::get_dim)
			.def("get_N", &ET::Grid<double>::get_N)
			.def("get_grid", &ET::Grid<double>::get_grid)
			.def("get_name", &ET::Grid<double>::get_name)
			.def("get_neighbors", &ET::Grid<double>::get_neighbors)
			.def("get_distances", &ET::Grid<double>::get_distances)
			.def("get_neighbors_radius", &ET::Grid<double>::get_neighbors_radius)
			.def("get_distances_radius", &ET::Grid<double>::get_distances_radius)
			.def("set_dim", &ET::Grid<double>::set_dim)
			.def("set_N", &ET::Grid<double>::set_N)
			.def("set_grid", &ET::Grid<double>::set_grid)
			.def("set_name", &ET::Grid<double>::set_name)
			//	access operators
			.def("__getitem__", [](const ET::Grid<double> &self,
						std::tuple<unsigned int, unsigned int> ij)
			{
					unsigned int i, j;
					std::tie(i, j) = ij;
					if (i < 0 || i >= self.get_N())
						throw py::index_error("Index " + std::to_string(i) +
																	" out of bounds for array with "
																	+ std::to_string(self.get_N()) + " points!");
					if (j < 0 || j >= self.get_dim())
						throw py::index_error("Index " + std::to_string(j) +
																  " out of bounds for array with dimension "
																  + std::to_string(self.get_dim()));
					return self(i,j);
			}, py::is_operator())

			.def("__setitem__", [](ET::Grid<double> &self,
						std::tuple<unsigned int, unsigned int> ij, const double& val)
			{
					unsigned int i, j;
					std::tie(i, j) = ij;
					if (i < 0 || i >= self.get_N())
						throw py::index_error("Index " + std::to_string(i) +
																	" out of bounds for array with "
																	+ std::to_string(self.get_N()) + " points!");
					if (j < 0 || j >= self.get_dim())
						throw py::index_error("Index " + std::to_string(j) +
																  " out of bounds for array with dimension "
																  + std::to_string(self.get_dim()));
					self(i,j) = val;
			}, py::is_operator())

			.def("query", &ET::Grid<double>::find_neighbors)
			.def("query_radius", &ET::Grid<double>::radius_search)
			;

}
