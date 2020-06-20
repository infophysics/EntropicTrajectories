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
#include "scalar.h"
#include "utils.h"
#include "approximator.h"




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
		.def(py::init([](py::buffer const b) {
            py::buffer_info info = b.request();
            if (info.format != py::format_descriptor<float>::format() || info.ndim != 2)
                throw std::runtime_error("Incompatible buffer format!");

            auto v = new ET::Matrix<double>(info.shape[0], info.shape[1]);
            memcpy(v->getArray().data(), info.ptr, sizeof(double) * (size_t) (v->getNumRows() * v->getNumCols()));
            return v;
        }))
		//	getters
		.def("get_num_rows", &ET::Matrix<double>::getNumRows)
		.def("get_num_cols", &ET::Matrix<double>::getNumCols)
		.def("get_name", &ET::Matrix<double>::getName)
		.def("get_array", &ET::Matrix<double>::getArray)
		.def("get_row", &ET::Matrix<double>::getRow)
		.def("get_col", &ET::Matrix<double>::getCol)
		//	setters
		.def("set_name", &ET::Matrix<double>::setName)
		.def("set_row", &ET::Matrix<double>::setRow)
		.def("set_col", &ET::Matrix<double>::setCol)
		.def("set_array", (void (ET::Matrix<double>::*)
							(uint64_t,std::vector<double>)) &ET::Matrix<double>::setArray)
		.def("set_array", (void (ET::Matrix<double>::*)
							(std::vector<std::vector<double>>)) &ET::Matrix<double>::setArray)

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
		//	to write x = m[i,j] to get elements.
    .def("__getitem__", [](const ET::Matrix<double> &self,
					std::tuple<unsigned int, unsigned int> ij)
		{
				unsigned int i, j;
				std::tie(i, j) = ij;
				if (i < 0 || i >= self.getNumRows())
					throw py::index_error("Index " + std::to_string(i) +
																" out of bounds for array with "
																+ std::to_string(self.getNumRows()) + " rows!");
				if (j < 0 || j >= self.getNumCols())
					throw py::index_error("Index " + std::to_string(j) +
															  " out of bounds for array with "
															  + std::to_string(self.getNumCols()) + " columns!");
				return self(i,j);
		}, py::is_operator())
		//	now for the setter of the same type.  To set an element to a Matrix
		//	at index i,j, write - m[i,j] = x.
		.def("__setitem__", [](ET::Matrix<double> &self,
					std::tuple<unsigned int, unsigned int> ij, const double& val)
		{
				unsigned int i, j;
				std::tie(i, j) = ij;
				if (i < 0 || i >= self.getNumRows())
					throw py::index_error("Index " + std::to_string(i) +
																" out of bounds for array with "
																+ std::to_string(self.getNumRows()) + " rows!");
				if (j < 0 || j >= self.getNumCols())
					throw py::index_error("Index " + std::to_string(j) +
															  " out of bounds for array with "
															  + std::to_string(self.getNumCols()) + " columns!");
				self(i,j) = val;
		}, py::is_operator())
		//	this will allow the user to get entire rows.
		.def("__getitem__", [](const ET::Matrix<double> &self, int i) {
			if (i < 0 || i >= self.getNumRows())
				throw py::index_error("Index " + std::to_string(i) +
															" out of bounds for array with "
															+ std::to_string(self.getNumRows()) + " rows!");
			return self.getRow(i);
		}, py::is_operator())

		//	print functionality
		.def("__repr__", [](const ET::Matrix<double> &m) {
				return m.summary();
		})

		//	various functions
		.def("print", &ET::Matrix<double>::print)

		.def_buffer([](ET::Matrix<double> &m) -> py::buffer_info {
        return py::buffer_info(
					m.data(),                               /* Pointer to buffer */
					sizeof(float),                          /* Size of one scalar */
					py::format_descriptor<float>::format(), /* Python struct-style format descriptor */
					2,                                      /* Number of dimensions */
					{ m.getNumRows(), m.getNumCols() },                 /* Buffer dimensions */
					{ sizeof(float) * m.getNumCols(),             /* Strides (in bytes) for each index */
						sizeof(float) }
					);
		   })

		.def("inverse", &ET::Matrix<double>::inverse)
		.def("permutation_matrix", &ET::Matrix<double>::permutationMatrix)
		.def("LU", [](ET::Matrix<double> &self) {
			if (self.getNumRows() != self.getNumCols())
				throw py::value_error("Expecting square matrix!");
			return self.LU();
		})
		.def("T", (ET::Matrix<double> (ET::Matrix<double>::*)()) &ET::Matrix<double>::transpose)
		.def("T", (void  (ET::Matrix<double>::*)(bool)) &ET::Matrix<double>::transpose)
		.def("find_singular_values", &ET::Matrix<double>::findSingularValues)
		.def("get_singular_values", &ET::Matrix<double>::getSingularValues)
		.def("SVD", &ET::Matrix<double>::SVD)
    ;

		py::class_<ET::Grid<double>>(m, "Grid")
			.def(py::init<>())
			.def(py::init<unsigned int>())
			.def(py::init<std::string, unsigned int>())
			.def(py::init<unsigned int, unsigned int>())
			.def(py::init<std::string, unsigned int, unsigned int>())
			//	getters
			.def("get_dim", &ET::Grid<double>::getDim)
			.def("get_N", &ET::Grid<double>::getN)
			.def("get_grid", &ET::Grid<double>::getGrid)
			.def("get_name", &ET::Grid<double>::getName)
			.def("get_neighbors", (std::vector<std::vector<size_t> > (ET::Grid<double>::*)
								()) &ET::Grid<double>::getNeighbors)
			.def("get_neighbors", (std::vector<size_t>* (ET::Grid<double>::*)
								(uint64_t)) &ET::Grid<double>::getNeighbors)
			.def("get_distances", &ET::Grid<double>::getDistances)
			.def("get_neighbors_radius", &ET::Grid<double>::getNeighborsRadius)
			.def("get_distances_radius", &ET::Grid<double>::getDistancesRadius)
			.def("set_dim", &ET::Grid<double>::setDim)
			.def("set_N", &ET::Grid<double>::setN)
			.def("set_grid", &ET::Grid<double>::setGrid)
			.def("set_name", &ET::Grid<double>::setName)
			//	access operators
			.def("__getitem__", [](const ET::Grid<double> &self,
						std::tuple<unsigned int, unsigned int> ij)
			{
					unsigned int i, j;
					std::tie(i, j) = ij;
					if (i < 0 || i >= self.getN())
						throw py::index_error("Index " + std::to_string(i) +
																	" out of bounds for array with "
																	+ std::to_string(self.getN()) + " points!");
					if (j < 0 || j >= self.getDim())
						throw py::index_error("Index " + std::to_string(j) +
																  " out of bounds for array with dimension "
																  + std::to_string(self.getDim()));
					return self(i,j);
			}, py::is_operator())

			.def("__setitem__", [](ET::Grid<double> &self,
						std::tuple<unsigned int, unsigned int> ij, const double& val)
			{
					unsigned int i, j;
					std::tie(i, j) = ij;
					if (i < 0 || i >= self.getN())
						throw py::index_error("Index " + std::to_string(i) +
																	" out of bounds for array with "
																	+ std::to_string(self.getN()) + " points!");
					if (j < 0 || j >= self.getDim())
						throw py::index_error("Index " + std::to_string(j) +
																  " out of bounds for array with dimension "
																  + std::to_string(self.getDim()));
					self(i,j) = val;
			}, py::is_operator())

			.def("query_neighbors", &ET::Grid<double>::queryNeighbors)
			.def("query_radius", &ET::Grid<double>::queryRadius)
			;

		py::class_<ET::ScalarField<double>>(m, "ScalarField")
			.def(py::init<>())
			.def(py::init<ET::Grid<double>*,std::vector<double>>())
			.def("set_derivative", &ET::ScalarField<double>::setDerivative)
			.def("get_approximator", &ET::ScalarField<double>::getApproximator)
			;

		py::class_<ET::Approximator<double>>(m, "Approximator")
			.def(py::init<>())
			.def("set_derivative", &ET::Approximator<double>::setDerivative)
			.def("set_k", &ET::Approximator<double>::set_k)
			.def("set_n", &ET::Approximator<double>::set_n)
			.def("gradient", &ET::Approximator<double>::gradient)
			.def("gradient_mls", &ET::Approximator<double>::gradientMLS)
			.def("construct_B_matrix", &ET::Approximator<double>::constructBMatrix)
			;

		//	various functions
		m.def("monomial_n", &ET::monomial_n, py::return_value_policy::reference);
		m.def("taylor_polynomial", &ET::taylorPolynomial);
		m.def("taylor_monomial_expansion", &ET::taylorMonomialExpansion);
		m.def("taylor_monomial_factors", &ET::taylorMonomialFactors);
}
