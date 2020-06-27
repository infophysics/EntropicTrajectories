//------------------------------------------------------------------------------
//  bindings.cpp
//  The Entropic Trajectories Framework
//  -----------------------------------
//  Copyright (C) [2020] by [N. Carrara, F. Costa, P. Pessoa]
//  [ncarrara@albany.edu,felipecosta.physics@gmail.com,
//    pedroh.pessoa100@gmail.com]
//
//  Permission to use, copy, modify, and/or distribute this software for any
//  purpose with or without fee is hereby granted.
//
//  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
//  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
//  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
//  SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
//  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
//  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR
//  IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
//------------------------------------------------------------------------------
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

  py::class_<ET::Vector<double>>(m, "Vector", py::buffer_protocol())
    .def(py::init<>())
		.def(py::init<ET::Vector<double>>())
    .def(py::init<uint32_t>())
		.def(py::init<std::string, uint32_t>())
    .def(py::init<std::vector<double>>())
		.def(py::init<std::string, std::vector<double>>())
		.def(py::init<uint32_t, const double&>())
		.def(py::init<std::string, uint32_t, const double&>())
		.def(py::init([](py::buffer const b)
		{
    	py::buffer_info info = b.request();
    	if (info.format != py::format_descriptor<float>::format()
					|| info.ndim != 1)
			{
      	throw std::runtime_error("Incompatible buffer format!");
			}
      auto v = new ET::Vector<double>(info.shape[0]);
      memcpy(v->getVec().data(), info.ptr,
						 sizeof(double) * (size_t) (v->getDim()));
      return v;
    }))
		.def("get_dim", &ET::Vector<double>::getDim)
		.def("get_vec", &ET::Vector<double>::getVec)
		.def("get_name", &ET::Vector<double>::getName)
		.def("set_dim", &ET::Vector<double>::setDim)
		.def("set_vec", &ET::Vector<double>::setVec)
		.def("set_name", &ET::Vector<double>::setName)
		.def("dot", &ET::Vector<double>::dot)
		//	operator overloads
		.def(py::self == py::self)
		.def(py::self != py::self)
    .def(py::self + py::self)
    .def(py::self += py::self)
    .def(py::self - py::self)
    .def(py::self -= py::self)
    .def(py::self + double())
    .def(py::self += double())
    .def(py::self - double())
    .def(py::self -= double())
		.def("__mul__", [](const ET::Vector<double>& v,
											 const ET::Vector<double>& u)
		{
			return v * u;
		}, py::is_operator())
    .def(py::self * double())
    .def(py::self *= double())
    .def(py::self / double())
    .def(py::self /= double())
		.def(-py::self)
		.def(double() + py::self)
		.def(double() - py::self)
		.def(double() * py::self)
		.def(double() / py::self)
		.def("__getitem__", [](const ET::Vector<double> &self, int i)
		{
			if (i < 0 || i >= self.getDim())
			{
				throw py::index_error("Index " + std::to_string(i) +
															" out of bounds for vector with dimension "
															+ std::to_string(self.getDim()) + "!");
			}
			return self(i);
		}, py::is_operator())
		.def("__setitem__", [](ET::Vector<double> &self,
													 int i, const double& val)
		{
			if (i < 0 || i >= self.getDim())
			{
				throw py::index_error("Index " + std::to_string(i) +
														" out of bounds for vector with dimension "
														+ std::to_string(self.getDim()) + "!");
			}
			self(i) = val;
		}, py::is_operator())
		//	print functionality
		.def("__repr__", [](const ET::Vector<double> &vector)
		{
				return vector.summary();
		})
		;
		//--------------------------------------------------------------------------

		//--------------------------------------------------------------------------
		//	Instantiations of matrices
		//--------------------------------------------------------------------------
		//m.def("zeros", (ET::Vector<double> (*)(uint32_t)) &ET::zeros);
		//m.def("zeros", (ET::Vector<double> (*)(uint32_t,uint32_t)) &ET::zeros);
		//m.def("ones", (ET::Vector<double> (*)(uint32_t)) &ET::ones);
		//m.def("ones", (ET::Vector<double> (*)(uint32_t,uint32_t)) &ET::ones);
		//--------------------------------------------------------------------------

		//--------------------------------------------------------------------------
		//	Level 1 BLAS methods
		//--------------------------------------------------------------------------
		m.def("dswap", &ET::DSWAP);
		m.def("dscal", &ET::DSCAL);
		m.def("dcopy", (void (*)(ET::Vector<double>&,
				  ET::Vector<double>&)) &ET::DCOPY);
		m.def("dcopy", (ET::Vector<double> (*)(ET::Vector<double>&)) &ET::DCOPY);
		m.def("daxpy", &ET::DAXPY);
		m.def("ddot",  &ET::DDOT);
		m.def("dnrm2", &ET::DNRM2);
		m.def("dasum", &ET::DASUM);
		m.def("idamax",&ET::IDAMAX);
		m.def("idamin",&ET::IDAMIN);
		//--------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Matrix class
	//----------------------------------------------------------------------------
	py::class_<ET::Matrix<double>>(m, "Matrix", py::buffer_protocol())
    .def(py::init<>())
    .def(py::init<ET::Matrix<double>>())
    .def(py::init<uint32_t>())
    .def(py::init<std::string, uint32_t>())
    .def(py::init<uint32_t, uint32_t>())
    .def(py::init<std::string, uint32_t, uint32_t>())
    .def(py::init<uint32_t, uint32_t, const double&>())
    .def(py::init<std::string, uint32_t, uint32_t, const double&>())
    .def(py::init<uint32_t, std::vector<double>>())
    .def(py::init<std::string, uint32_t, std::vector<double>>())
    .def(py::init<uint32_t, uint32_t, std::vector<double>>())
    .def(py::init<std::string, uint32_t, uint32_t, std::vector<double>>())
		.def(py::init<std::vector<std::vector<double>>>())
		.def(py::init<std::string, std::vector<std::vector<double>>>())
		.def(py::init([](py::buffer const b)
		{
    	py::buffer_info info = b.request();
      if (info.format != py::format_descriptor<float>::format()
					|| info.ndim != 2)
			{
      	throw std::runtime_error("Incompatible buffer format!");
			}
      auto v = new ET::Matrix<double>(info.shape[0], info.shape[1]);
      memcpy(v->getArray().data(), info.ptr,
						 sizeof(double) * (size_t) (v->getNumRows() * v->getNumCols()));
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
				 (uint32_t,std::vector<double>)) &ET::Matrix<double>::setArray)
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
													 std::tuple<uint32_t, uint32_t> ij)
		{
			uint32_t i, j;
			std::tie(i, j) = ij;
			if (i < 0 || i >= self.getNumRows())
			{
				throw py::index_error("Index " + std::to_string(i) +
																" out of bounds for array with "
																+ std::to_string(self.getNumRows())
																+ " rows!");
			}
			if (j < 0 || j >= self.getNumCols())
			{
				throw py::index_error("Index " + std::to_string(j) +
															  " out of bounds for array with "
															  + std::to_string(self.getNumCols())
																+ " columns!");
			}
			return self(i,j);
		}, py::is_operator())
		//	now for the setter of the same type.  To set an element to a Matrix
		//	at index i,j, write - m[i,j] = x.
		.def("__setitem__", [](ET::Matrix<double> &self,
													 std::tuple<uint32_t, uint32_t> ij,
													 const double& val)
		{
			uint32_t i, j;
			std::tie(i, j) = ij;
			if (i < 0 || i >= self.getNumRows())
			{
				throw py::index_error("Index " + std::to_string(i) +
																" out of bounds for array with "
																+ std::to_string(self.getNumRows())
																+ " rows!");
			}
			if (j < 0 || j >= self.getNumCols())
			{
				throw py::index_error("Index " + std::to_string(j) +
															  " out of bounds for array with "
															  + std::to_string(self.getNumCols())
																+ " columns!");
			}
			self(i,j) = val;
		}, py::is_operator())
		//	this will allow the user to get entire rows.
		.def("__getitem__", [](const ET::Matrix<double> &self, int i)
		{
			if (i < 0 || i >= self.getNumRows())
			{
				throw py::index_error("Index " + std::to_string(i) +
															" out of bounds for array with "
															+ std::to_string(self.getNumRows())
															+ " rows!");
			}
			return self.getRow(i);
		}, py::is_operator())
		//	print functionality
		.def("__repr__", [](const ET::Matrix<double> &m)
		{
				return m.summary();
		})
		//	Buffer definition allows one to cast a Matrix as a numpy array
		.def_buffer([](ET::Matrix<double> &m) -> py::buffer_info
		{
      return py::buffer_info(
				m.data(),                               /* Pointer to buffer */
				sizeof(float),                          /* Size of one scalar */
				py::format_descriptor<float>::format(), /* Python struct-style format descriptor */
				2,                                      /* Number of dimensions */
				{ m.getNumRows(), m.getNumCols() },     /* Buffer dimensions */
				{ sizeof(float) * m.getNumCols(),       /* Strides (in bytes) for each index */
					sizeof(float) });
		})
		//--------------------------------------------------------------------------

		//--------------------------------------------------------------------------
		//	Various methods
		//--------------------------------------------------------------------------
		.def("print", &ET::Matrix<double>::print)
		.def("transpose", (ET::Matrix<double> (ET::Matrix<double>::*)())
				 &ET::Matrix<double>::transpose)
		.def("transpose", (void (ET::Matrix<double>::*)(bool))
				 &ET::Matrix<double>::transpose_inplace)
		.def("trace", &ET::Matrix<double>::trace)
		.def("find_singular_values", &ET::Matrix<double>::findSingularValues)
		.def("get_singular_values", &ET::Matrix<double>::getSingularValues)
		//--------------------------------------------------------------------------

		//--------------------------------------------------------------------------
		//	Linear algebra tools
		//--------------------------------------------------------------------------
		.def("inverse", &ET::Matrix<double>::inverse)
		.def("pseudo_inverse", &ET::Matrix<double>::pseudoInverse)
		.def("LU", [](ET::Matrix<double> &self)
		{
			if (self.getNumRows() != self.getNumCols())
			{
				throw py::value_error("Expecting square matrix!");
			}
			return self.LU();
		})
		.def("SVD", &ET::Matrix<double>::SVD)
		//--------------------------------------------------------------------------
		;
		//--------------------------------------------------------------------------

		//--------------------------------------------------------------------------
		//	Instantiations of matrices
		//--------------------------------------------------------------------------
		m.def("identity", (ET::Matrix<double> (*)(uint32_t)) &ET::identity_d);
		//m.def("zeros", (ET::Matrix<double> (*)(uint32_t)) &ET::zeros);
		//m.def("zeros", (ET::Matrix<double> (*)(uint32_t,uint32_t)) &ET::zeros);
		//m.def("ones", (ET::Matrix<double> (*)(uint32_t)) &ET::ones);
		//m.def("ones", (ET::Matrix<double> (*)(uint32_t,uint32_t)) &ET::ones);
		m.def("permutation_matrix", (ET::Matrix<double> (*)(const uint32_t&,
		      const std::vector<uint32_t>)) &ET::permutationMatrix_d);
		//--------------------------------------------------------------------------
		//	Level 2 BLAS functions
		//--------------------------------------------------------------------------
		m.def("dgemv", (ET::Vector<double> (*)(double&,ET::Matrix<double>&,
					ET::Vector<double>&)) &ET::DGEMV);
		m.def("dgemv", (void (*)(double&,ET::Matrix<double>&,
					ET::Vector<double>&, double&, ET::Vector<double>&)) &ET::DGEMV);
		m.def("dger", (ET::Matrix<double> (*)(double&,ET::Vector<double>&,
					ET::Vector<double>&)) &ET::DGER);
		m.def("dger", (void (*)(double&,ET::Vector<double>&,
					ET::Vector<double>&, ET::Matrix<double>&)) &ET::DGER);
		//--------------------------------------------------------------------------

		//--------------------------------------------------------------------------
		//  Level 3 BLAS functions
		//--------------------------------------------------------------------------
		m.def("dgemm", (ET::Matrix<double> (*)(const double&,
			    const ET::Matrix<double>&, const ET::Matrix<double>&))
					&ET::DGEMM);
		m.def("dgemm", (void (*)(const double&, const ET::Matrix<double>&,
			    const ET::Matrix<double>&, const double&,
					ET::Matrix<double>&)) &ET::DGEMM);
		//	DGEMM with no scalar "alpha"
		m.def("dgemm", (ET::Matrix<double> (*)(const ET::Matrix<double>&,
			    const ET::Matrix<double>&)) &ET::DGEMM);
		//--------------------------------------------------------------------------

		//--------------------------------------------------------------------------
		//	LAPACK functions
		//--------------------------------------------------------------------------
		m.def("dgels", (ET::Matrix<double> (*)(const ET::Matrix<double>&,
		      const ET::Matrix<double>&)) &ET::DGELS);
		m.def("dgels", (ET::Vector<double> (*)(const ET::Matrix<double>&,
		      const ET::Vector<double>&)) &ET::DGELS);
		m.def("dgetrf", &ET::DGETRF);
		m.def("dgetrf_l_u", &ET::DGETRF_L_U);
		m.def("dgetrf_lu", &ET::DGETRF_LU);
		m.def("dgetrf_plu", &ET::DGETRF_PLU);
		m.def("dgeqrf", &ET::DGEQRF);
		m.def("dorgqr", &ET::DORGQR);
		m.def("dgeqrf_qr", &ET::DGEQRF_QR);
		m.def("dgesvd", &ET::DGESVD);
		m.def("dgesvd_svd", &ET::DGESVD_SVD);
		m.def("dgesdd", &ET::DGESDD);
		m.def("dgesdd_svd", &ET::DGESDD_SVD);
		//--------------------------------------------------------------------------

		py::class_<ET::Grid<double>>(m, "Grid")
			.def(py::init<>())
			.def(py::init<uint64_t>())
			.def(py::init<std::string, uint64_t>())
			.def(py::init<uint64_t, uint64_t>())
			.def(py::init<std::string, uint64_t, uint64_t>())
			//	getters
			.def("get_dim", &ET::Grid<double>::getDim)
			.def("get_N", &ET::Grid<double>::getN)
			.def("get_grid", &ET::Grid<double>::getGrid)
			.def("get_name", &ET::Grid<double>::getName)
			.def("get_neighbors", (std::vector<std::vector<size_t> >
					 (ET::Grid<double>::*)()) &ET::Grid<double>::getNeighbors)
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
						std::tuple<uint64_t, uint64_t> ij)
			{
				uint32_t i, j;
				std::tie(i, j) = ij;
				if (i < 0 || i >= self.getN())
				{
					throw py::index_error("Index " + std::to_string(i) +
																" out of bounds for array with "
																+ std::to_string(self.getN())
																+ " points!");
				}
				if (j < 0 || j >= self.getDim())
				{
					throw py::index_error("Index " + std::to_string(j) +
															  " out of bounds for array with dimension "
															  + std::to_string(self.getDim()));
				}
				return self(i,j);
			}, py::is_operator())
			.def("__setitem__", [](ET::Grid<double> &self,
						std::tuple<uint32_t, uint32_t> ij, const double& val)
			{
				uint32_t i, j;
				std::tie(i, j) = ij;
				if (i < 0 || i >= self.getN())
				{
					throw py::index_error("Index " + std::to_string(i) +
																" out of bounds for array with "
																+ std::to_string(self.getN())
																+ " points!");
				}
				if (j < 0 || j >= self.getDim())
				{
					throw py::index_error("Index " + std::to_string(j) +
															  " out of bounds for array with dimension "
															  + std::to_string(self.getDim()));
				}
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
