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
#include <memory>
#include <cblas.h>
#include "linalg/vector.h"
#include "linalg/matrix.h"
#include "ugrid.h"
#include "scalarfield.h"
#include "utils.h"
#include "approximator.h"

namespace py = pybind11;

PYBIND11_MODULE(etraj, m) {

  py::class_<ET::Vector<double>, std::shared_ptr<ET::Vector<double>>>(m, "Vector", py::buffer_protocol())
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
		.def("get_dim", &ET::Vector<double>::getDim, py::return_value_policy::reference)
		.def("get_vec", &ET::Vector<double>::getVec, py::return_value_policy::reference)
		.def("get_name", &ET::Vector<double>::getName, py::return_value_policy::reference)
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
	py::class_<ET::Matrix<double>,std::shared_ptr<ET::Matrix<double>>>(m, "Matrix", py::buffer_protocol())
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
		.def("get_num_rows", &ET::Matrix<double>::getNumRows, py::return_value_policy::reference)
		.def("get_num_cols", &ET::Matrix<double>::getNumCols, py::return_value_policy::reference)
		.def("get_name", &ET::Matrix<double>::getName, py::return_value_policy::reference)
		.def("get_array", &ET::Matrix<double>::getArray, py::return_value_policy::reference)
		.def("get_row", &ET::Matrix<double>::getRow, py::return_value_policy::reference)
		.def("get_col", &ET::Matrix<double>::getCol, py::return_value_policy::reference)
		.def("get_info", &ET::Matrix<double>::getInfo, py::return_value_policy::reference)
		.def("get_flag", &ET::Matrix<double>::getFlag, py::return_value_policy::reference)
		.def("get_rank", &ET::Matrix<double>::getRank, py::return_value_policy::reference)
		.def("get_singular_values", &ET::Matrix<double>::getSingularValues, py::return_value_policy::reference)
		//	setters
		.def("set_name", &ET::Matrix<double>::setName)
		.def("set_row", &ET::Matrix<double>::setRow)
		.def("set_col", &ET::Matrix<double>::setCol)
		.def("set_array", (void (ET::Matrix<double>::*)
				 (uint32_t,std::vector<double>)) &ET::Matrix<double>::setArray)
		.def("set_array", (void (ET::Matrix<double>::*)
				 (std::vector<std::vector<double>>)) &ET::Matrix<double>::setArray)
		.def("set_singular_values", &ET::Matrix<double>::setSingularValues)
		.def("set_info", &ET::Matrix<double>::setInfo)
		.def("set_flag", &ET::Matrix<double>::setFlag)
		.def("set_rank", &ET::Matrix<double>::setRank)
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
			    const ET::Matrix<double>&)) &ET::DGEMM, py::keep_alive<1, 2>());
		//--------------------------------------------------------------------------

		//--------------------------------------------------------------------------
		//	LAPACK functions
		//--------------------------------------------------------------------------
		m.def("dgels", (ET::Matrix<double> (*)(const ET::Matrix<double>&,
		      const ET::Matrix<double>&)) &ET::DGELS, py::keep_alive<0, 1>());
		m.def("dgels", (ET::Vector<double> (*)(const ET::Matrix<double>&,
		      const ET::Vector<double>&)) &ET::DGELS, py::keep_alive<0, 1>());
		m.def("dgelsy", (ET::Matrix<double> (*)(const ET::Matrix<double>&,
		      const ET::Matrix<double>&)) &ET::DGELSY);
		m.def("dgelsy", (ET::Vector<double> (*)(const ET::Matrix<double>&,
		      const ET::Vector<double>&)) &ET::DGELSY);
		m.def("dgelsd", (ET::Matrix<double> (*)(const ET::Matrix<double>&,
		      const ET::Matrix<double>&)) &ET::DGELSD);
		m.def("dgelsd", (ET::Vector<double> (*)(const ET::Matrix<double>&,
		      const ET::Vector<double>&)) &ET::DGELSD);
		m.def("dgelss", (ET::Matrix<double> (*)(const ET::Matrix<double>&,
		      const ET::Matrix<double>&)) &ET::DGELSS);
		m.def("dgelss", (ET::Vector<double> (*)(const ET::Matrix<double>&,
		      const ET::Vector<double>&)) &ET::DGELSS);
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

		py::class_<ET::UGrid<double>, std::shared_ptr<ET::UGrid<double>>>(m, "UGrid")
			.def(py::init<>())
			.def(py::init<uint64_t>())
			.def(py::init<std::string, uint64_t>())
			.def(py::init<uint64_t, uint64_t>())
			.def(py::init<std::string, uint64_t, uint64_t>())
			.def(py::init<std::vector<double>>(), py::keep_alive<1, 2>())
			.def(py::init<std::vector<std::vector<double>>>(), py::keep_alive<1, 2>())
			//	getters
			.def("get_dim", &ET::UGrid<double>::getDim, py::return_value_policy::reference)
			.def("get_N", &ET::UGrid<double>::getN, py::return_value_policy::reference)
			.def("get_grid", &ET::UGrid<double>::getUGrid, py::return_value_policy::reference)
			.def("get_name", &ET::UGrid<double>::getName, py::return_value_policy::reference)
			.def("get_neighbors", (std::vector<std::vector<size_t>>
					 (ET::UGrid<double>::*)()) &ET::UGrid<double>::getNeighbors, py::return_value_policy::reference)
			.def("get_neighbors", (std::vector<size_t> (ET::UGrid<double>::*)
					 (uint64_t)) &ET::UGrid<double>::getNeighbors, py::return_value_policy::reference)
			.def("get_distances", &ET::UGrid<double>::getDistances, py::return_value_policy::reference)
			.def("get_neighbors_radius", &ET::UGrid<double>::getNeighborsRadius, py::return_value_policy::reference)
			.def("get_distances_radius", &ET::UGrid<double>::getDistancesRadius, py::return_value_policy::reference)
			.def("set_dim", &ET::UGrid<double>::setDim)
			.def("set_N", &ET::UGrid<double>::setN)
			.def("set_grid", &ET::UGrid<double>::setUGrid)
			.def("set_name", &ET::UGrid<double>::setName)
			//	access operators
			.def("__getitem__", [](const ET::UGrid<double> &self,
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
			.def("__setitem__", [](ET::UGrid<double> &self,
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
			.def("__getitem__", [](const ET::UGrid<double> &self, int i)
			{
				if (i < 0 || i >= self.getN())
				{
					throw py::index_error("Index " + std::to_string(i) +
																" out of bounds for array with "
																+ std::to_string(self.getN())
																+ " rows!");
				}
				return self.getPoint(i);
			}, py::is_operator())
			.def("query_neighbors", &ET::UGrid<double>::queryNeighbors)
			.def("query_radius", &ET::UGrid<double>::queryRadius)
			;

		py::class_<ET::ScalarField<double>, std::shared_ptr<ET::ScalarField<double>>>(m, "ScalarField")
			.def(py::init<>())
			.def(py::init<std::shared_ptr<ET::UGrid<double>>>())
			.def(py::init<std::string,std::shared_ptr<ET::UGrid<double>>>())
			.def(py::init<std::shared_ptr<ET::UGrid<double>>,std::vector<double>>())
			.def(py::init<std::string,std::shared_ptr<ET::UGrid<double>>,std::vector<double>>())
			//	Getters and Setters
			.def("get_field", &ET::ScalarField<double>::getField, py::return_value_policy::reference)
			.def("access_field", &ET::ScalarField<double>::accessField, py::return_value_policy::reference)
			.def("data", &ET::ScalarField<double>::data, py::return_value_policy::reference)
			.def("get_name", &ET::ScalarField<double>::getName, py::return_value_policy::reference)
			.def("get_N", &ET::ScalarField<double>::getN, py::return_value_policy::reference)
			.def("get_approximator", &ET::ScalarField<double>::getApproximator, py::return_value_policy::reference)
			.def("set_ugrid", &ET::ScalarField<double>::setUGrid)
			.def("set_field", &ET::ScalarField<double>::setField)
			.def("set_name", &ET::ScalarField<double>::setName)
			.def("set_approx_type", &ET::ScalarField<double>::setApproxType)
			//	Operator overloads
			.def("__getitem__", [](const ET::ScalarField<double> &self, int i)
			{
				if (i < 0 || i >= self.getN())
				{
					throw py::index_error("Index " + std::to_string(i) +
																" out of bounds for scalar field with "
																+ std::to_string(self.getN()) + " points!");
				}
				return self(i);
			}, py::is_operator())
			.def("__setitem__", [](ET::ScalarField<double> &self,
														 int i, const double& val)
			{
				if (i < 0 || i >= self.getN())
				{
					throw py::index_error("Index " + std::to_string(i) +
															" out of bounds for scalar field with "
															+ std::to_string(self.getN()) + " points!");
				}
				self(i) = val;
			}, py::is_operator())
			;

		py::class_<ET::Approximator<double>, std::shared_ptr<ET::Approximator<double>>>(m, "Approximator")
			.def(py::init<>())
			.def(py::init<int>())
			.def(py::init<std::string>())
			//	Getters and setters
			.def("get_approx_type", &ET::Approximator<double>::getApproxType, py::return_value_policy::reference)
			.def("get_approx_params", &ET::Approximator<double>::getApproxParams, py::return_value_policy::reference)
			.def("get_lsdriver", &ET::Approximator<double>::getLSDriver, py::return_value_policy::reference)
			.def("get_flag", &ET::Approximator<double>::getFlag, py::return_value_policy::reference)
			.def("get_info", &ET::Approximator<double>::getInfo, py::return_value_policy::reference)
			.def("set_approx_type", &ET::Approximator<double>::setApproxType)
			.def("set_approx_params", &ET::Approximator<double>::setApproxParams)
			.def("set_lsdriver", &ET::Approximator<double>::setLSDriver)
			.def("set_k", &ET::Approximator<double>::set_k)
			.def("set_n", &ET::Approximator<double>::set_n)
			.def("set_flag", &ET::Approximator<double>::setFlag)
			.def("set_info", &ET::Approximator<double>::setInfo)
			//	scalar field methods
			.def("scalar_gradient", &ET::Approximator<double>::scalarGradient)
			.def("scalar_gradient_mls", &ET::Approximator<double>::scalarGradientMLS)
			.def("construct_taylor_matrix", (ET::Matrix<double> (ET::Approximator<double>::*)(const std::shared_ptr<ET::UGrid<double>>,
	                                    const std::vector<uint64_t>,uint64_t,uint64_t)) &ET::Approximator<double>::constructTaylorMatrix)
			.def("construct_taylor_matrix", (ET::Matrix<double> (ET::Approximator<double>::*)(const std::shared_ptr<ET::UGrid<double>>,
	                                    const std::vector<uint64_t>,uint64_t,ET::Monomial&)) &ET::Approximator<double>::constructTaylorMatrix)
			//	print functionality
			.def("__repr__", [](const ET::Approximator<double> &app)
			{
					return app.summary();
			})
			;

		py::class_<ET::ApproxParams>(m, "ApproxParams")
			.def(py::init<>())
			.def_readwrite("k", &ET::ApproxParams::k)
			.def_readwrite("n", &ET::ApproxParams::n)
			;

		py::class_<ET::Monomial>(m, "Monomial")
			.def(py::init<>())
			.def(py::init<uint32_t>())
			.def(py::init<uint32_t,uint32_t>())
			.def("get_dim", &ET::Monomial::getDim, py::return_value_policy::reference)
			.def("get_deg", &ET::Monomial::getDeg, py::return_value_policy::reference)
			.def("get_mono", (std::vector<std::vector<uint32_t>>
				   (ET::Monomial::*)()) &ET::Monomial::getMono, py::return_value_policy::reference)
			.def("get_mono", (std::vector<std::vector<uint32_t>>
					 (ET::Monomial::*)(uint32_t)) &ET::Monomial::getMono, py::return_value_policy::reference)
			.def("get_mono_factors", &ET::Monomial::getMonoFactors, py::return_value_policy::reference)
			.def("get_multiset_coeff", &ET::Monomial::getMultisetCoefficient, py::return_value_policy::reference)
			.def("get_taylor_index", &ET::Monomial::getTaylorIndex, py::return_value_policy::reference)
			.def("set_dim", &ET::Monomial::setDim)
			.def("set_deg", &ET::Monomial::setDeg)
			.def("generate_monomial", (void (ET::Monomial::*)())
			     &ET::Monomial::generateMonomial)
			.def("generate_monomial", (void (ET::Monomial::*)(uint32_t))
			 		 &ET::Monomial::generateMonomial)
			.def("taylor_monomial_expansion", (std::vector<double>
				   (ET::Monomial::*)(const std::vector<double>&,const std::vector<double>&))
				   &ET::Monomial::taylorMonomialExpansion)
			.def("taylor_monomial_expansion", (std::vector<double>
		 		   (ET::Monomial::*)(const std::vector<double>&,const std::vector<double>&,uint32_t))
		 		   &ET::Monomial::taylorMonomialExpansion)
		  //	print functionality
			.def("__repr__", [](const ET::Monomial &mon)
			{
					return mon.summary();
			})
			;

		//	various functions
		m.def("taylor_polynomial", &ET::taylorPolynomial);
}
