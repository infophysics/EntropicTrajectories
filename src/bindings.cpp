//------------------------------------------------------------------------------
//  bindings.cpp
//  The Entropic Trajectories Framework
//  -----------------------------------
//  Copyright (C) [2020] by [N. Carrara]
//  [ncarrara@albany.edu]

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

#include "vector.h"
#include "matrix.h"
#include "grid.h"
#include "ugrid.h"
#include "log.h"
#include "kdtree.h"
#include "scalarfield.h"
#include "vectorfield.h"
#include "random.h"
#include "utilities.h"
#include "radial_basis_interpolation.h"
#include "local_taylor_interpolation.h"
#include "interpolator.h"
#include "interpolant.h"
#include "local_taylor_interpolant.h"
#include "radial_basis_interpolant.h"
#include "diffeq.h"
#include "dynamicalsystem.h"
#include "scalarfieldexample.h"

namespace py = pybind11;
using namespace ET;

PYBIND11_MODULE(etraj, m) {
  m.doc() = R"pbdoc(
          Pybind11 example plugin
          -----------------------
          .. currentmodule:: etraj.etraj
          .. autosummary::
             :toctree: _generate
      )pbdoc";
	//----------------------------------------------------------------------------
	//	Vector class
	//----------------------------------------------------------------------------
  py::class_<Vector<double>,
	           std::shared_ptr<Vector<double>>>
						 (m, "Vector", py::buffer_protocol())
    .def(py::init<>())
		.def(py::init<Vector<double>>())
    .def(py::init<size_t>(),
         py::arg("dim"))
		.def(py::init<std::string, size_t>(),
         py::arg("name"), py::arg("dim"))
    .def(py::init<std::vector<double>>(),
         py::arg("vec"))
		.def(py::init<std::string, std::vector<double>>(),
         py::arg("name"), py::arg("vec"))
		.def(py::init<size_t, const double&>(),
         py::arg("dim"), py::arg("init"))
		.def(py::init<std::string, size_t, const double&>(),
         py::arg("name"), py::arg("dim"), py::arg("init"))
		.def(py::init([](py::buffer const b)
		{
    	py::buffer_info info = b.request();
    	if (info.format != py::format_descriptor<float>::format()
					|| info.ndim != 1) {
      	throw std::runtime_error("Incompatible buffer format!");
			}
      auto v = new Vector<double>(info.shape[0]);
      memcpy(v->getVec().data(), info.ptr,
						 sizeof(double) * (size_t) (v->getDim()));
      return v;
    }))
    //  Getters and Setters
		.def("get_dim", &Vector<double>::getDim)
		.def("get_vec", &Vector<double>::getVec)
		.def("get_name", &Vector<double>::getName)
    .def("get_flag", &Vector<double>::getFlag)
    .def("get_info", &Vector<double>::getInfo)
		.def("set_dim", &Vector<double>::setDim)
		.def("set_vec", &Vector<double>::setVec)
		.def("set_name", &Vector<double>::setName)
    .def("set_flag", &Vector<double>::setFlag)
    .def("set_info", &Vector<double>::setInfo)
    .def_property("name", &Vector<double>::getName, &Vector<double>::setName)
    .def_property("dim", &Vector<double>::getDim, &Vector<double>::setDim)
    .def_property("flag", &Vector<double>::getFlag, &Vector<double>::setFlag)
    .def_property("info", &Vector<double>::getInfo, &Vector<double>::setInfo)
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
		.def("__mul__", [](const Vector<double>& v,
											 const Vector<double>& u)
		{
			return v * u;
		}, py::is_operator())
    .def("dot", &Vector<double>::dot)
    .def(py::self * double())
    .def(py::self *= double())
    .def(py::self / double())
    .def(py::self /= double())
		.def(-py::self)
		.def(double() + py::self)
		.def(double() - py::self)
		.def(double() * py::self)
		.def(double() / py::self)
		.def("__getitem__", [](const Vector<double> &self, int i)
		{
			if (i < 0 || i >= self.getDim())
			{
				throw py::index_error("Index " + std::to_string(i) +
															" out of bounds for vector with dimension "
															+ std::to_string(self.getDim()) + "!");
			}
			return self(i);
		}, py::is_operator())
		.def("__setitem__", [](Vector<double> &self,
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
    //    __len__
    .def("__len__", &Vector<double>::getDim)
    //    __repr__
    .def("__repr__", [](const Vector<double> &v) {
      return "<etraj.Vector<double> named '" + v.getName() + "'>";
    })
		//	  __str__
		.def("__str__", [](Vector<double> &v)
		{
				return v.summary();
		})
		;
		//--------------------------------------------------------------------------

		//--------------------------------------------------------------------------
		//	Level 1 BLAS methods
		//--------------------------------------------------------------------------
		m.def("dswap", &DSWAP);
		m.def("dscal", &DSCAL);
		m.def("dcopy", (void (*)(Vector<double>&, Vector<double>&)) &DCOPY);
		m.def("dcopy", (Vector<double> (*)(Vector<double>&)) &DCOPY);
		m.def("daxpy", &DAXPY);
		m.def("ddot",  &DDOT);
		m.def("dnrm2", &DNRM2);
		m.def("dasum", &DASUM);
		m.def("idamax",&IDAMAX);
		m.def("idamin",&IDAMIN);
	//--------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Matrix class
	//----------------------------------------------------------------------------
	py::class_<Matrix<double>,
	           std::shared_ptr<Matrix<double>>>
						(m, "Matrix", py::buffer_protocol())
    .def(py::init<>())
    .def(py::init<Matrix<double>>())
    .def(py::init<size_t>(),
         py::arg("m"))
    .def(py::init<std::string, size_t>(),
         py::arg("name"), py::arg("m"))
    .def(py::init<size_t, size_t>(),
         py::arg("m"), py::arg("n"))
    .def(py::init<std::string, size_t, size_t>(),
         py::arg("name"), py::arg("m"), py::arg("n"))
    .def(py::init<size_t, size_t, const double&>(),
         py::arg("m"), py::arg("n"), py::arg("init"))
    .def(py::init<std::string, size_t, size_t, const double&>(),
         py::arg("name"), py::arg("m"), py::arg("n"), py::arg("init"))
    .def(py::init<size_t, std::vector<double>>(),
         py::arg("m"), py::arg("flat"))
    .def(py::init<std::string, size_t, std::vector<double>>(),
         py::arg("name"), py::arg("m"), py::arg("flat"))
    .def(py::init<size_t, size_t, std::vector<double>>(),
         py::arg("m"), py::arg("n"), py::arg("flat"))
    .def(py::init<std::string, size_t, size_t, std::vector<double>>(),
         py::arg("name"), py::arg("m"), py::arg("n"), py::arg("flat"))
		.def(py::init<std::vector<std::vector<double>>>(),
         py::arg("array"))
		.def(py::init<std::string, std::vector<std::vector<double>>>(),
         py::arg("name"), py::arg("array"))
		.def(py::init([](py::buffer const b)
		{
    	py::buffer_info info = b.request();
      if (info.format != py::format_descriptor<float>::format()
					|| info.ndim != 2) {
      	throw std::runtime_error("Incompatible buffer format!");
			}
      auto v = new Matrix<double>(info.shape[0], info.shape[1]);
      memcpy(v->getArray().data(), info.ptr,
						 sizeof(double) * (size_t) (v->getNumRows() * v->getNumCols()));
      return v;
    }))
		//	getters
		.def("get_num_rows", &Matrix<double>::getNumRows,
		     py::return_value_policy::reference)
		.def("get_num_cols", &Matrix<double>::getNumCols,
		     py::return_value_policy::reference)
    .def("get_m", &Matrix<double>::getNumRows,
         py::return_value_policy::reference)
    .def("get_n", &Matrix<double>::getNumCols,
         py::return_value_policy::reference)
		.def("get_name", &Matrix<double>::getName,
		     py::return_value_policy::reference)
		.def("get_array", &Matrix<double>::getArray,
		     py::return_value_policy::reference)
		.def("get_row", &Matrix<double>::getRow,
		     py::return_value_policy::reference)
		.def("get_col", &Matrix<double>::getCol,
		     py::return_value_policy::reference)
		.def("get_info", &Matrix<double>::getInfo,
		     py::return_value_policy::reference)
		.def("get_flag", &Matrix<double>::getFlag,
		     py::return_value_policy::reference)
		.def("get_rank", &Matrix<double>::getRank,
		     py::return_value_policy::reference)
		.def("set_name", &Matrix<double>::setName)
		.def("set_row", &Matrix<double>::setRow)
		.def("set_col", &Matrix<double>::setCol)
		.def("set_array", (void (Matrix<double>::*)
				 (size_t,std::vector<double>)) &Matrix<double>::setArray)
		.def("set_array", (void (Matrix<double>::*)
				 (std::vector<std::vector<double>>)) &Matrix<double>::setArray)
		.def("set_singular_values", &Matrix<double>::setSingularValues)
		.def("set_info", &Matrix<double>::setInfo)
		.def("set_flag", &Matrix<double>::setFlag)
		.def("set_rank", &Matrix<double>::setRank)
    .def_property("name", &Matrix<double>::getName,
                  &Matrix<double>::setName)
    .def_property_readonly("m", &Matrix<double>::getNumRows)
    .def_property_readonly("n", &Matrix<double>::getNumCols)
    .def_property("flag", &Matrix<double>::getFlag,
                  &Matrix<double>::setFlag)
    .def_property("info", &Matrix<double>::getInfo,
                  &Matrix<double>::setInfo)
    .def("get_singular_values", &Matrix<double>::getSingularValues,
		     py::return_value_policy::reference)
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
    .def("__getitem__", [](const Matrix<double> &self,
													 std::tuple<size_t, size_t> ij)
		{
			size_t i, j;
			std::tie(i, j) = ij;
			if (i < 0 || i >= self.getNumRows()) {
				throw py::index_error("Index " + std::to_string(i) +
																" out of bounds for array with "
																+ std::to_string(self.getNumRows())
																+ " rows!");
			}
			if (j < 0 || j >= self.getNumCols()) {
				throw py::index_error("Index " + std::to_string(j) +
															  " out of bounds for array with "
															  + std::to_string(self.getNumCols())
																+ " columns!");
			}
			return self(i,j);
		}, py::is_operator())
		//	now for the setter of the same type.  To set an element to a Matrix
		//	at index i,j, write - m[i,j] = x.
		.def("__setitem__", [](Matrix<double> &self,
													 std::tuple<size_t, size_t> ij,
													 const double& val)
		{
			size_t i, j;
			std::tie(i, j) = ij;
			if (i < 0 || i >= self.getNumRows()) {
				throw py::index_error("Index " + std::to_string(i) +
																" out of bounds for array with "
																+ std::to_string(self.getNumRows())
																+ " rows!");
			}
			if (j < 0 || j >= self.getNumCols()) {
				throw py::index_error("Index " + std::to_string(j) +
															  " out of bounds for array with "
															  + std::to_string(self.getNumCols())
																+ " columns!");
			}
			self(i,j) = val;
		}, py::is_operator())
		//	this will allow the user to get entire rows.
		.def("__getitem__", [](Matrix<double> &self, int i)
		{
			if (i < 0 || i >= self.getNumRows()) {
				throw py::index_error("Index " + std::to_string(i) +
															" out of bounds for array with "
															+ std::to_string(self.getNumRows())
															+ " rows!");
			}
			return self.getRow(i);
		}, py::is_operator())
    //    __repr__
    .def("__repr__", [](const Matrix<double> &m) {
      return "<etraj.Matrix<double> named '" + m.getName() + "'>";
    })
		//	  __str__
		.def("__str__", [](Matrix<double> &m)
		{
				return m.summary();
		})
		//	Buffer definition allows one to cast a Matrix as a numpy array
		.def_buffer([](Matrix<double> &m) -> py::buffer_info
		{
      return py::buffer_info(
				m.data(),                               /* Pointer to buffer */
				sizeof(float),                          /* Size of one scalar */
				py::format_descriptor<float>::format(), /* format descriptor */
				2,                                      /* Number of dimensions */
				{ m.getNumRows(), m.getNumCols() },     /* Buffer dimensions */
				{ sizeof(float) * m.getNumCols(),       /* Strides for each index */
					sizeof(float) });
		})
		//--------------------------------------------------------------------------
		//	Various methods
		//--------------------------------------------------------------------------
		.def("print", &Matrix<double>::print)
		.def("transpose", (Matrix<double> (Matrix<double>::*)())
				 &Matrix<double>::transpose)
		.def("transpose", (void (Matrix<double>::*)(bool))
				 &Matrix<double>::transpose_inplace)
		.def("trace", &Matrix<double>::trace)
		;
	//--------------------------------------------------------------------------

	//--------------------------------------------------------------------------
	//	Instantiations of matrices
	//--------------------------------------------------------------------------
	m.def("identity", (Matrix<double> (*)(size_t)) &identity_d);
	m.def("permutation_matrix", (Matrix<double> (*)(const size_t&,
	      const std::vector<size_t>)) &permutationMatrix_d);
	//--------------------------------------------------------------------------
	//	Level 2 BLAS functions
	//--------------------------------------------------------------------------
	m.def("dgemv", (Vector<double> (*)(double&,Matrix<double>&,
				Vector<double>&)) &DGEMV);
	m.def("dgemv", (void (*)(double&,Matrix<double>&,
				Vector<double>&, double&, Vector<double>&)) &DGEMV);
	m.def("dger", (Matrix<double> (*)(double&,Vector<double>&,
				Vector<double>&)) &DGER);
	m.def("dger", (void (*)(double&,Vector<double>&,
				Vector<double>&, Matrix<double>&)) &DGER);
	//--------------------------------------------------------------------------

	//--------------------------------------------------------------------------
	//  Level 3 BLAS functions
	//--------------------------------------------------------------------------
	m.def("dgemm", (Matrix<double> (*)(const double&,
		    const Matrix<double>&, const Matrix<double>&))
				&DGEMM);
	m.def("dgemm", (void (*)(const double&, const Matrix<double>&,
		    const Matrix<double>&, const double&,
				Matrix<double>&)) &DGEMM);
	//	DGEMM with no scalar "alpha"
	m.def("dgemm", (Matrix<double> (*)(const Matrix<double>&,
		    const Matrix<double>&)) &DGEMM, py::keep_alive<1, 2>());
	//--------------------------------------------------------------------------

	//--------------------------------------------------------------------------
	//	LAPACK functions
	//--------------------------------------------------------------------------
	m.def("dgels", (Matrix<double> (*)(const Matrix<double>&,
	      const Matrix<double>&)) &DGELS, py::keep_alive<0, 1>());
	m.def("dgels", (Vector<double> (*)(const Matrix<double>&,
	      const Vector<double>&)) &DGELS, py::keep_alive<0, 1>());
	m.def("dgelsy", (Matrix<double> (*)(const Matrix<double>&,
	      const Matrix<double>&)) &DGELSY);
	m.def("dgelsy", (Vector<double> (*)(const Matrix<double>&,
	      const Vector<double>&)) &DGELSY);
	m.def("dgelsd", (Matrix<double> (*)(const Matrix<double>&,
	      const Matrix<double>&)) &DGELSD);
	m.def("dgelsd", (Vector<double> (*)(const Matrix<double>&,
	      const Vector<double>&)) &DGELSD);
	m.def("dgelss", (Matrix<double> (*)(const Matrix<double>&,
	      const Matrix<double>&)) &DGELSS);
	m.def("dgelss", (Vector<double> (*)(const Matrix<double>&,
	      const Vector<double>&)) &DGELSS);
	m.def("dgetrf", &DGETRF);
	m.def("dgetrf_l_u", &DGETRF_L_U);
	m.def("dgetrf_lu", &DGETRF_LU);
	m.def("dgetrf_plu", &DGETRF_PLU);
	m.def("dgeqrf", &DGEQRF);
	m.def("dorgqr", &DORGQR);
	m.def("dgeqrf_qr", &DGEQRF_QR);
	m.def("dgesvd", &DGESVD);
	m.def("dgesvd_svd", &DGESVD_SVD);
	m.def("dgesdd", &DGESDD);
	m.def("dgesdd_svd", &DGESDD_SVD);
	m.def("dgetri", &DGETRI);
	//--------------------------------------------------------------------------

  //--------------------------------------------------------------------------
  //  GridType enum
  //--------------------------------------------------------------------------
  py::enum_<GridType>(m, "GridType")
		.value("default", GridType::DEFAULT)
    .value("structured", GridType::STRUCTURED)
    .value("unstructured", GridType::UNSTRUCTURED)
    .value("meshless", GridType::MESHLESS)
		;
  //--------------------------------------------------------------------------

  //----------------------------------------------------------------------------
	//	Grid Base class
	//----------------------------------------------------------------------------
  py::class_<Grid<double>, std::shared_ptr<Grid<double>>>(m, "Grid")
    .def(py::init<>())
    .def(py::init<std::shared_ptr<Log>>(),
         py::arg("log"), R"pbdoc(
           Constructor with a shared Log instance.
         )pbdoc")
    .def(py::init<std::string>(),
         py::arg("name"))
    .def(py::init<std::string,std::shared_ptr<Log>>(),
         py::arg("name"), py::arg("log"))
    .def(py::init<size_t>(),
         py::arg("dim"))
    .def(py::init<std::string,size_t>(),
         py::arg("name"), py::arg("dim"))
    .def(py::init<size_t,std::shared_ptr<Log>>(),
         py::arg("dim"), py::arg("log"))
    .def(py::init<std::string,size_t,std::shared_ptr<Log>>(),
         py::arg("name"), py::arg("dim"), py::arg("log"))
    .def(py::init<size_t,size_t>(),
         py::arg("dim"), py::arg("N"))
    .def(py::init<size_t,size_t,std::shared_ptr<Log>>(),
         py::arg("dim"), py::arg("N"), py::arg("log"))
    .def(py::init<std::string,size_t,size_t>(),
         py::arg("name"), py::arg("dim"), py::arg("N"))
    .def(py::init<std::string,size_t,size_t,std::shared_ptr<Log>>(),
         py::arg("name"), py::arg("dim"), py::arg("N"), py::arg("log"))
    .def(py::init<std::vector<std::vector<double>>,bool>(),
         py::arg("grid"), py::arg("move_grid")=false, py::keep_alive<0, 1>())
    .def(py::init<std::vector<std::vector<double>>,std::shared_ptr<Log>,bool>(),
         py::arg("grid"), py::arg("log"), py::arg("move_grid")=false,
         py::keep_alive<0, 1>())
    .def(py::init<std::string,std::vector<std::vector<double>>,bool>(),
         py::arg("name"), py::arg("grid"), py::arg("move_grid")=false,
         py::keep_alive<0, 2>())
    .def(py::init<std::string,std::vector<std::vector<double>>,
                  std::shared_ptr<Log>,bool>(),
        py::arg("name"), py::arg("grid"), py::arg("log"), py::arg("move_grid")=false,
        py::keep_alive<0, 2>())
    //  Getters and Setters
    .def("get_name", &Grid<double>::getName)
    .def("get_dim", &Grid<double>::getDim)
    .def("get_N", &Grid<double>::getN)
    .def("get_grid", &Grid<double>::getGrid,
         py::return_value_policy::reference)
    .def("get_coords", &Grid<double>::getCoords)
    .def("get_log", &Grid<double>::getLog)
    .def("get_type", &Grid<double>::getType)
    .def("get_point", &Grid<double>::getPoint,
         py::return_value_policy::reference)
    .def("set_name", &Grid<double>::setName)
    .def("set_dim", &Grid<double>::setDim)
    .def("set_N", &Grid<double>::setN)
    .def("set_grid", &Grid<double>::setGrid)
    .def("set_coords", &Grid<double>::setCoords)
    .def("set_log", &Grid<double>::setLog)
    .def("move_grid", &Grid<double>::moveGrid)
    .def_property("name", &Grid<double>::getName,
                  &Grid<double>::setName)
    .def_property("dim", &Grid<double>::getDim,
                  &Grid<double>::setDim)
    .def_property("N", &Grid<double>::getN,
                  &Grid<double>::setN)
    .def_property("grid", &Grid<double>::getGrid,
                  &Grid<double>::setGrid)
    .def_property("coords", &Grid<double>::getCoords,
                  &Grid<double>::setCoords)
    .def_property_readonly("type", &Grid<double>::getType)
    .def_property("log", &Grid<double>::getLog, &Grid<double>::setLog)
    //  Operator overloads
    .def(py::self == py::self)
		.def(py::self != py::self)
    .def(py::self + py::self)
    .def(py::self += py::self)
    .def(py::self - py::self)
    .def(py::self -= py::self)
    .def(-py::self)
    //	Matrix array access.  This first method allows the user
		//	to write x = m[i,j] to get elements.
    .def("__getitem__", [](const Grid<double> &self,
													 std::tuple<size_t, size_t> ij)
		{
			size_t i, j;
			std::tie(i, j) = ij;
			if (i < 0 || i >= self.getN()) {
				throw py::index_error("Index " + std::to_string(i) +
																" out of bounds for array with "
																+ std::to_string(self.getN())
																+ " rows!");
			}
			if (j < 0 || j >= self.getDim()) {
				throw py::index_error("Index " + std::to_string(j) +
															  " out of bounds for array with "
															  + std::to_string(self.getDim())
																+ " columns!");
			}
			return self(i,j);
		}, py::is_operator())
		//	now for the setter of the same type.  To set an element to a Grid
		//	at index i,j, write - m[i,j] = x.
		.def("__setitem__", [](Grid<double> &self,
													 std::tuple<size_t, size_t> ij,
													 const double& val)
		{
			size_t i, j;
			std::tie(i, j) = ij;
			if (i < 0 || i >= self.getN()) {
				throw py::index_error("Index " + std::to_string(i) +
																" out of bounds for array with "
																+ std::to_string(self.getN())
																+ " rows!");
			}
			if (j < 0 || j >= self.getDim()) {
				throw py::index_error("Index " + std::to_string(j) +
															  " out of bounds for array with "
															  + std::to_string(self.getDim())
																+ " columns!");
			}
			self(i,j) = val;
		}, py::is_operator())
		//	this will allow the user to get a point.
		.def("__getitem__", [](Grid<double> &self, int i)
		{
			if (i < 0 || i >= self.getN()) {
				throw py::index_error("Index " + std::to_string(i) +
															" out of bounds for array with "
															+ std::to_string(self.getN())
															+ " rows!");
			}
			return self(i);
		}, py::is_operator())
    //    __len__
    .def("__len__", &Grid<double>::getN)
    //    __repr__
    .def("__repr__", [](const Grid<double> &g) {
      return "<etraj.Grid<double> named '" + g.getName() + "'>";
    })
    //    __str__
    .def("__str__", [](const Grid<double> &g) {
      std::stringstream s;
			s << &g;
      std::string sum;
      sum += "++++++++++++++++++++++++++++++++++++++++++++++++++++";
			sum += "\n<etraj.Grid<double> ref at " + s.str() + ">";
      sum += "\n---------------------------------------------------";
      sum += "\n<ET::Grid<double> object at " + address_to_string(g) + ">";
      sum += "\n---------------------------------------------------";
      sum += "\n   name: '" + g.getName() + "'";
      sum += "\n    dim: " + std::to_string(g.getDim());
      sum += "\n      N: " + std::to_string(g.getN());
      sum += "\n---------------------------------------------------";
      sum += "\n Logger at: " + address_to_string(*g.getLog()) + ",";
      sum += "\n    ref at: " + address_to_string(g.getLog());
      sum += "\n++++++++++++++++++++++++++++++++++++++++++++++++++++";
      return sum;
    })
    //  Special functions
    .def("proj", (std::vector<double> (Grid<double>::*)(const size_t))
         &Grid<double>::proj)
    .def("proj", (std::vector<std::vector<double>> (Grid<double>::*)
         (const std::vector<size_t>)) &Grid<double>::proj)
    .def("log_output", [](Grid<double>& self)
		{
			std::string out = self.getLog()->getOutput();
			py::print(out);
		})
		.def("log_output", [](Grid<double>& self, size_t i)
		{
			std::string out = self.getLog()->getOutput(i);
			py::print(out);
		})
    ;
  //----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	UGrid class
	//----------------------------------------------------------------------------
	py::class_<UGrid<double>, Grid<double>,
             std::shared_ptr<UGrid<double>>>(m, "UGrid")
    .def(py::init<>())
    .def(py::init<std::shared_ptr<Log>>(),
        py::arg("log"), R"pbdoc(
          Constructor with a shared Log instance.
        )pbdoc")
    .def(py::init<std::string>(),
        py::arg("name"))
    .def(py::init<std::string,std::shared_ptr<Log>>(),
        py::arg("name"), py::arg("log"))
    .def(py::init<size_t>(),
        py::arg("dim"))
    .def(py::init<std::string,size_t>(),
        py::arg("name"), py::arg("dim"))
    .def(py::init<size_t,std::shared_ptr<Log>>(),
        py::arg("dim"), py::arg("log"))
    .def(py::init<std::string,size_t,std::shared_ptr<Log>>(),
        py::arg("name"), py::arg("dim"), py::arg("log"))
    .def(py::init<size_t,size_t>(),
        py::arg("dim"), py::arg("N"))
    .def(py::init<size_t,size_t,std::shared_ptr<Log>>(),
        py::arg("dim"), py::arg("N"), py::arg("log"))
    .def(py::init<std::string,size_t,size_t>(),
        py::arg("name"), py::arg("dim"), py::arg("N"))
    .def(py::init<std::string,size_t,size_t,std::shared_ptr<Log>>(),
        py::arg("name"), py::arg("dim"), py::arg("N"), py::arg("log"))
    .def(py::init<std::vector<std::vector<double>>,bool>(),
         py::arg("grid"), py::arg("move_grid")=false, py::keep_alive<0, 1>())
    .def(py::init<std::vector<std::vector<double>>,std::shared_ptr<Log>,bool>(),
         py::arg("grid"), py::arg("log"), py::arg("move_grid")=false,
         py::keep_alive<0, 1>())
    .def(py::init<std::string,std::vector<std::vector<double>>,bool>(),
         py::arg("name"), py::arg("grid"), py::arg("move_grid")=false,
         py::keep_alive<0, 2>())
    .def(py::init<std::string,std::vector<std::vector<double>>,
                  std::shared_ptr<Log>,bool>(),
        py::arg("name"), py::arg("grid"), py::arg("log"), py::arg("move_grid")=false,
        py::keep_alive<0, 2>())
    //  Getters and Setters
    .def("get_type", &UGrid<double>::getType)
    .def("get_kdtree", &UGrid<double>::getKDTree)
    .def_property_readonly("kdtree", &UGrid<double>::getKDTree)
    .def_property_readonly("type", &UGrid<double>::getType)
    //    __str__
    .def("__str__", [](const UGrid<double> &g) {
      std::stringstream s;
			s << &g;
      std::string sum;
      sum += "++++++++++++++++++++++++++++++++++++++++++++++++++++";
			sum += "\n<etraj.UGrid<double> ref at " + s.str() + ">";
      sum += "\n---------------------------------------------------";
      sum += "\n<ET::UGrid<double> object at " + address_to_string(g) + ">";
      sum += "\n---------------------------------------------------";
      sum += "\n   name: '" + g.getName() + "'";
      sum += "\n    dim: " + std::to_string(g.getDim());
      sum += "\n      N: " + std::to_string(g.getN());
      sum += "\n---------------------------------------------------";
      sum += "\n KDTree at: " + address_to_string(g.getKDTree()) + ",";
      sum += "\n---------------------------------------------------";
      sum += "\n Logger at: " + address_to_string(*g.getLog()) + ",";
      sum += "\n    ref at: " + address_to_string(g.getLog());
      sum += "\n++++++++++++++++++++++++++++++++++++++++++++++++++++";
      return sum;
    })
    //  KDTree overloads
    .def("get_current_neighbor_indices", (std::vector<std::vector<size_t>>
         (UGrid<double>::*)()) &UGrid<double>::getCurrentNeighborIndices)
    .def("get_current_neighbor_distances", (std::vector<std::vector<double>>
         (UGrid<double>::*)()) &UGrid<double>::getCurrentNeighborDistances)
    .def("get_current_neighbor_indices",
         (std::vector<size_t> (UGrid<double>::*)(const size_t))
         &UGrid<double>::getCurrentNeighborIndices,
         py::arg("index"))
    .def("get_current_neighbor_distances",
         (std::vector<double> (UGrid<double>::*)(const size_t))
         &UGrid<double>::getCurrentNeighborDistances,
         py::arg("index"))
    .def("get_current_global_k", &UGrid<double>::getCurrentGlobalK)
    .def("get_current_global_radius", &UGrid<double>::getCurrentGlobalRadius)
    .def("set_current_global_k", &UGrid<double>::setCurrentGlobalK)
    .def("set_current_global_radius", &UGrid<double>::setCurrentGlobalRadius)
    .def("setup_tree", &UGrid<double>::setupTree)
    //  Overloads of query_neighbors
    .def("query_neighbors", (void (UGrid<double>::*)(size_t))
         &UGrid<double>::queryNeighbors,
         py::arg("k"))
    .def("query_neighbors", (void (UGrid<double>::*)(double))
         &UGrid<double>::queryNeighbors,
         py::arg("radius"))
    .def("query_neighbors", (std::vector<size_t> (UGrid<double>::*)
         (const std::vector<double>&,size_t)) &UGrid<double>::queryNeighbors,
         py::arg("point"), py::arg("k"))
    .def("query_distances", (std::vector<double> (UGrid<double>::*)
         (const std::vector<double>&,size_t)) &UGrid<double>::queryDistances,
         py::arg("point"), py::arg("k"))
    .def("query", (std::tuple<std::vector<size_t>,std::vector<double>>
         (UGrid<double>::*)(const std::vector<double>&,size_t))
         &UGrid<double>::query,
         py::arg("point"), py::arg("k"))
    .def("query_neighbors", (std::vector<size_t> (UGrid<double>::*)
         (const std::vector<double>&,double)) &UGrid<double>::queryNeighbors,
         py::arg("point"), py::arg("radius"))
    .def("query_distances", (std::vector<double> (UGrid<double>::*)
         (const std::vector<double>&,double)) &UGrid<double>::queryDistances,
         py::arg("point"), py::arg("radius"))
    .def("query", (std::tuple<std::vector<size_t>,std::vector<double>>
         (UGrid<double>::*)(const std::vector<double>&,double))
         &UGrid<double>::query,
         py::arg("point"), py::arg("radius"))
    .def("query_neighbors", (std::vector<std::vector<size_t>>
         (UGrid<double>::*)(const std::vector<std::vector<double>>&,size_t))
         &UGrid<double>::queryNeighbors,
         py::arg("points"), py::arg("k"))
    .def("query_distances", (std::vector<std::vector<double>>
         (UGrid<double>::*)(const std::vector<std::vector<double>>&,size_t))
         &UGrid<double>::queryDistances,
         py::arg("points"), py::arg("k"))
    .def("query", (std::tuple<std::vector<std::vector<size_t>>,
                              std::vector<std::vector<double>>>
        (UGrid<double>::*)(const std::vector<std::vector<double>>&,size_t))
        &UGrid<double>::query,
        py::arg("points"), py::arg("k"))
    .def("query_neighbors", (std::vector<std::vector<size_t>>
         (UGrid<double>::*)(const std::vector<std::vector<double>>&,double))
         &UGrid<double>::queryNeighbors,
         py::arg("points"), py::arg("radius"))
    .def("query_distances", (std::vector<std::vector<double>>
         (UGrid<double>::*)(const std::vector<std::vector<double>>&,double))
         &UGrid<double>::queryDistances,
         py::arg("points"), py::arg("radius"))
    .def("query", (std::tuple<std::vector<std::vector<size_t>>,
                             std::vector<std::vector<double>>>
       (UGrid<double>::*)(const std::vector<std::vector<double>>&,double))
       &UGrid<double>::query,
       py::arg("points"), py::arg("radius"))
		;
	//----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Enums from KDTree
  //----------------------------------------------------------------------------
  py::enum_<KDTreeBackend>(m, "KDTreeBackend")
    .value("nanoflann", KDTreeBackend::NANOFLANN)
    .value("flann", KDTreeBackend::FLANN)
    ;
  py::enum_<KDTreeSearchFlags>(m, "KDTreeSearchFlags")
    .value("outdated", KDTreeSearchFlags::OUTDATED)
    .value("built", KDTreeSearchFlags::BUILT)
    .value("new_k", KDTreeSearchFlags::NEW_K)
    .value("new_radius", KDTreeSearchFlags::NEW_RADIUS)
    .value("new_points", KDTreeSearchFlags::NEW_POINTS)
    .value("current_k", KDTreeSearchFlags::CURRENT_K)
    .value("current_radius", KDTreeSearchFlags::CURRENT_RADIUS)
    ;
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //	KDTree class
  //----------------------------------------------------------------------------
  py::class_<KDTree<double>,
             std::shared_ptr<KDTree<double>>>(m, "KDTree")
    .def(py::init<>())
    .def(py::init([](const std::vector<std::vector<double>> points) {
      return std::shared_ptr<KDTree<double>>(new KDTree<double>
      (std::make_shared<std::vector<std::vector<double>>>(points)));
    }))
    .def(py::init([](const std::string name,
                     const std::vector<std::vector<double>> points) {
      return std::shared_ptr<KDTree<double>>(new KDTree<double>
      (name,std::make_shared<std::vector<std::vector<double>>>(points)));
    }))
    .def("get_backend", &KDTree<double>::getBackend)
    .def("get_name", &KDTree<double>::getName)
    .def("get_dim", &KDTree<double>::getDim)
    .def("get_N", &KDTree<double>::getN)
    .def("get_points", [](const KDTree<double>& self) {
      return *self.getPoints();
    }, py::return_value_policy::reference)
    .def("get_current_neighbor_indices", (std::vector<std::vector<size_t>>
         (KDTree<double>::*)()) &KDTree<double>::getCurrentNeighborIndices)
    .def("get_current_neighbor_distances", (std::vector<std::vector<double>>
         (KDTree<double>::*)()) &KDTree<double>::getCurrentNeighborDistances)
    .def("get_current_neighbor_indices",
         (std::vector<size_t> (KDTree<double>::*)(const size_t))
         &KDTree<double>::getCurrentNeighborIndices,
         py::arg("index"))
    .def("get_current_neighbor_distances",
         (std::vector<double> (KDTree<double>::*)(const size_t))
         &KDTree<double>::getCurrentNeighborDistances,
         py::arg("index"))
    .def("get_current_global_k", &KDTree<double>::getCurrentGlobalK)
    .def("get_current_global_radius", &KDTree<double>::getCurrentGlobalRadius)
    .def("get_log", &KDTree<double>::getLog)
    .def("get_searchflag", &KDTree<double>::getSearchFlag)
    .def("set_backend", &KDTree<double>::setBackend)
    .def("set_name", &KDTree<double>::setName)
    .def("set_dim", &KDTree<double>::setDim)
    .def("set_N", &KDTree<double>::setN)
    .def("set_points", [](KDTree<double>& self,
                          const std::vector<std::vector<double>> points) {
      self.setPoints(std::make_shared<std::vector<std::vector<double>>>(points));
    })
    .def("set_current_global_k", &KDTree<double>::setCurrentGlobalK)
    .def("set_current_global_radius", &KDTree<double>::setCurrentGlobalRadius)
    .def("set_log", &KDTree<double>::setLog)
    .def_property("backend", &KDTree<double>::getBackend,
                  &KDTree<double>::setBackend)
    .def_property("name", &KDTree<double>::getName, &KDTree<double>::setName)
    .def_property("dim", &KDTree<double>::getDim, &KDTree<double>::setDim)
    .def_property("N", &KDTree<double>::getN, &KDTree<double>::setN)
    .def_property("log", &KDTree<double>::getLog, &KDTree<double>::setLog)
    //  KDTree methods
    .def("setup_tree", &KDTree<double>::setupTree)
    //  Overloads of query_neighbors
    .def("query_neighbors", (void (KDTree<double>::*)(size_t))
         &KDTree<double>::queryNeighbors,
         py::arg("k"))
    .def("query_neighbors", (void (KDTree<double>::*)(double))
         &KDTree<double>::queryNeighbors,
         py::arg("radius"))
    .def("query_neighbors", (std::vector<size_t> (KDTree<double>::*)
         (const std::vector<double>&,size_t)) &KDTree<double>::queryNeighbors,
         py::arg("point"), py::arg("k"))
    .def("query_distances", (std::vector<double> (KDTree<double>::*)
         (const std::vector<double>&,size_t)) &KDTree<double>::queryDistances,
         py::arg("point"), py::arg("k"))
    .def("query", (std::tuple<std::vector<size_t>,std::vector<double>>
         (KDTree<double>::*)(const std::vector<double>&,size_t))
         &KDTree<double>::query,
         py::arg("point"), py::arg("k"))
    .def("query_neighbors", (std::vector<size_t> (KDTree<double>::*)
         (const std::vector<double>&,double)) &KDTree<double>::queryNeighbors,
         py::arg("point"), py::arg("radius"))
    .def("query_distances", (std::vector<double> (KDTree<double>::*)
         (const std::vector<double>&,double)) &KDTree<double>::queryDistances,
         py::arg("point"), py::arg("radius"))
    .def("query", (std::tuple<std::vector<size_t>,std::vector<double>>
         (KDTree<double>::*)(const std::vector<double>&,double))
         &KDTree<double>::query,
         py::arg("point"), py::arg("radius"))
    .def("query_neighbors", (std::vector<std::vector<size_t>>
         (KDTree<double>::*)(const std::vector<std::vector<double>>&,size_t))
         &KDTree<double>::queryNeighbors,
         py::arg("points"), py::arg("k"))
    .def("query_distances", (std::vector<std::vector<double>>
         (KDTree<double>::*)(const std::vector<std::vector<double>>&,size_t))
         &KDTree<double>::queryDistances,
         py::arg("points"), py::arg("k"))
    .def("query", (std::tuple<std::vector<std::vector<size_t>>,
                              std::vector<std::vector<double>>>
        (KDTree<double>::*)(const std::vector<std::vector<double>>&,size_t))
        &KDTree<double>::query,
        py::arg("points"), py::arg("k"))
    .def("query_neighbors", (std::vector<std::vector<size_t>>
         (KDTree<double>::*)(const std::vector<std::vector<double>>&,double))
         &KDTree<double>::queryNeighbors,
         py::arg("points"), py::arg("radius"))
    .def("query_distances", (std::vector<std::vector<double>>
         (KDTree<double>::*)(const std::vector<std::vector<double>>&,double))
         &KDTree<double>::queryDistances,
         py::arg("points"), py::arg("radius"))
    .def("query", (std::tuple<std::vector<std::vector<size_t>>,
                             std::vector<std::vector<double>>>
       (KDTree<double>::*)(const std::vector<std::vector<double>>&,double))
       &KDTree<double>::query,
       py::arg("points"), py::arg("radius"))
    .def("log_output", [](KDTree<double>& self)
    {
    	std::string out = self.getLog()->getOutput();
    	py::print(out);
    })
    .def("log_output", [](KDTree<double>& self, size_t i)
    {
    	std::string out = self.getLog()->getOutput(i);
    	py::print(out);
    })
    //    __str__
    .def("__str__", [](const KDTree<double> &kdt) {
     std::stringstream s;
     s << &kdt;
     std::string sum;
     sum += "++++++++++++++++++++++++++++++++++++++++++++++++++++";
     sum += "\n<etraj.KDTree<double> ref at " + s.str() + ">";
     sum += "\n---------------------------------------------------";
     sum += "\n<ET::KDTree<double> object at " + address_to_string(kdt) + ">";
     sum += "\n---------------------------------------------------";
     sum += "\n   name: '" + kdt.getName() + "'";
     sum += "\n    dim: " + std::to_string(kdt.getDim());
     sum += "\n      N: " + std::to_string(kdt.getN());
     sum += "\nBackend: " + KDTreeBackendNameMap[kdt.getBackend()];
     sum += "\n---------------------------------------------------";
     sum += "\n Logger at: " + address_to_string(*kdt.getLog()) + ",";
     sum += "\n    ref at: " + address_to_string(kdt.getLog());
     sum += "\n++++++++++++++++++++++++++++++++++++++++++++++++++++";
     return sum;
    })
    ;
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Field Type enum
  //----------------------------------------------------------------------------
  py::enum_<FieldType>(m, "FieldType")
    .value("default", FieldType::DEFAULT)
    .value("scalar", FieldType::SCALAR)
    .value("vector", FieldType::VECTOR)
    .value("frame", FieldType::FRAME)
    ;
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Field Base Class
  //----------------------------------------------------------------------------
  py::class_<Field<double>, std::shared_ptr<Field<double>>>(m, "Field")
    .def(py::init<>())
    .def(py::init<std::shared_ptr<Log>>(),
         py::arg("log"))
    .def(py::init<std::shared_ptr<Grid<double>>>(),
         py::arg("Grid"))
    .def(py::init<std::shared_ptr<Interpolator<double>>>(),
         py::arg("Interpolator"))
    .def(py::init<std::shared_ptr<Grid<double>>,std::shared_ptr<Log>>(),
         py::arg("Grid"), py::arg("log"))
    .def(py::init<std::string,std::shared_ptr<Log>>(),
         py::arg("name"), py::arg("log"))
    .def(py::init<std::string,std::shared_ptr<Grid<double>>>(),
         py::arg("name"), py::arg("Grid"))
    .def(py::init<std::string,std::shared_ptr<Interpolator<double>>>(),
         py::arg("name"), py::arg("Interpolator"))
    .def(py::init<std::string,std::shared_ptr<Grid<double>>,std::shared_ptr<Log>>(),
         py::arg("name"), py::arg("Grid"), py::arg("log"))
    //  Getters and Setters
    .def("get_name", &Field<double>::getName)
    .def("get_dim", &Field<double>::getDim)
    .def("get_N", &Field<double>::getN)
    .def("get_log", &Field<double>::getLog)
    .def("get_Grid", &Field<double>::getGrid)
    .def("get_Interpolator", &Field<double>::getInterpolator)
    .def("get_DiffEQ", &Field<double>::getDiffEQ)
    .def("get_Integrator", &Field<double>::getIntegrator)
    .def("get_flag", &Field<double>::getFlag)
    .def("get_info", &Field<double>::getInfo)
    .def("set_name", &Field<double>::setName)
    .def("set_dim", &Field<double>::setDim)
    .def("set_N", &Field<double>::setN)
    .def("set_log", &Field<double>::setLog)
    .def("set_Grid", &Field<double>::setGrid)
    .def("set_Interpolator", &Field<double>::setInterpolator)
    .def("set_DiffEQ", &Field<double>::setDiffEQ)
    .def("set_Integrator", &Field<double>::setIntegrator)
    .def("set_flag", &Field<double>::setFlag)
    .def("set_info", &Field<double>::setInfo)
    .def("get_type", &Field<double>::getType)
    .def_property("name", &Field<double>::getName, &Field<double>::setName)
    .def_property("dim", &Field<double>::getDim, &Field<double>::setDim)
    .def_property("N", &Field<double>::getN, &Field<double>::setN)
    .def_property("log", &Field<double>::getLog, &Field<double>::setLog)
    .def_property("Grid", &Field<double>::getGrid, &Field<double>::setGrid)
    .def_property("Interpolator", &Field<double>::getInterpolator,
                  &Field<double>::setInterpolator)
    .def_property("DiffEQ", &Field<double>::getDiffEQ,
                  &Field<double>::setDiffEQ)
    .def_property("Integrator", &Field<double>::getIntegrator,
                  &Field<double>::setIntegrator)
    .def_property("flag", &Field<double>::getFlag, &Field<double>::setFlag)
    .def_property("info", &Field<double>::getInfo, &Field<double>::setInfo)
    .def_property_readonly("type", &Field<double>::getType)
    .def("log_output", [](Field<double>& self)
    {
    	std::string out = self.getLog()->getOutput();
    	py::print(out);
    })
    .def("log_output", [](Field<double>& self, size_t i)
    {
    	std::string out = self.getLog()->getOutput(i);
    	py::print(out);
    })
    //    __len__
    .def("__len__", &Field<double>::getN)
    //    __repr__
    .def("__repr__", [](const Field<double> &f) {
      return "<etraj.Field<double> named '" + f.getName() + "'>";
    })
    //    __str__
    .def("__str__", [](const Field<double> &f) {
      std::stringstream s;
			s << &f;
      std::string sum;
      sum += "++++++++++++++++++++++++++++++++++++++++++++++++++++";
			sum += "\n<etraj.Field<double> ref at " + s.str() + ">";
      sum += "\n---------------------------------------------------";
      sum += "\n<ET::Field<double> object at " + address_to_string(f) + ">";
      sum += "\n---------------------------------------------------";
      sum += "\n   name: '" + f.getName() + "'";
      sum += "\n    dim: " + std::to_string(f.getDim());
      sum += "\n      N: " + std::to_string(f.getN());
      sum += "\n---------------------------------------------------";
      sum += "\n        Grid at: " + address_to_string(*f.getGrid()) + ",";
      sum += "\n         ref at: " + address_to_string(f.getGrid());
      sum += "\n---------------------------------------------------";
      sum += "\nInterpolator at: " + address_to_string(*f.getInterpolator()) + ",";
      sum += "\n         ref at: " + address_to_string(f.getInterpolator());
      sum += "\n---------------------------------------------------";
      sum += "\n      DiffEQ at: " + address_to_string(*f.getDiffEQ()) + ",";
      sum += "\n         ref at: " + address_to_string(f.getDiffEQ());
      sum += "\n---------------------------------------------------";
      sum += "\n  Integrator at: " + address_to_string(*f.getIntegrator()) + ",";
      sum += "\n         ref at: " + address_to_string(f.getIntegrator());
      sum += "\n---------------------------------------------------";
      sum += "\n      Logger at: " + address_to_string(*f.getLog()) + ",";
      sum += "\n         ref at: " + address_to_string(f.getLog());
      sum += "\n++++++++++++++++++++++++++++++++++++++++++++++++++++";
      return sum;
    })
    ;
  //----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	ScalarField class
	//----------------------------------------------------------------------------
	py::class_<ScalarField<double>, Field<double>,
	           std::shared_ptr<ScalarField<double>>>(m, "ScalarField")
    .def(py::init<>())
    .def(py::init<std::string>(),
         py::arg("name"))
    .def(py::init<std::shared_ptr<Log>>(),
        py::arg("log"))
    .def(py::init<std::shared_ptr<Grid<double>>>(),
        py::arg("Grid"))
    .def(py::init<std::shared_ptr<Interpolator<double>>>(),
        py::arg("Interpolator"))
    .def(py::init<std::shared_ptr<Grid<double>>,std::shared_ptr<Log>>(),
        py::arg("Grid"), py::arg("log"))
    .def(py::init<std::vector<double>>(),
         py::arg("field"))
    .def(py::init<std::vector<double>,std::shared_ptr<Log>>(),
         py::arg("field"), py::arg("log"))
    .def(py::init<std::vector<double>,std::shared_ptr<Grid<double>>>(),
         py::arg("field"), py::arg("Grid"))
    .def(py::init<std::vector<double>,std::shared_ptr<Interpolator<double>>>(),
         py::arg("field"), py::arg("Interpolator"))
    .def(py::init<std::vector<double>,std::shared_ptr<Grid<double>>,
                  std::shared_ptr<Log>>(),
         py::arg("field"), py::arg("Grid"), py::arg("log"))
    .def(py::init<std::string,std::shared_ptr<Log>>(),
         py::arg("name"), py::arg("log"))
    .def(py::init<std::string,std::shared_ptr<Grid<double>>>(),
         py::arg("name"), py::arg("Grid"))
    .def(py::init<std::string,std::shared_ptr<Interpolator<double>>>(),
         py::arg("name"), py::arg("Interpolator"))
    .def(py::init<std::string,std::shared_ptr<Grid<double>>,
                  std::shared_ptr<Log>>(),
         py::arg("name"), py::arg("Grid"), py::arg("log"))
    .def(py::init<std::string,std::vector<double>>(),
          py::arg("name"), py::arg("field"))
    .def(py::init<std::string,std::vector<double>,std::shared_ptr<Log>>(),
          py::arg("name"), py::arg("field"), py::arg("log"))
    .def(py::init<std::string,std::vector<double>,
                  std::shared_ptr<Grid<double>>>(),
          py::arg("name"), py::arg("field"), py::arg("Grid"))
    .def(py::init<std::string,std::vector<double>,
                  std::shared_ptr<Interpolator<double>>>(),
          py::arg("name"), py::arg("field"), py::arg("Interpolator"))
    .def(py::init<std::string,std::vector<double>,
                  std::shared_ptr<Grid<double>>,
                  std::shared_ptr<Log>>(),
          py::arg("name"), py::arg("field"), py::arg("Grid"), py::arg("log"))
		//	Getters and Setters
    .def("get_type", &ScalarField<double>::getType)
		.def("get_field", &ScalarField<double>::getField,
		     py::return_value_policy::reference)
		.def("access_field", &ScalarField<double>::accessField,
		     py::return_value_policy::reference)
		.def("data", &ScalarField<double>::data,
		     py::return_value_policy::reference)
		.def("set_field", &ScalarField<double>::setField)
    .def("set_Interpolator", &ScalarField<double>::setInterpolator)
    .def("set_DiffEQ", &ScalarField<double>::setDiffEQ)
    .def("set_Integrator", &ScalarField<double>::setIntegrator)
    .def_property("Interpolator", &Field<double>::getInterpolator,
                                  &ScalarField<double>::setInterpolator)
    .def_property("DiffEQ", &Field<double>::getDiffEQ,
                                  &ScalarField<double>::setDiffEQ)
    .def_property("Integrator", &Field<double>::getIntegrator,
                                  &ScalarField<double>::setIntegrator)
    .def_property_readonly("type", &ScalarField<double>::getType)
		//	Operator overloads
		.def(py::self + py::self)
		.def(py::self - py::self)
		.def(py::self * py::self)
		.def(py::self / py::self)
		.def(py::self += py::self)
		.def(py::self -= py::self)
		.def(py::self *= py::self)
		.def(py::self /= py::self)
		.def("__getitem__", [](const ScalarField<double> &self, int i)
		{
			if (i < 0 || i >= self.getN()) {
				throw py::index_error("Index " + std::to_string(i) +
															" out of bounds for scalar field with "
															+ std::to_string(self.getN()) + " points!");
			}
			return self(i);
		}, py::is_operator())
		.def("__setitem__", [](ScalarField<double> &self,
													 int i, const double& val)
		{
			if (i < 0 || i >= self.getN()) {
				throw py::index_error("Index " + std::to_string(i) +
														" out of bounds for scalar field with "
														+ std::to_string(self.getN()) + " points!");
			}
			self(i) = val;
		}, py::is_operator())
		//  Special functions
    .def("construct_local_field_values", (Vector<double>
         (ScalarField<double>::*)(size_t))
         &ScalarField<double>::constructLocalFieldValues)
    .def("construct_local_field_values", (Vector<double>
        (ScalarField<double>::*)(const std::vector<double>&,size_t))
        &ScalarField<double>::constructLocalFieldValues)
    .def("derivative", (Vector<double> (ScalarField<double>::*)
         (const size_t, const size_t))
         &ScalarField<double>::derivative)
    .def("derivative", (double (ScalarField<double>::*)
         (const size_t, const size_t, const size_t))
         &ScalarField<double>::derivative)
    .def("derivative", (Vector<double> (ScalarField<double>::*)
         (const std::vector<double>&, const size_t))
         &ScalarField<double>::derivative)
    .def("derivative", (double (ScalarField<double>::*)
         (const std::vector<double>&, const size_t, const size_t))
         &ScalarField<double>::derivative)
    .def("field_derivative", (std::vector<Vector<double>>
         (ScalarField<double>::*)(const size_t))
         &ScalarField<double>::fieldDerivative)
    .def("field_derivative", (std::vector<double>
         (ScalarField<double>::*)(const size_t, const size_t))
         &ScalarField<double>::fieldDerivative)
    .def("log_output", [](ScalarField<double>& self)
    {
    	std::string out = self.getLog()->getOutput();
    	py::print(out);
    })
    .def("log_output", [](ScalarField<double>& self, size_t i)
    {
    	std::string out = self.getLog()->getOutput(i);
    	py::print(out);
    })
    //    __len__
    .def("__len__", &ScalarField<double>::getN)
    //    __repr__
    .def("__repr__", [](const ScalarField<double> &f) {
      return "<etraj.ScalarField<double> named '" + f.getName() + "'>";
    })
    //    __str__
    .def("__str__", [](const ScalarField<double> &f) {
      std::stringstream s;
			s << &f;
      std::string sum;
      sum += "++++++++++++++++++++++++++++++++++++++++++++++++++++";
			sum += "\n<etraj.ScalarField<double> ref at " + s.str() + ">";
      sum += "\n---------------------------------------------------";
      sum += "\n<ET::ScalarField<double> object at " + address_to_string(f) + ">";
      sum += "\n---------------------------------------------------";
      sum += "\n   name: '" + f.getName() + "'";
      sum += "\n    dim: " + std::to_string(f.getDim());
      sum += "\n      N: " + std::to_string(f.getN());
      sum += "\n---------------------------------------------------";
      sum += "\n        Grid at: " + address_to_string(*f.getGrid()) + ",";
      sum += "\n         ref at: " + address_to_string(f.getGrid());
      sum += "\n---------------------------------------------------";
      sum += "\nInterpolator at: " + address_to_string(*f.getInterpolator()) + ",";
      sum += "\n         ref at: " + address_to_string(f.getInterpolator());
      sum += "\n---------------------------------------------------";
      sum += "\n      DiffEQ at: " + address_to_string(*f.getDiffEQ()) + ",";
      sum += "\n         ref at: " + address_to_string(f.getDiffEQ());
      sum += "\n---------------------------------------------------";
      sum += "\n  Integrator at: " + address_to_string(*f.getIntegrator()) + ",";
      sum += "\n         ref at: " + address_to_string(f.getIntegrator());
      sum += "\n---------------------------------------------------";
      sum += "\n      Logger at: " + address_to_string(*f.getLog()) + ",";
      sum += "\n         ref at: " + address_to_string(f.getLog());
      sum += "\n++++++++++++++++++++++++++++++++++++++++++++++++++++";
      return sum;
    })
    ;
	//----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  InterpolatorType enum
  //----------------------------------------------------------------------------
  py::enum_<InterpolatorType>(m, "InterpolatorType")
    .value("default", InterpolatorType::DEFAULT)
    .value("local_taylor", InterpolatorType::LTE)
    .value("radial_basis", InterpolatorType::RBF)
    ;
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
	//	LSDriver enum
	//----------------------------------------------------------------------------
	py::enum_<LSDriver>(m, "LSDriver")
		.value("xGELS", LSDriver::xGELS)
		.value("xGELSY", LSDriver::xGELSY)
		.value("xGELSD", LSDriver::xGELSD)
		.value("xGELSS", LSDriver::xGELSS)
		;
	//----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  SolverType enum
  //----------------------------------------------------------------------------
  py::enum_<SolverType>(m, "SolverType")
    .value("LS", SolverType::LS)
    .value("MLS", SolverType::MLS)
    .value("WMLS", SolverType::WMLS)
    ;
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  WeightMatrixType enum
  //----------------------------------------------------------------------------
  py::enum_<WeightMatrixType>(m, "WeightMatrixType")
    .value("gaussian", WeightMatrixType::GAUSSIAN)
    ;
  //----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Interpolator class
	//----------------------------------------------------------------------------
	py::class_<Interpolator<double>,
	           std::shared_ptr<Interpolator<double>>>(m, "Interpolator")
		.def(py::init<>())
		.def(py::init<std::string>(),
         py::arg("name"))
		.def(py::init<std::shared_ptr<Log>>(),
         py::arg("log"))
    .def(py::init<std::string,std::shared_ptr<Log>>(),
         py::arg("name"), py::arg("log"))
		.def(py::init<std::shared_ptr<Grid<double>>>(),
         py::arg("Grid"))
    .def(py::init<std::string,std::shared_ptr<Grid<double>>>(),
         py::arg("name"), py::arg("Grid"))
    .def(py::init<std::shared_ptr<Grid<double>>,std::shared_ptr<Log>>(),
         py::arg("Grid"), py::arg("log"))
    .def(py::init<std::string,std::shared_ptr<Grid<double>>,
                  std::shared_ptr<Log>>(),
         py::arg("name"), py::arg("Grid"), py::arg("log"))
    //  Getters and Setters
    .def("get_name", &Interpolator<double>::getName)
    .def("get_Grid", &Interpolator<double>::getGrid)
    .def("get_Field", &Interpolator<double>::getField)
    .def("get_log", &Interpolator<double>::getLog)
    .def("get_flag", &Interpolator<double>::getFlag)
    .def("get_info", &Interpolator<double>::getInfo)
    .def("get_type", &Interpolator<double>::getInterpolatorType)
		.def("get_lsdriver", &Interpolator<double>::getLSDriver)
    .def("get_solvertype", &Interpolator<double>::getSolverType)
		.def("set_name", &Interpolator<double>::setName)
    .def("set_Grid", &Interpolator<double>::setGrid)
    .def("set_Field", &Interpolator<double>::setField)
    .def("set_log", &Interpolator<double>::setLog)
    .def("set_flag", &Interpolator<double>::setFlag)
    .def("set_info", &Interpolator<double>::setInfo)
		.def("set_lsdriver", &Interpolator<double>::setLSDriver)
    .def("set_solvertype", &Interpolator<double>::setSolverType)
    .def_property("name", &Interpolator<double>::getName,
                          &Interpolator<double>::setName)
    .def_property("Grid", &Interpolator<double>::getGrid,
                          &Interpolator<double>::setGrid)
    .def_property("Field", &Interpolator<double>::getField,
                           &Interpolator<double>::setField)
    .def_property("log", &Interpolator<double>::getLog,
                         &Interpolator<double>::setLog)
    .def_property("flag", &Interpolator<double>::setFlag,
                          &Interpolator<double>::setFlag)
    .def_property("info", &Interpolator<double>::getInfo,
                          &Interpolator<double>::setInfo)
    .def_property("lsdriver", &Interpolator<double>::getLSDriver,
                              &Interpolator<double>::setLSDriver)
    .def_property("solvertype", &Interpolator<double>::getSolverType,
                                &Interpolator<double>::setSolverType)
    .def_property_readonly("type", &Interpolator<double>::getInterpolatorType)
		.def("output", [](Interpolator<double>& self)
		{
			std::string out = self.getLog()->getOutput();
			py::print(out);
		})
		.def("output", [](Interpolator<double>& self, size_t i)
		{
			std::string out = self.getLog()->getOutput(i);
			py::print(out);
		})
    //    __str__
    .def("__str__", [](const Interpolator<double> &inter) {
     std::stringstream s;
     s << &inter;
     std::string sum;
     sum += "++++++++++++++++++++++++++++++++++++++++++++++++++++";
     sum += "\n<etraj.Inteprolator<double> ref at " + s.str() + ">";
     sum += "\n---------------------------------------------------";
     sum += "\n<ET::Inteprolator<double> object at " + address_to_string(inter) + ">";
     sum += "\n---------------------------------------------------";
     sum += "\n      name: '" + inter.getName() + "'";
     sum += "\nSolverType: " + SolverTypeNameMap[inter.getSolverType()];
     sum += "\n  LSDriver: " + LSDriverNameMap[inter.getLSDriver()];
     sum += "\n---------------------------------------------------";
     sum += "\n        Grid at: " + address_to_string(*inter.getGrid()) + ",";
     sum += "\n         ref at: " + address_to_string(inter.getGrid());
     sum += "\n---------------------------------------------------";
     sum += "\n       Field at: " + address_to_string(*inter.getGrid()) + ",";
     sum += "\n         ref at: " + address_to_string(inter.getGrid());
     sum += "\n---------------------------------------------------";
     sum += "\n      Logger at: " + address_to_string(*inter.getLog()) + ",";
     sum += "\n         ref at: " + address_to_string(inter.getLog());
     sum += "\n++++++++++++++++++++++++++++++++++++++++++++++++++++";
     return sum;
    })
    ;
	//----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Search Scheme enum
  //----------------------------------------------------------------------------
  py::enum_<SearchScheme>(m, "SearchScheme")
    .value("nearest_neighbors", SearchScheme::NEAREST_NEIGHBORS)
    .value("radius", SearchScheme::RADIUS)
    ;
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Local Taylor Interpolator
  //----------------------------------------------------------------------------
  py::class_<LocalTaylorInterpolator<double>, Interpolator<double>,
             std::shared_ptr<LocalTaylorInterpolator<double>>>
             (m, "LocalTaylorInterpolator")
    .def(py::init<>())
    .def(py::init<std::string>(),
        py::arg("name"))
    .def(py::init<std::shared_ptr<Log>>(),
        py::arg("log"))
    .def(py::init<std::string,std::shared_ptr<Log>>(),
        py::arg("name"), py::arg("log"))
    .def(py::init<std::shared_ptr<UGrid<double>>>(),
        py::arg("UGrid"))
    .def(py::init<std::string,std::shared_ptr<UGrid<double>>>(),
        py::arg("name"), py::arg("UGrid"))
    .def(py::init<std::shared_ptr<UGrid<double>>,std::shared_ptr<Log>>(),
        py::arg("UGrid"), py::arg("log"))
    .def(py::init<std::string,std::shared_ptr<UGrid<double>>,
                 std::shared_ptr<Log>>(),
        py::arg("name"), py::arg("UGrid"), py::arg("log"))
    //  Getters and Setters
    .def("get_k", &LocalTaylorInterpolator<double>::getK)
    .def("get_n", &LocalTaylorInterpolator<double>::getN)
    .def("get_radius", &LocalTaylorInterpolator<double>::getRadius)
    .def("get_search_scheme", &LocalTaylorInterpolator<double>::getSearchScheme)
    .def("set_k", &LocalTaylorInterpolator<double>::setK)
    .def("set_n", &LocalTaylorInterpolator<double>::setN)
    .def("set_radius", &LocalTaylorInterpolator<double>::setRadius)
    .def("set_search_scheme", &LocalTaylorInterpolator<double>::setSearchScheme)
    .def_property("k", &LocalTaylorInterpolator<double>::getK,
                       &LocalTaylorInterpolator<double>::setK)
    .def_property("n", &LocalTaylorInterpolator<double>::getN,
                       &LocalTaylorInterpolator<double>::setN)
    .def_property("radius", &LocalTaylorInterpolator<double>::getRadius,
                            &LocalTaylorInterpolator<double>::setRadius)
    .def_property("search_scheme", &LocalTaylorInterpolator<double>::getSearchScheme,
                                   &LocalTaylorInterpolator<double>::setSearchScheme)
    //  Taylor Matrices
    .def("construct_local_taylor_matrix", (Matrix<double>
         (LocalTaylorInterpolator<double>::*)
         (const size_t))
         &LocalTaylorInterpolator<double>::constructLocalTaylorMatrix)
    .def("construct_local_taylor_matrix", (Matrix<double>
         (LocalTaylorInterpolator<double>::*)
         (const std::vector<double>&))
         &LocalTaylorInterpolator<double>::constructLocalTaylorMatrix)
    .def("construct_local_taylor_matrix", (Matrix<double>
         (LocalTaylorInterpolator<double>::*)
         (const size_t, const size_t))
         &LocalTaylorInterpolator<double>::constructLocalTaylorMatrix)
    .def("construct_local_taylor_matrix", (Matrix<double>
         (LocalTaylorInterpolator<double>::*)
         (const std::vector<double>&, size_t))
         &LocalTaylorInterpolator<double>::constructLocalTaylorMatrix)
    .def("construct_local_taylor_matrix", (Matrix<double>
         (LocalTaylorInterpolator<double>::*)
         (const size_t, const double))
         &LocalTaylorInterpolator<double>::constructLocalTaylorMatrix)
    .def("construct_local_taylor_matrix", (Matrix<double>
         (LocalTaylorInterpolator<double>::*)
         (const std::vector<double>&, const double))
         &LocalTaylorInterpolator<double>::constructLocalTaylorMatrix)
    .def("construct_local_taylor_matrix", (Matrix<double>
         (LocalTaylorInterpolator<double>::*)
         (const size_t, const size_t, const size_t))
         &LocalTaylorInterpolator<double>::constructLocalTaylorMatrix)
    .def("construct_local_taylor_matrix", (Matrix<double>
         (LocalTaylorInterpolator<double>::*)
         (const std::vector<double>&, const size_t, const size_t))
         &LocalTaylorInterpolator<double>::constructLocalTaylorMatrix)
    .def("construct_local_taylor_matrix", (Matrix<double>
         (LocalTaylorInterpolator<double>::*)
         (const size_t, const double, const size_t))
         &LocalTaylorInterpolator<double>::constructLocalTaylorMatrix)
    .def("construct_local_taylor_matrix", (Matrix<double>
         (LocalTaylorInterpolator<double>::*)
         (const std::vector<double>&, const double, const size_t))
         &LocalTaylorInterpolator<double>::constructLocalTaylorMatrix)
    .def("derivative", (Vector<double> (LocalTaylorInterpolator<double>::*)
        (const size_t, const size_t))
        &LocalTaylorInterpolator<double>::derivative)
    .def("derivative", (double (LocalTaylorInterpolator<double>::*)
        (const size_t, const size_t, const size_t))
        &LocalTaylorInterpolator<double>::derivative)
    .def("derivative", (Vector<double> (LocalTaylorInterpolator<double>::*)
        (const std::vector<double>&, const size_t))
        &LocalTaylorInterpolator<double>::derivative)
    .def("derivative", (double (LocalTaylorInterpolator<double>::*)
        (const std::vector<double>&, const size_t, const size_t))
        &LocalTaylorInterpolator<double>::derivative)
    .def("field_derivative", (std::vector<Vector<double>>
         (LocalTaylorInterpolator<double>::*)(const size_t))
         &LocalTaylorInterpolator<double>::fieldDerivative)
    .def("field_derivative", (std::vector<double>
         (LocalTaylorInterpolator<double>::*)(const size_t, const size_t))
         &LocalTaylorInterpolator<double>::fieldDerivative)
    .def("scalarFieldDerivative", (Vector<double>
        (LocalTaylorInterpolator<double>::*)
        (const size_t, const size_t))
        &LocalTaylorInterpolator<double>::scalarFieldDerivative)
    .def("scalarFieldDerivative", (double
        (LocalTaylorInterpolator<double>::*)
        (const size_t, const size_t, const size_t))
        &LocalTaylorInterpolator<double>::scalarFieldDerivative)
    .def("scalarFieldDerivative", (Vector<double>
        (LocalTaylorInterpolator<double>::*)
        (const std::vector<double>&, const size_t))
        &LocalTaylorInterpolator<double>::scalarFieldDerivative)
    .def("scalarFieldDerivative", (double
        (LocalTaylorInterpolator<double>::*)
        (const std::vector<double>&, const size_t, const size_t))
        &LocalTaylorInterpolator<double>::scalarFieldDerivative)
    ;
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
	//  RBF class
	//----------------------------------------------------------------------------
  py::class_<RadialBasisInterpolator<double>, Interpolator<double>,
             std::shared_ptr<RadialBasisInterpolator<double>>>
             (m, "RadialBasisInterpolator")
    .def(py::init<>())
    .def("get_type", &RadialBasisInterpolator<double>::getType)
    .def("get_shape", &RadialBasisInterpolator<double>::getShape)
    .def("get_logger", &RadialBasisInterpolator<double>::getLogger)
    .def("set_type", (void (RadialBasisInterpolator<double>::*)
        (std::string)) &RadialBasisInterpolator<double>::setType)
    .def("set_shape", &RadialBasisInterpolator<double>::setShape)
    .def("set_logger", &RadialBasisInterpolator<double>::setLogger)
    .def("construct_rbf_matrix",
         (Matrix<double> (RadialBasisInterpolator<double>::*)
    		  (const std::shared_ptr<UGrid<double>>,
          const std::vector<size_t>,size_t))
    		 &RadialBasisInterpolator<double>::constructRBFMatrix)
    .def("construct_rbfd_matrix",
        (Matrix<double> (RadialBasisInterpolator<double>::*)
    		  (const std::shared_ptr<UGrid<double>>,
         const std::vector<size_t>,size_t,size_t))
    		 &RadialBasisInterpolator<double>::constructRBFFirstDerivativeMatrix)
    .def("construct_rbf_matrix",
        (Matrix<double> (RadialBasisInterpolator<double>::*)
    		  (const std::shared_ptr<UGrid<double>>))
    		 &RadialBasisInterpolator<double>::constructRBFMatrix)
    .def("construct_rbf_vector",
       (Vector<double> (RadialBasisInterpolator<double>::*)
    		  (const std::shared_ptr<UGrid<double>>,std::vector<double>))
    		 &RadialBasisInterpolator<double>::constructRBFVector)
    .def("construct_rbfd_matrix",
       (Matrix<double> (RadialBasisInterpolator<double>::*)
    		  (const std::shared_ptr<UGrid<double>>,size_t))
    		 &RadialBasisInterpolator<double>::constructRBFFirstDerivativeMatrix)
    .def("construct_rbfd_vector",
      (Vector<double> (RadialBasisInterpolator<double>::*)
    		  (const std::shared_ptr<UGrid<double>>,std::vector<double>))
    		 &RadialBasisInterpolator<double>::constructRBFFirstDerivativeVector)
    .def("output", [](RadialBasisInterpolator<double>& self)
		{
			std::string out = self.getLogger()->getOutput();
			py::print(out);
		})
		.def("output", [](RadialBasisInterpolator<double>& self, size_t i)
		{
			std::string out = self.getLogger()->getOutput(i);
			py::print(out);
		})
    //	print functionality
		.def("__repr__", [](RadialBasisInterpolator<double> &field)
		{
			std::stringstream s;
			s << &field;
			std::string res = "++++++++++++++++++++++++++++++++++++++++++++++++++++";
			res += "\n<etraj.RadialBasisInterpolator ref at " + s.str() + ">\n";
			return res + field.summary();
		})
    ;
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Interpolant Base class
  //----------------------------------------------------------------------------
  py::class_<Interpolant<double>,
             std::shared_ptr<Interpolant<double>>>(m, "Interpolant")
    .def(py::init<>())
    .def("get_name", &Interpolant<double>::getName,
         "Get the name of the Interpolant")
    .def("get_dim", &Interpolant<double>::getDim)
    .def("get_ranges", &Interpolant<double>::getRanges)
    .def("get_range", &Interpolant<double>::getRange)
    .def("get_range_min", &Interpolant<double>::getRangeMin)
    .def("get_range_max", &Interpolant<double>::getRangeMax)
    .def("set_name", &Interpolant<double>::setName)
    .def("set_dim", &Interpolant<double>::setDim)
    .def("set_ranges", &Interpolant<double>::setRanges)
    .def("set_range", &Interpolant<double>::setRange)
    .def("set_range_min", &Interpolant<double>::setRangeMin)
    .def("set_range_max", &Interpolant<double>::setRangeMax)
    .def_property("name", &Interpolant<double>::getName,
                          &Interpolant<double>::setName)
    .def_property("dim", &Interpolant<double>::getDim,
                         &Interpolant<double>::setDim)
    .def_property("ranges", &Interpolant<double>::getRanges,
                            &Interpolant<double>::setRanges)
    //  Access operators
    .def("__call__", [](Interpolant<double> &self, std::vector<double>& p)
    {
      return self(p);
    })
    .def("__call__", [](Interpolant<double> &self,
                        std::vector<std::vector<double>>& p)
    {
      return self(p);
    })
    //	this will allow the user to get a range.
		.def("__getitem__", [](Interpolant<double> &self, int i)
		{
			if (i < 0 || i >= self.getDim()) {
				throw py::index_error("Index " + std::to_string(i) +
															" out of bounds for array with "
															+ std::to_string(self.getDim())
															+ " rows!");
			}
			return self[i];
		}, py::is_operator())
    //	this will allow the user to get a range.
		.def("__setitem__", [](Interpolant<double> &self, int i,
                           std::vector<double> range)
		{
			if (i < 0 || i >= self.getDim()) {
				throw py::index_error("Index " + std::to_string(i) +
															" out of bounds for array with "
															+ std::to_string(self.getDim())
															+ " rows!");
			}
			self[i] = range;
		}, py::is_operator())
    //  Derivatives
    .def("d", (double (Interpolant<double>::*)
         (const std::vector<double>&))
         &Interpolant<double>::d)
    .def("d", (std::vector<double> (Interpolant<double>::*)
        (const std::vector<std::vector<double>>&))
        &Interpolant<double>::d)
    .def("dd", (double (Interpolant<double>::*)
         (const std::vector<double>&))
         &Interpolant<double>::dd)
    .def("dd", (std::vector<double> (Interpolant<double>::*)
        (const std::vector<std::vector<double>>&))
        &Interpolant<double>::dd)
    ;
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Local Taylor Interpolant Base class
  //----------------------------------------------------------------------------
  py::class_<LocalTaylorInterpolant<double>, Interpolant<double>,
             std::shared_ptr<LocalTaylorInterpolant<double>>>
             (m, "LocalTaylorInterpolant")
    .def(py::init<>())
    .def("get_n", &LocalTaylorInterpolant<double>::get_n)
    .def("get_expansion_point",
         &LocalTaylorInterpolant<double>::getExpansionPoint)
    .def("get_expansion_coefficients",
         &LocalTaylorInterpolant<double>::getExpansionCoefficients)
    .def("set_n", &LocalTaylorInterpolant<double>::set_n)
    .def("set_expansion_point",
         &LocalTaylorInterpolant<double>::setExpansionPoint)
    .def("set_expansion_coefficients",
         &LocalTaylorInterpolant<double>::setExpansionCoefficients)
    .def_property("n", &LocalTaylorInterpolant<double>::get_n,
                       &LocalTaylorInterpolant<double>::set_n)
    .def_property("expansion_point",
                  &LocalTaylorInterpolant<double>::getExpansionPoint,
                  &LocalTaylorInterpolant<double>::setExpansionPoint)
    .def_property("expansion_coefficients",
                  &LocalTaylorInterpolant<double>::getExpansionCoefficients,
                  &LocalTaylorInterpolant<double>::setExpansionCoefficients)
    ;
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Inteprolant generating functions
  //----------------------------------------------------------------------------
  m.def("create_local_taylor_interpolant",
        &createLocalTaylorInterpolant<double>);
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
	//	DiffEQ Base class
	//----------------------------------------------------------------------------
  py::class_<DiffEQ<double>,
             std::shared_ptr<DiffEQ<double>>>(m, "DiffEQ")
    .def(py::init<>())
    ;
  //----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  Integrator class
	//----------------------------------------------------------------------------
	py::class_<Integrator<double>,
	           std::shared_ptr<Integrator<double>>>(m, "Integrator")
		.def(py::init<>())
		//.def("scalar_RK4_step", &Integrator<double>::scalarRK4Step)
		;
	//----------------------------------------------------------------------------


	// //----------------------------------------------------------------------------
	// //  Appproximator enum
	// //----------------------------------------------------------------------------
	// py::class_<VectorField<double>,
	//            std::shared_ptr<VectorField<double>>>(m, "VectorField")
	//   .def(py::init<>())
	// 	;
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //  Appproximator enum
	// //----------------------------------------------------------------------------
	// py::enum_<InterpolatorType>(m, "InterpolatorType")
	// 	.value("LS", InterpolatorType::LS)
	// 	.value("RBF", InterpolatorType::RBF)
	// 	;
	// //----------------------------------------------------------------------------



	// //----------------------------------------------------------------------------
	// //	Interpolator Parameters Struct
	// //----------------------------------------------------------------------------
	// py::class_<InterpolatorParams>(m, "InterpolatorParams")
	// 	.def(py::init<>())
	// 	.def_readwrite("k", &InterpolatorParams::k)
	// 	.def_readwrite("n", &InterpolatorParams::n)
	// 	;
	// //----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Monomial class
	//----------------------------------------------------------------------------
	py::class_<Monomial, std::shared_ptr<Monomial>>(m, "Monomial")
		.def(py::init<>())
		.def(py::init<size_t>())
		.def(py::init<size_t,size_t>())
		.def("get_dim", &Monomial::getDim,
		     py::return_value_policy::reference)
		.def("get_deg", &Monomial::getDeg,
		     py::return_value_policy::reference)
		.def("get_mono", (std::vector<std::vector<size_t>>
			   (Monomial::*)()) &Monomial::getMono,
				 py::return_value_policy::reference)
		.def("get_mono", (std::vector<std::vector<size_t>>
				 (Monomial::*)(size_t)) &Monomial::getMono,
				 py::return_value_policy::reference)
		.def("get_mono_factors", &Monomial::getMonoFactors,
		     py::return_value_policy::reference)
		.def("get_multiset_coeff", &Monomial::getMultisetCoefficient,
		     py::return_value_policy::reference)
		.def("get_taylor_index", &Monomial::getTaylorIndex,
		     py::return_value_policy::reference)
		.def("set_dim", &Monomial::setDim)
		.def("set_deg", &Monomial::setDeg)
		.def("generate_monomial", (void (Monomial::*)())
		     &Monomial::generateMonomial)
		.def("generate_monomial", (void (Monomial::*)(size_t))
		 		 &Monomial::generateMonomial)
		.def("taylor_monomial_expansion", (std::vector<double>
			   (Monomial::*)(const std::vector<double>&,
					                 const std::vector<double>&))
			   &Monomial::taylorMonomialExpansion)
		.def("taylor_monomial_expansion", (std::vector<double>
	 		   (Monomial::*)(const std::vector<double>&,
					                 const std::vector<double>&,size_t))
	 		   &Monomial::taylorMonomialExpansion)
	  //	print functionality
		.def("__repr__", [](Monomial &mon)
		{
				return mon.summary();
		})
		;
	//----------------------------------------------------------------------------


  //
  // //----------------------------------------------------------------------------
	// //	FirstOrderODE Base class
	// //----------------------------------------------------------------------------
  // py::class_<FirstOrderODE<double>, DiffEQ<double>,
  //            std::shared_ptr<FirstOrderODE<double>>>(m, "FirstOrderODE")
  //   .def(py::init<>())
  //   .def("set_scalarfield", &FirstOrderODE<double>::setScalarField)
  //   .def("set_vectorfield", &FirstOrderODE<double>::setVectorField)
  //   .def("get_scalarfield", &FirstOrderODE<double>::getScalarField)
  //   .def("get_vectorfield", &FirstOrderODE<double>::getVectorField)
  //   //.def("dt", &FirstOrderODE<double>::dt)
  //   ;
  // //----------------------------------------------------------------------------
  //
  // //----------------------------------------------------------------------------
	// //	SecondOrderODE Base class
	// //----------------------------------------------------------------------------
  // py::class_<SecondOrderODE<double>, DiffEQ<double>,
  //            std::shared_ptr<SecondOrderODE<double>>>(m, "SecondOrderODE")
  //   .def(py::init<>())
  //   .def("set_scalarfield", &SecondOrderODE<double>::setScalarField)
  //   .def("set_vectorfield", &SecondOrderODE<double>::setVectorField)
  //   .def("get_scalarfield", &SecondOrderODE<double>::getScalarField)
  //   .def("get_vectorfield", &SecondOrderODE<double>::getVectorField)
  //   //.def("dt", &SecondOrderODE<double>::dt)
  //   //.def("dt", &SecondOrderODE<double>::d2t)
  //   ;
  // //----------------------------------------------------------------------------
  //
  // //----------------------------------------------------------------------------
	// //	FirstOrderODE Examples
	// //----------------------------------------------------------------------------
  // //----------------------------------------------------------------------------
	// //	heat equation
	// //----------------------------------------------------------------------------
  // py::class_<HeatEquation<double>, FirstOrderODE<double>,
  //            std::shared_ptr<HeatEquation<double>>>(m, "HeatEquation")
  //   .def(py::init<>())
  //   .def("set_alpha", &HeatEquation<double>::setAlpha)
  //   .def("get_alpha", &HeatEquation<double>::getAlpha)
  //   .def("dt", (double (HeatEquation<double>::*)
  //        (std::vector<double>,double,double))
  //        &HeatEquation<double>::dt)
  //   .def("dt", (std::vector<double> (HeatEquation<double>::*)
  //       (std::vector<std::vector<double>>,double,std::vector<double>))
  //       &HeatEquation<double>::dt)
  //   ;
  // //----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Logger class
	//----------------------------------------------------------------------------
	py::class_<Log, std::shared_ptr<Log>>(m, "Log")
		.def(py::init<>())
		.def("get_output", (std::string (Log::*)(size_t)) &Log::getOutput)
		.def("get_output", (std::string (Log::*)()) &Log::getOutput)
		;
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Dynamical System class
	//----------------------------------------------------------------------------
	py::class_<DynamicalSystem<double>>(m, "DynamicalSystem")
		.def(py::init<>())
		.def(py::init<std::string>())
		;
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	various functions
	//----------------------------------------------------------------------------
	m.def("taylor_polynomial", &taylorPolynomial);
	//----------------------------------------------------------------------------



}
