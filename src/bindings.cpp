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
#include "utils.h"
#include "radialbasis.h"
#include "localtaylor.h"
#include "interpolator.h"
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
    .def(py::init<size_t>())
		.def(py::init<std::string, size_t>())
    .def(py::init<std::vector<double>>())
		.def(py::init<std::string, std::vector<double>>())
		.def(py::init<size_t, const double&>())
		.def(py::init<std::string, size_t, const double&>())
		.def(py::init([](py::buffer const b)
		{
    	py::buffer_info info = b.request();
    	if (info.format != py::format_descriptor<float>::format()
					|| info.ndim != 1)
			{
      	throw std::runtime_error("Incompatible buffer format!");
			}
      auto v = new Vector<double>(info.shape[0]);
      memcpy(v->getVec().data(), info.ptr,
						 sizeof(double) * (size_t) (v->getDim()));
      return v;
    }))
		.def("get_dim", &Vector<double>::getDim,
		     py::return_value_policy::reference)
		.def("get_vec", &Vector<double>::getVec,
		     py::return_value_policy::reference)
		.def("get_name", &Vector<double>::getName,
		     py::return_value_policy::reference)
		.def("set_dim", &Vector<double>::setDim)
		.def("set_vec", &Vector<double>::setVec)
		.def("set_name", &Vector<double>::setName)
		.def("dot", &Vector<double>::dot)
    .def_property("name", &Vector<double>::getName,
                  &Vector<double>::setName)
    .def_property("dim", &Vector<double>::getDim,
                  &Vector<double>::setDim)
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
		//	print functionality
		.def("__repr__", [](Vector<double> &vector)
		{
				return vector.summary();
		})
		;
		//--------------------------------------------------------------------------

		//--------------------------------------------------------------------------
		//	Instantiations of vectors
		//--------------------------------------------------------------------------
		//m.def("zeros", (Vector<double> (*)(size_t)) &zeros);
		//m.def("zeros", (Vector<double> (*)(size_t,size_t)) &zeros);
		//m.def("ones", (Vector<double> (*)(size_t)) &ones);
		//m.def("ones", (Vector<double> (*)(size_t,size_t)) &ones);
		//--------------------------------------------------------------------------

		//--------------------------------------------------------------------------
		//	Level 1 BLAS methods
		//--------------------------------------------------------------------------
		m.def("dswap", &DSWAP);
		m.def("dscal", &DSCAL);
		m.def("dcopy", (void (*)(Vector<double>&,
				  Vector<double>&)) &DCOPY);
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
    .def(py::init<size_t>())
    .def(py::init<std::string, size_t>())
    .def(py::init<size_t, size_t>())
    .def(py::init<std::string, size_t, size_t>())
    .def(py::init<size_t, size_t, const double&>())
    .def(py::init<std::string, size_t, size_t, const double&>())
    .def(py::init<size_t, std::vector<double>>())
    .def(py::init<std::string, size_t, std::vector<double>>())
    .def(py::init<size_t, size_t, std::vector<double>>())
    .def(py::init<std::string, size_t, size_t, std::vector<double>>())
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
		.def("get_singular_values", &Matrix<double>::getSingularValues,
		     py::return_value_policy::reference)
		//	setters
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
		.def("__setitem__", [](Matrix<double> &self,
													 std::tuple<size_t, size_t> ij,
													 const double& val)
		{
			size_t i, j;
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
		.def("__getitem__", [](Matrix<double> &self, int i)
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
		.def("__repr__", [](Matrix<double> &m)
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
    //  Getters and Setters
    .def("get_name", &Grid<double>::getName)
    .def("get_dim", &Grid<double>::getDim)
    .def("get_N", &Grid<double>::getN)
    .def("get_log", &Grid<double>::getLog)
    .def("set_name", &Grid<double>::setName)
    .def("set_dim", &Grid<double>::setDim)
    .def("set_N", &Grid<double>::setN)
    .def("set_log", &Grid<double>::setLog)
    .def_property("name", &Grid<double>::getName,
                  &Grid<double>::setName)
    .def_property("dim", &Grid<double>::getDim,
                  &Grid<double>::setDim)
    .def_property("N", &Grid<double>::getN,
                  &Grid<double>::setN)
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
      sum += "\n<ET::Grid<double> object at " + getMem(g) + ">";
      sum += "\n---------------------------------------------------";
      sum += "\n   name: '" + g.getName() + "'";
      sum += "\n    dim: " + std::to_string(g.getDim());
      sum += "\n      N: " + std::to_string(g.getN());
      sum += "\n---------------------------------------------------";
      sum += "\n Logger at: " + getMem(*g.getLog()) + ",";
      sum += "\n    ref at: " + getMem(g.getLog());
      sum += "\n++++++++++++++++++++++++++++++++++++++++++++++++++++";
      return sum;
    })
    ;
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
		.def(py::init<std::vector<double>>(), py::keep_alive<1, 2>())
		.def(py::init<std::vector<std::vector<double>>>(), py::keep_alive<1, 2>())

		.def(py::init<std::vector<double>, std::shared_ptr<Log>>(),
		     py::keep_alive<1, 2>())
		.def(py::init<std::vector<std::vector<double>>,
			            std::shared_ptr<Log>>(), py::keep_alive<1, 2>())
    //  Getters and Setters
		.def("get_ugrid", &UGrid<double>::getUGrid,
		     py::return_value_policy::reference)
		.def("get_neighbors", (std::vector<std::vector<size_t>>
				 (UGrid<double>::*)()) &UGrid<double>::getNeighbors,
				 py::return_value_policy::reference)
		.def("get_neighbors", [](UGrid<double>& self, size_t i)
		{
			if (i < 0 || i >= self.getN())
			{
				//######################################################################
				self.getLog()->ERROR("UGrid " + self.getName()
										+ ": Attempted to access neighbors array of size "
										+ std::to_string(self.getN()) + " with index "
										+ std::to_string(i));
				//######################################################################
				throw py::index_error("Index " + std::to_string(i) +
															" out of bounds for array with "
															+ std::to_string(self.getN())
															+ " points!");
			}
			return self.getNeighbors(i);
		}, py::return_value_policy::reference)
		.def("get_distances", &UGrid<double>::getDistances,
		     py::return_value_policy::reference)
		.def("get_neighbors_radius", &UGrid<double>::getNeighborsRadius,
		     py::return_value_policy::reference)
		.def("get_distances_radius", &UGrid<double>::getDistancesRadius,
		     py::return_value_policy::reference)
		// .def("get_logger", &UGrid<double>::getLogger,
		//      py::return_value_policy::reference)
    // .def_property("name", &UGrid<double>::getName,
    //               &UGrid<double>::setName)
    // .def_property("dim", &UGrid<double>::getDim,
    //               &UGrid<double>::setDim)
    // .def_property("N", &UGrid<double>::getN,
    //               &UGrid<double>::setN)
		.def("output", [](UGrid<double>& self)
		{
			std::string out = self.getLog()->getOutput();
			py::print(out);
		})
		.def("output", [](UGrid<double>& self, size_t i)
		{
			std::string out = self.getLog()->getOutput(i);
			py::print(out);
		})

		.def("set_ugrid", &UGrid<double>::setUGrid)
		//  Access operators
		.def("__getitem__", [](UGrid<double> &self,
					std::tuple<int, int> ij)
		{
			int i, j;
			std::tie(i, j) = ij;
			if (i < 0 || i >= self.getN())
			{
				if (i < 0)
				{
					//####################################################################
					self.getLog()->ERROR("UGrid " + self.getName()
											+ ": Attempted to access _ugrid array of size "
											+ std::to_string(self.getN()) + " with invalid value "
											+ std::to_string(i));
					//####################################################################
				}
				else if (i >= self.getN())
				{
					//####################################################################
					self.getLog()->ERROR("UGrid " + self.getName()
											+ ": Attempted to access _ugrid array of size "
											+ std::to_string(self.getN()) + " with index "
											+ std::to_string(i));
					//####################################################################
				}
				throw py::index_error("Index " + std::to_string(i) +
															" out of bounds for array with "
															+ std::to_string(self.getN())
															+ " points!");
			}
			if (j < 0 || j >= self.getDim())
			{
				if (j < 0)
				{
					//####################################################################
					self.getLog()->ERROR("UGrid " + self.getName()
											+ ": Attempted to access _ugrid array of dimension "
											+ std::to_string(self.getDim()) + " with invalid value "
											+ std::to_string(j));
					//####################################################################
				}
				else if (j >= self.getDim())
				{
					//####################################################################
					self.getLog()->ERROR("UGrid " + self.getName()
											+ ": Attempted to access _ugrid array of dimension "
											+ std::to_string(self.getDim()) + " with index "
											+ std::to_string(j));
					//####################################################################
				}
				throw py::index_error("Index " + std::to_string(j) +
														  " out of bounds for array with dimension "
														  + std::to_string(self.getDim()));
			}
			return self(i,j);
		}, py::is_operator())
		//--------------------------------------------------------------------------
		.def("__setitem__", [](UGrid<double> &self,
					std::tuple<int, int> ij, const double& val)
		{
			int i, j;
			std::tie(i, j) = ij;
			if (i < 0 || i >= self.getN())
			{
				if (i < 0)
				{
					//####################################################################
					self.getLog()->ERROR("UGrid " + self.getName()
											+ ": Attempted to access _ugrid array of size "
											+ std::to_string(self.getN()) + " with invalid value "
											+ std::to_string(i));
					//####################################################################
				}
				else if (i >= self.getN())
				{
					//####################################################################
					self.getLog()->ERROR("UGrid " + self.getName()
											+ ": Attempted to access _ugrid array of size "
											+ std::to_string(self.getN()) + " with index "
											+ std::to_string(i));
					//####################################################################
				}
				throw py::index_error("Index " + std::to_string(i) +
															" out of bounds for array with "
															+ std::to_string(self.getN())
															+ " points!");
			}
			if (j < 0 || j >= self.getDim())
			{
				if (j < 0)
				{
					//####################################################################
					self.getLog()->ERROR("UGrid " + self.getName()
											+ ": Attempted to access _ugrid array of dimension "
											+ std::to_string(self.getDim()) + " with invalid value "
											+ std::to_string(j));
					//####################################################################
				}
				else if (j >= self.getDim())
				{
					//####################################################################
					self.getLog()->ERROR("UGrid " + self.getName()
											+ ": Attempted to access _ugrid array of dimension "
											+ std::to_string(self.getDim()) + " with index "
											+ std::to_string(j));
					//####################################################################
				}
				throw py::index_error("Index " + std::to_string(j) +
														  " out of bounds for array with dimension "
														  + std::to_string(self.getDim()));
			}
			self(i,j) = val;
		}, py::is_operator())
		//--------------------------------------------------------------------------
		.def("__getitem__", [](UGrid<double> &self, int i)
		{
			if (i < 0 || i >= self.getN())
			{
				if (i < 0)
				{
					//####################################################################
					self.getLog()->ERROR("UGrid " + self.getName()
											+ ": Attempted to access _ugrid array of size "
											+ std::to_string(self.getN()) + " with invalid value "
											+ std::to_string(i));
					//####################################################################
				}
				else if (i >= self.getN())
				{
					//####################################################################
					self.getLog()->ERROR("UGrid " + self.getName()
											+ ": Attempted to access _ugrid array of size "
											+ std::to_string(self.getN()) + " with index "
											+ std::to_string(i));
					//####################################################################
				}
				throw py::index_error("Index " + std::to_string(i) +
															" out of bounds for array with "
															+ std::to_string(self.getN())
															+ " rows!");
			}
			else
			{
				return self(i);
			}
		}, py::is_operator())
		//--------------------------------------------------------------------------
		.def("__setitem__", [](UGrid<double> &self, int i,
			                     std::vector<double> x)
		{
			if (i < 0 || i >= self.getN())
			{
				if (i < 0)
				{
					//####################################################################
					self.getLog()->ERROR("UGrid " + self.getName()
											+ ": Attempted to access _ugrid array of size "
											+ std::to_string(self.getN()) + " with invalid value "
											+ std::to_string(i));
					//####################################################################
				}
				else if (i >= self.getN())
				{
					//####################################################################
					self.getLog()->ERROR("UGrid " + self.getName()
											+ ": Attempted to access _ugrid array of size "
											+ std::to_string(self.getN()) + " with index "
											+ std::to_string(i));
					//####################################################################
				}
				throw py::index_error("Index " + std::to_string(i) +
															" out of bounds for array with "
															+ std::to_string(self.getN())
															+ " rows!");
			}
			else
			{
				self(i) = x;
			}
		}, py::is_operator())
		//--------------------------------------------------------------------------
		//  KDTree functions
		//--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
		//  Query neighbors for points on the grid
		//--------------------------------------------------------------------------
		.def("query_neighbors", (void
		     (UGrid<double>::*)(size_t k))
				 &UGrid<double>::queryNeighbors)
		.def("query_neighbors", (std::vector<std::vector<size_t>>
 		     (UGrid<double>::*)(const std::vector<std::vector<double>>&,
					                      size_t k))
 				 &UGrid<double>::queryNeighbors)
		.def("query_radius", &UGrid<double>::queryRadius)
    //--------------------------------------------------------------------------
		//  Query neighbors for arbitrary points
		//--------------------------------------------------------------------------
    .def("query_neighbors", (std::vector<size_t>
		     (UGrid<double>::*)(const std::vector<double>&, size_t))
				 &UGrid<double>::queryNeighbors)
    .def("query_distances", (std::vector<size_t>
         (UGrid<double>::*)(const std::vector<double>&, size_t))
         &UGrid<double>::queryDistances)
    .def("query_neighbors", (std::vector<std::vector<size_t>>
         (UGrid<double>::*)(const std::vector<std::vector<double>>&,
          size_t))
    		 &UGrid<double>::queryNeighbors)
		//--------------------------------------------------------------------------
		;
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	ScalarField class
	//----------------------------------------------------------------------------
	py::class_<ScalarField<double>,
	           std::shared_ptr<ScalarField<double>>>(m, "ScalarField")
		//--------------------------------------------------------------------------
		//	Constructors
		//--------------------------------------------------------------------------
		.def(py::init<>())
		.def(py::init<std::shared_ptr<UGrid<double>>>())
		//.def(py::init<std::string,std::shared_ptr<UGrid<double>>>())
		//.def(py::init<std::shared_ptr<UGrid<double>>,std::vector<double>>())
		//.def(py::init<std::string,std::shared_ptr<UGrid<double>>,
		//	            std::vector<double>>())
		//--------------------------------------------------------------------------
		//	Constructors with shared loggers
		//--------------------------------------------------------------------------
		.def(py::init<std::shared_ptr<Log>>())
		.def(py::init<std::shared_ptr<UGrid<double>>,
			            std::shared_ptr<Log>>())
		//.def(py::init<std::string,std::shared_ptr<UGrid<double>>,
		//	            std::shared_ptr<Log>>())
		//.def(py::init<std::shared_ptr<UGrid<double>>,std::vector<double>,
		//	            std::shared_ptr<Log>>())
		//.def(py::init<std::string,std::shared_ptr<UGrid<double>>,
		//	            std::vector<double>,std::shared_ptr<Log>>())
		//--------------------------------------------------------------------------
		//	Getters and Setters
		//--------------------------------------------------------------------------
		.def("get_field", &ScalarField<double>::getField,
		     py::return_value_policy::reference)
		.def("access_field", &ScalarField<double>::accessField,
		     py::return_value_policy::reference)
		.def("data", &ScalarField<double>::data,
		     py::return_value_policy::reference)
		.def("get_name", &ScalarField<double>::getName,
		     py::return_value_policy::reference)
		.def("get_N", &ScalarField<double>::getN,
		     py::return_value_policy::reference)
	  .def("get_dim", &ScalarField<double>::getDim,
		     py::return_value_policy::reference)
		.def("get_interpolator", &ScalarField<double>::getInterpolator,
		     py::return_value_policy::reference)
		.def("get_integrator", &ScalarField<double>::getIntegrator,
	       py::return_value_policy::reference)
		// .def("get_logger", &ScalarField<double>::getLogger,
 		//      py::return_value_policy::reference)
		.def("set_ugrid", &ScalarField<double>::setUGrid)
		.def("set_field", &ScalarField<double>::setField)
		.def("set_name", &ScalarField<double>::setName)
		//.def("set_approx_type", &ScalarField<double>::setInterpolatorType)
		// .def("output", [](ScalarField<double>& self)
		// {
		// 	std::string out = self.getLogger()->getOutput();
		// 	py::print(out);
		// })
		// .def("output", [](ScalarField<double>& self, size_t i)
		// {
		// 	std::string out = self.getLogger()->getOutput(i);
		// 	py::print(out);
		// })
		//--------------------------------------------------------------------------
		//  Direct exposure
		//--------------------------------------------------------------------------
		//.def_readwrite("name", &ScalarField<double>::getName)
		//--------------------------------------------------------------------------

		//--------------------------------------------------------------------------
		//	Operator overloads
		//--------------------------------------------------------------------------
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
			if (i < 0 || i >= self.getN())
			{
				throw py::index_error("Index " + std::to_string(i) +
															" out of bounds for scalar field with "
															+ std::to_string(self.getN()) + " points!");
			}
			return self(i);
		}, py::is_operator())
		.def("__setitem__", [](ScalarField<double> &self,
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
		//--------------------------------------------------------------------------
		//  Derivatives of the field
		//--------------------------------------------------------------------------
		// .def("gradient", &ScalarField<double>::gradient,
		// 		 py::return_value_policy::reference)
    // //--------------------------------------------------------------------------
 		// //  Derivatives along the entire grid
 		// //--------------------------------------------------------------------------
		// .def("derivative",
		//      (std::vector<std::vector<double>> (ScalarField<double>::*)
		// 	    (size_t)) &ScalarField<double>::derivative,
	  //      py::return_value_policy::reference)
	 	// .def("derivative",
		//      (std::vector<double> (ScalarField<double>::*)
		// 	    (size_t, size_t)) &ScalarField<double>::derivative,
	  //      py::return_value_policy::reference)
    // //--------------------------------------------------------------------------
 		// //  Derivatives at points on the grid
 		// //--------------------------------------------------------------------------
	  // .def("derivative_point",
		//      (std::vector<double> (ScalarField<double>::*)
		// 	    (size_t, size_t))
    //       &ScalarField<double>::derivativePoint,
	  //      py::return_value_policy::reference)
	  // .def("derivative_point",
		//      (double (ScalarField<double>::*)
		// 	    (size_t, size_t, size_t))
		// 		 &ScalarField<double>::derivativePoint,
	  //      py::return_value_policy::reference)
    // .def("derivative_point",
		//      (double (ScalarField<double>::*)
		// 	    (size_t, std::vector<size_t>))
		// 		 &ScalarField<double>::derivativePoint,
	  //      py::return_value_policy::reference)
    // //--------------------------------------------------------------------------
 		// //  Derivatives at arbitrary points
 		// //--------------------------------------------------------------------------
    // .def("derivative_point",
		//      (std::vector<double> (ScalarField<double>::*)
		// 	    (std::vector<double>, size_t))
    //       &ScalarField<double>::derivativePoint,
	  //      py::return_value_policy::reference)
	  // .def("derivative_point",
		//      (double (ScalarField<double>::*)
		// 	    (std::vector<double>, size_t, size_t))
		// 		 &ScalarField<double>::derivativePoint,
	  //      py::return_value_policy::reference)
    //  .def("derivative_point",
 		//      (double (ScalarField<double>::*)
 		// 	    (std::vector<double>, std::vector<size_t>))
 		// 		 &ScalarField<double>::derivativePoint,
 	  //      py::return_value_policy::reference)
		//--------------------------------------------------------------------------
		//	print functionality
		.def("__repr__", [](ScalarField<double> &field)
		{
			std::stringstream s;
			s << &field;
			std::string res = "++++++++++++++++++++++++++++++++++++++++++++++++++++";
			res += "\n<etraj.ScalarField ref at " + s.str() + ">\n";
			return res + field.summary();
		})
		;
	//----------------------------------------------------------------------------

	// //----------------------------------------------------------------------------
	// //	Scalarfield examples
	// //----------------------------------------------------------------------------
	// //----------------------------------------------------------------------------
	// //	Wave EQ in 1D
	// //----------------------------------------------------------------------------
	// py::class_<WaveEQ1D<double>,
	//            ScalarField<double>,
	// 					 std::shared_ptr<WaveEQ1D<double>>>(m,"WaveEQ1D")
	// 	.def(py::init<>())
	// 	.def(py::init<std::shared_ptr<UGrid<double>>>())
	// 	.def(py::init<std::shared_ptr<UGrid<double>>,double,double,double>())
	// 	.def("get_A", &WaveEQ1D<double>::getA)
	// 	.def("set_A", &WaveEQ1D<double>::setA)
	// 	.def("get_k", &WaveEQ1D<double>::getk)
	// 	.def("set_k", &WaveEQ1D<double>::setk)
	// 	.def("get_w", &WaveEQ1D<double>::getw)
	// 	.def("set_w", &WaveEQ1D<double>::setw)
	// 	;
	// //----------------------------------------------------------------------------
	// //----------------------------------------------------------------------------
	// //	Klein-Gordon field in 1D
	// //----------------------------------------------------------------------------
	// py::class_<KleinGordon1D<double>,
	//            ScalarField<double>,
	// 					 std::shared_ptr<KleinGordon1D<double>>>(m,"KleinGordon1D")
	// 	.def(py::init<>())
	// 	.def(py::init<std::shared_ptr<UGrid<double>>>())
	// 	.def(py::init<std::shared_ptr<UGrid<double>>,double>())
	// 	.def("get_mass", &KleinGordon1D<double>::getMass)
	// 	.def("set_mass", &KleinGordon1D<double>::setMass)
	// 	;
	// //----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Interpolator class
	//----------------------------------------------------------------------------
	py::class_<Interpolator<double>,
	           std::shared_ptr<Interpolator<double>>>(m, "Interpolator")
		//--------------------------------------------------------------------------
		//	Constructors
		//--------------------------------------------------------------------------
		.def(py::init<>())
		//.def(py::init<int>())
		//.def(py::init<std::string>())
		//--------------------------------------------------------------------------
		//	Constructors with shared loggers
		//--------------------------------------------------------------------------
		.def(py::init<std::shared_ptr<Log>>())
		//.def(py::init<int,std::shared_ptr<Log>>())
		//.def(py::init<std::string,std::shared_ptr<Log>>())
		//--------------------------------------------------------------------------
		//	Getters and setters
		//--------------------------------------------------------------------------
		// .def("get_approx_type", &Interpolator<double>::getInterpolatorType,
		//      py::return_value_policy::reference)
		// .def("get_approx_params", &Interpolator<double>::getInterpolatorParams,
		//      py::return_value_policy::reference)
		.def("get_lsdriver", &Interpolator<double>::getLSDriver,
		     py::return_value_policy::reference)
		.def("get_flag", &Interpolator<double>::getFlag,
		     py::return_value_policy::reference)
		.def("get_info", &Interpolator<double>::getInfo,
		     py::return_value_policy::reference)
		.def("get_logger", &Interpolator<double>::getLogger,
  	     py::return_value_policy::reference)
		// .def("set_approx_type", &Interpolator<double>::setInterpolatorType)
		// .def("set_approx_params", &Interpolator<double>::setInterpolatorParams)
		.def("set_lsdriver", &Interpolator<double>::setLSDriver)
		// .def("set_k", [](Interpolator<double>& self, int k)
		// {
		// 	if (k <= 0)
		// 	{
		// 		throw py::value_error("Invalid value " + std::to_string(k)
		// 		                      + " passed to Interpolator::set_k");
		// 	}
		// 	else
		// 	{
		// 		return self.set_k(k);
		// 	}
		// })
		// .def("set_n", [](Interpolator<double>& self, int n)
		// {
		// 	if (n <= 0)
		// 	{
		// 		throw py::value_error("Invalid value " + std::to_string(n)
		// 		                      + " passed to Interpolator::set_n");
		// 	}
		// 	else
		// 	{
		// 		return self.set_n(n);
		// 	}
		// })
    // .def("set_shape", &Interpolator<double>::set_shape)
		.def("set_flag", &Interpolator<double>::setFlag)
		.def("set_info", &Interpolator<double>::setInfo)
		.def("output", [](Interpolator<double>& self)
		{
			std::string out = self.getLogger()->getOutput();
			py::print(out);
		})
		.def("output", [](Interpolator<double>& self, size_t i)
		{
			std::string out = self.getLogger()->getOutput(i);
			py::print(out);
		})
		// //--------------------------------------------------------------------------
		// //	scalar field methods
		// //--------------------------------------------------------------------------
		// .def("scalar_gradient",
		//      (std::vector<std::vector<double>> (Interpolator<double>::*)
		// 		  (const std::shared_ptr<UGrid<double>>,
		// 		   const std::shared_ptr<ScalarField<double>>))
		// 		 &Interpolator<double>::scalarGradient,
	  //      py::return_value_policy::reference)
	  // .def("scalar_gradient",
		//      (std::vector<std::vector<double>> (Interpolator<double>::*)
		// 		  (const std::shared_ptr<UGrid<double>>,
		// 			 const ScalarField<double>&))
		// 		 &Interpolator<double>::scalarGradient,
	  //      py::return_value_policy::reference)
		//  .def("scalar_gradient_ls",
 		//      (std::vector<std::vector<double>> (Interpolator<double>::*)
 		// 		  (const std::shared_ptr<UGrid<double>>,
 		// 		   const std::shared_ptr<ScalarField<double>>))
 		// 		 &Interpolator<double>::scalarGradientLS,
 	  //      py::return_value_policy::reference)
 	  // .def("scalar_gradient_ls",
 		//      (std::vector<std::vector<double>> (Interpolator<double>::*)
 		// 		  (const std::shared_ptr<UGrid<double>>,
 		// 			 const ScalarField<double>&))
 		// 		 &Interpolator<double>::scalarGradientLS,
 	  //      py::return_value_policy::reference)
	  //--------------------------------------------------------------------------
		//  Scalar derivative (order n)
		//--------------------------------------------------------------------------
	  // .def("scalar_derivative",
		// 		 (std::vector<std::vector<double>> (Interpolator<double>::*)
		// 		 	(const std::shared_ptr<UGrid<double>>,
		//  			 const std::shared_ptr<ScalarField<double>>,
		// 			 size_t))
		//  		 &Interpolator<double>::scalarDerivative,
	  //      py::return_value_policy::reference)
		// //--------------------------------------------------------------------------
	 	// //  Scalar derivative (direction dir, order n)
	 	// //--------------------------------------------------------------------------
		// .def("scalar_derivative",
		// 		 (std::vector<double> (Interpolator<double>::*)
		// 		 	(const std::shared_ptr<UGrid<double>>,
		// 			 const std::shared_ptr<ScalarField<double>>,
		// 			 size_t, size_t))
		// 		 &Interpolator<double>::scalarDerivative,
		// 		 py::return_value_policy::reference)
		// //--------------------------------------------------------------------------
 	 	// //  Scalar derivative (direction dir, order n)
 	 	// //--------------------------------------------------------------------------
 		// .def("scalar_derivative",
 		// 		 (std::vector<double> (Interpolator<double>::*)
 		// 		 	(const std::shared_ptr<UGrid<double>>,
 		// 			 const std::shared_ptr<ScalarField<double>>,
 		// 			 std::vector<size_t>))
 		// 		 &Interpolator<double>::scalarDerivative,
 		// 		 py::return_value_policy::reference)
		// //--------------------------------------------------------------------------
 	 	// //  Scalar derivative (index i, order n)
 	 	// //--------------------------------------------------------------------------
		// .def("scalar_derivative_point",
		// 		 (std::vector<double> (Interpolator<double>::*)
		// 		 	(const std::shared_ptr<UGrid<double>>,
		// 			 const std::shared_ptr<ScalarField<double>>,
		// 			 size_t, size_t))
		// 		 &Interpolator<double>::scalarDerivativePoint,
		// 		 py::return_value_policy::reference)
    // //--------------------------------------------------------------------------
	 	// //  Scalar derivative (point p, order n)
	 	// //--------------------------------------------------------------------------
		// .def("scalar_derivative_point",
		// 		 (std::vector<double> (Interpolator<double>::*)
		// 		 	(const std::shared_ptr<UGrid<double>>,
		// 			 const std::shared_ptr<ScalarField<double>>,
		// 			 std::vector<double>, size_t))
		// 		 &Interpolator<double>::scalarDerivativePoint,
		// 		 py::return_value_policy::reference)
		// //--------------------------------------------------------------------------
 	 	// //  Scalar derivative (index i, direction dir, order n)
 	 	// //--------------------------------------------------------------------------
	  // .def("scalar_derivative_point",
	 	// 		 (double (Interpolator<double>::*)
		// 		 	(const std::shared_ptr<UGrid<double>>,
		// 			 const std::shared_ptr<ScalarField<double>>,
		// 			 size_t, size_t, size_t))
		// 		 &Interpolator<double>::scalarDerivativePoint,
		// 		 py::return_value_policy::reference)
    // //--------------------------------------------------------------------------
	 	// //  Scalar derivative (point p, direction dir, order n)
	 	// //--------------------------------------------------------------------------
	  // .def("scalar_derivative_point",
	 	// 		 (double (Interpolator<double>::*)
		// 		 	(const std::shared_ptr<UGrid<double>>,
		// 			 const std::shared_ptr<ScalarField<double>>,
		// 			 std::vector<double>, size_t, size_t))
		// 		 &Interpolator<double>::scalarDerivativePoint,
		// 		 py::return_value_policy::reference)
		// //--------------------------------------------------------------------------
  	// //  Scalar derivative (index i, direction dir, order n)
  	// //--------------------------------------------------------------------------
 	  // .def("scalar_derivative_point",
 	 	// 		 (double (Interpolator<double>::*)
 		// 		 	(const std::shared_ptr<UGrid<double>>,
 		// 			 const std::shared_ptr<ScalarField<double>>,
 		// 			 size_t, std::vector<size_t>))
 		// 		 &Interpolator<double>::scalarDerivativePoint,
 		// 		 py::return_value_policy::reference)
    // //--------------------------------------------------------------------------
   	// //  Scalar derivative (point p, direction dir, order n)
   	// //--------------------------------------------------------------------------
	  // .def("scalar_derivative_point",
	 	// 		 (double (Interpolator<double>::*)
		// 		 	(const std::shared_ptr<UGrid<double>>,
		// 			 const std::shared_ptr<ScalarField<double>>,
		// 			 std::vector<double>, std::vector<size_t>))
		// 		 &Interpolator<double>::scalarDerivativePoint,
		// 		 py::return_value_policy::reference)
		// //--------------------------------------------------------------------------
		// //  Passing Scalarfield by reference
		// //--------------------------------------------------------------------------
		// //--------------------------------------------------------------------------
 	 	// //  Scalar derivative (order n)
 	 	// //--------------------------------------------------------------------------
		// .def("scalar_derivative",
 		// 		 (std::vector<std::vector<double>> (Interpolator<double>::*)
 		// 		 	(const std::shared_ptr<UGrid<double>>,
 		//  			 const ScalarField<double>&,
 		// 			 size_t))
 		//  		 &Interpolator<double>::scalarDerivative,
 	  //      py::return_value_policy::reference)
		// //--------------------------------------------------------------------------
 	 	// //  Scalar derivative (direction dir, order n)
 	 	// //--------------------------------------------------------------------------
 	  // .def("scalar_derivative",
 		// 		 (std::vector<double> (Interpolator<double>::*)
 		// 		 	(const std::shared_ptr<UGrid<double>>,
 		// 			 const ScalarField<double>&,
 		// 			 size_t, size_t))
 		// 		 &Interpolator<double>::scalarDerivative,
 		// 		 py::return_value_policy::reference)
	  // //--------------------------------------------------------------------------
	 	// //  Scalar derivative (direction dir, order n)
	 	// //--------------------------------------------------------------------------
	  // .def("scalar_derivative",
		// 		 (std::vector<double> (Interpolator<double>::*)
		// 		 	(const std::shared_ptr<UGrid<double>>,
		// 			 const ScalarField<double>&,
		// 			 std::vector<size_t>))
		// 		 &Interpolator<double>::scalarDerivative,
		// 		 py::return_value_policy::reference)
		// //--------------------------------------------------------------------------
		// //  Scalar derivative (index i, order n)
		// //--------------------------------------------------------------------------
 		// .def("scalar_derivative_point",
 		// 		 (std::vector<double> (Interpolator<double>::*)
 		// 		 	(const std::shared_ptr<UGrid<double>>,
 		// 			 const ScalarField<double>&,
 		// 			 size_t, size_t))
 		// 		 &Interpolator<double>::scalarDerivativePoint,
 		// 		 py::return_value_policy::reference)
    // //--------------------------------------------------------------------------
		// //  Scalar derivative (point p, order n)
		// //--------------------------------------------------------------------------
		// .def("scalar_derivative_point",
		// 		 (std::vector<double> (Interpolator<double>::*)
		// 		 	(const std::shared_ptr<UGrid<double>>,
		// 			 const ScalarField<double>&,
		// 			 std::vector<double>, size_t))
		// 		 &Interpolator<double>::scalarDerivativePoint,
		// 		 py::return_value_policy::reference)
		// //--------------------------------------------------------------------------
		// //  Scalar derivative (index i, direction dir, order n)
		// //--------------------------------------------------------------------------
	  // .def("scalar_derivative_point",
	 	// 		 (double (Interpolator<double>::*)
		// 		 	(const std::shared_ptr<UGrid<double>>,
		// 			 const ScalarField<double>&,
		// 			 size_t, size_t, size_t))
		// 		 &Interpolator<double>::scalarDerivativePoint,
		// 		 py::return_value_policy::reference)
    // //--------------------------------------------------------------------------
 		// //  Scalar derivative (point p, direction dir, order n)
 		// //--------------------------------------------------------------------------
 	  // .def("scalar_derivative_point",
 	 	// 		 (double (Interpolator<double>::*)
 		// 		 	(const std::shared_ptr<UGrid<double>>,
 		// 			 const ScalarField<double>&,
 		// 			 std::vector<double>, size_t, size_t))
 		// 		 &Interpolator<double>::scalarDerivativePoint,
 		// 		 py::return_value_policy::reference)
		// //--------------------------------------------------------------------------
 		// //  Scalar derivative (index i, direction dir, order n)
 		// //--------------------------------------------------------------------------
 	  // .def("scalar_derivative_point",
 	 	// 		 (double (Interpolator<double>::*)
 		// 		 	(const std::shared_ptr<UGrid<double>>,
 		// 			 const ScalarField<double>&,
 		// 			 size_t, std::vector<size_t>))
 		// 		 &Interpolator<double>::scalarDerivativePoint,
 		// 		 py::return_value_policy::reference)
    // //--------------------------------------------------------------------------
		// //  Scalar derivative (point p, direction dir, order n)
		// //--------------------------------------------------------------------------
	  // .def("scalar_derivative_point",
	 	// 		 (double (Interpolator<double>::*)
		// 		 	(const std::shared_ptr<UGrid<double>>,
		// 			 const ScalarField<double>&,
		// 			 std::vector<double>, std::vector<size_t>))
		// 		 &Interpolator<double>::scalarDerivativePoint,
		// 		 py::return_value_policy::reference)
    // //--------------------------------------------------------------------------
    // //  Various matrices
    // //--------------------------------------------------------------------------
		// .def("construct_taylor_matrix",
		//      (Matrix<double> (Interpolator<double>::*)
		// 		  (const std::shared_ptr<UGrid<double>>,
    //        const std::vector<size_t>,size_t,size_t))
		// 		 &Interpolator<double>::constructTaylorMatrix)
    // .def("construct_taylor_matrix",
    //      (Matrix<double> (Interpolator<double>::*)
    // 		  (const std::shared_ptr<UGrid<double>>,
    //       const std::vector<size_t>,std::vector<double>,size_t))
    // 		 &Interpolator<double>::constructTaylorMatrix)
		// .def("construct_taylor_matrix",
		//      (Matrix<double> (Interpolator<double>::*)
		// 		  (const std::shared_ptr<UGrid<double>>,
    //        const std::vector<size_t>,size_t,Monomial&))
		// 		 &Interpolator<double>::constructTaylorMatrix)
		//--------------------------------------------------------------------------
		//	print functionality
		//--------------------------------------------------------------------------
    .def("__repr__", [](Interpolator<double> &field)
		{
			std::stringstream s;
			s << &field;
			std::string res = "++++++++++++++++++++++++++++++++++++++++++++++++++++";
			res += "\n<etraj.Interpolator ref at " + s.str() + ">\n";
			return res + field.summary();
		})
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

  //----------------------------------------------------------------------------
	//	DiffEQ Base class
	//----------------------------------------------------------------------------
  py::class_<DiffEQ<double>,
             std::shared_ptr<DiffEQ<double>>>(m, "DiffEQ")
    .def(py::init<>())
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

  //----------------------------------------------------------------------------
	//	KDTree class
	//----------------------------------------------------------------------------
  py::class_<KDTree<double>,
             std::shared_ptr<KDTree<double>>>(m, "KDTree")
    .def(py::init<>())
    ;
  //----------------------------------------------------------------------------

}
