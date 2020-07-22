//------------------------------------------------------------------------------
//  Interpolator.cpp
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
#include "interpolator.h"


namespace ET
{
  //  map for a string to enum of Interpolator type
  std::map<std::string, InterpolatorType> InterpolatorTypeMap =
  {
    { "DEFAULT", InterpolatorType::DEFAULT },
    { "LTE",     InterpolatorType::LTE },
    { "RBF",     InterpolatorType::RBF },
  };
  //  map for a  enum to string of Interpolator type
  std::map<InterpolatorType, std::string> InterpolatorTypeNameMap =
  {
    { InterpolatorType::DEFAULT, "DEFAULT" },
    { InterpolatorType::LTE,     "LTE" },
    { InterpolatorType::RBF,     "RBF" },
  };
  //  map for a string to enum of LSDriver type
  std::map<std::string, LSDriver> LSDriverMap =
  {
    { "xGELS",  LSDriver::xGELS },
    { "xGELSY", LSDriver::xGELSY },
    { "xGELSD", LSDriver::xGELSD },
    { "xGELSS", LSDriver::xGELSS }
  };
  //  map for a string to enum of LSDriver type
  std::map<LSDriver, std::string> LSDriverNameMap =
  {
    { LSDriver::xGELS,  "xGELS" },
    { LSDriver::xGELSY, "xGELSY" },
    { LSDriver::xGELSD, "xGELSD" },
    { LSDriver::xGELSS, "xGELSS" }
  };
  //  map for a string to enum of Solver type
  std::map<std::string, SolverType> SolverTypeMap =
  {
    { "LS",   SolverType::LS },
    { "MLS",  SolverType::MLS },
    { "WMLS", SolverType::WMLS },
  };
  //  map for a  enum to string of Solver type
  std::map<SolverType, std::string> SolverTypeNameMap =
  {
    { SolverType::LS,   "LS" },
    { SolverType::MLS,  "MLS" },
    { SolverType::WMLS, "WMLS" },
  };
  //  map for a string to enum of WeightMatrix type
  std::map<std::string, WeightMatrixType> WeightMatrixTypeMap =
  {
    { "GAUSSIAN",   WeightMatrixType::GAUSSIAN },
  };
  //  map for a  enum to string of WeightMatrix type
  std::map<WeightMatrixType, std::string> WeightMatrixTypeNameMap =
  {
    { WeightMatrixType::GAUSSIAN,   "GAUSSIAN" },
  };
  //----------------------------------------------------------------------------
	template<typename T>
  Interpolator<T>::Interpolator()
  : m_name("default")
  {
    m_lsdriver = LSDriver::xGELS;
    m_Grid = std::make_shared<Grid<T>>();
		m_log = std::make_shared<Log>();
		m_log->init("ET:Interpolator:" + m_name, ".logs/interpolator_" + m_name + ".txt");
		m_log->TRACE("Interpolator '" + m_name + "' created at location "
		            + address_to_string(*this));
  }
	//----------------------------------------------------------------------------
  template<typename T>
  Interpolator<T>::~Interpolator()
  {
		m_log->TRACE("Interpolator '" + m_name
								+ "' destroyed at location " + address_to_string(*this));
	}
  //----------------------------------------------------------------------------
  template<typename T>
  Interpolator<T>::Interpolator(std::string t_name)
  : m_name(t_name)
  {
    m_lsdriver = LSDriver::xGELS;
    m_Grid = std::make_shared<Grid<T>>();
		m_log = std::make_shared<Log>();
		m_log->init("ET:Interpolator:" + m_name, ".logs/interpolator_" + m_name + ".txt");
		m_log->TRACE("Interpolator '" + m_name + "' created at location "
		            + address_to_string(*this));
  }
	//----------------------------------------------------------------------------
  template<typename T>
  Interpolator<T>::Interpolator(std::shared_ptr<Grid<T>> t_Grid)
  : m_name("default"), m_Grid(t_Grid)
  {
    m_lsdriver = LSDriver::xGELS;
		m_log = std::make_shared<Log>();
    m_log->init("ET:Interpolator:" + m_name, ".logs/interpolator_" + m_name + ".txt");
		m_log->TRACE("Interpolator '" + m_name + "' created at location "
		            + address_to_string(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Interpolator<T>::Interpolator(std::string t_name, std::shared_ptr<Grid<T>> t_Grid)
  : m_name(t_name), m_Grid(t_Grid)
  {
    m_lsdriver = LSDriver::xGELS;
		m_log = std::make_shared<Log>();
    m_log->init("ET:Interpolator:" + m_name, ".logs/interpolator_" + m_name + ".txt");
		m_log->TRACE("Interpolator '" + m_name + "' created at location "
		            + address_to_string(*this));
  }
	//----------------------------------------------------------------------------
	template<typename T>
  Interpolator<T>::Interpolator(std::shared_ptr<Log> t_log)
  : m_name("default")
  {
    m_lsdriver = LSDriver::xGELS;
    m_Grid = std::make_shared<Grid<T>>();
		m_log = t_log;
		m_log->TRACE("Interpolator '" + m_name + "' created at location "
		            + address_to_string(*this));
		m_log->INFO("Log passed to Interpolator '" + m_name + "'");
  }
  //----------------------------------------------------------------------------
	template<typename T>
  Interpolator<T>::Interpolator(std::string t_name, std::shared_ptr<Log> t_log)
  : m_name(t_name)
  {
    m_lsdriver = LSDriver::xGELS;
    m_Grid = std::make_shared<Grid<T>>();
		m_log = t_log;
		m_log->TRACE("Interpolator '" + m_name + "' created at location "
		            + address_to_string(*this));
		m_log->INFO("Log passed to Interpolator '" + m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Interpolator<T>::Interpolator(std::shared_ptr<Grid<T>> t_Grid,
                                std::shared_ptr<Log> t_log)
  : m_name("default"), m_Grid(t_Grid), m_log(t_log)
  {
    m_lsdriver = LSDriver::xGELS;
    m_log->TRACE("Interpolator '" + m_name + "' created at location "
                + address_to_string(*this));
    m_log->INFO("Log passed to Interpolator '" + m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Interpolator<T>::Interpolator(std::string t_name,
                                std::shared_ptr<Grid<T>> t_Grid,
                                std::shared_ptr<Log> t_log)
  : m_name(t_name), m_Grid(t_Grid), m_log(t_log)
  {
    m_lsdriver = LSDriver::xGELS;
    m_log->TRACE("Interpolator '" + m_name + "' created at location "
                + address_to_string(*this));
    m_log->INFO("Log passed to Interpolator '" + m_name + "'");
  }
	//----------------------------------------------------------------------------
  template<typename T>
  std::string Interpolator<T>::getName() const
  {
    return m_name;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::shared_ptr<Grid<T>> Interpolator<T>::getGrid() const
  {
    return m_Grid;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::shared_ptr<Field<T>> Interpolator<T>::getField() const
  {
    return m_Field;
  }
  //----------------------------------------------------------------------------
  template<typename T>
	std::shared_ptr<Log> Interpolator<T>::getLog() const
	{
		return m_log;
	}
  //----------------------------------------------------------------------------
  template<typename T>
  LSDriver Interpolator<T>::getLSDriver() const
  {
    return m_lsdriver;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  SolverType Interpolator<T>::getSolverType() const
  {
    return m_solvertype;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  int Interpolator<T>::getFlag() const
  {
    return m_flag;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::string Interpolator<T>::getInfo() const
  {
    return m_info;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  InterpolatorType Interpolator<T>::getInterpolatorType() const
  {
    return m_interpolatortype;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Interpolator<T>::setName(std::string t_name)
  {
    m_name = t_name;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Interpolator<T>::setGrid(std::shared_ptr<Grid<T>> t_Grid)
  {
    m_Grid = t_Grid;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Interpolator<T>::setField(std::shared_ptr<Field<T>> t_Field)
  {
    m_Field = t_Field;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Interpolator<T>::setLog(std::shared_ptr<Log> t_log)
  {
    m_log = t_log;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Interpolator<T>::setLSDriver(std::string t_type)
  {
		auto res = LSDriverMap.find(t_type);
		if (res == LSDriverMap.end()) {
			m_log->ERROR("Interpolator " + m_name + ": Attempted to set LSDriver "
		              + "to " + t_type + " which is not a valid type");
		}
		else {
			m_lsdriver = LSDriverMap[t_type];
			m_log->INFO("Interpolator " + m_name + ": LSDriver set to " + t_type);
		}
	}
  //----------------------------------------------------------------------------
  template<typename T>
  void Interpolator<T>::setSolverType(std::string t_type)
  {
		auto res = SolverTypeMap.find(t_type);
		if (res == SolverTypeMap.end()) {
			m_log->ERROR("Interpolator " + m_name + ": Attempted to set SolverType "
		              + "to " + t_type + " which is not a valid type");
		}
		else {
			m_solvertype = SolverTypeMap[t_type];
			m_log->INFO("Interpolator " + m_name + ": SolverType set to " + t_type);
		}
	}
  //----------------------------------------------------------------------------
  template<typename T>
  void Interpolator<T>::setFlag(int t_flag)
  {
    m_flag = t_flag;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Interpolator<T>::setInfo(std::string t_info)
  {
    m_info = t_info;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Vector<T> Interpolator<T>::xLinearSolvex(Matrix<T> t_A, Vector<T> t_x,
                                           size_t t_index)
  {
    if (m_solvertype == SolverType::LS) {
      return xGELSx(t_A,t_x);
    }
    else if (m_solvertype == SolverType::MLS) {
      return xMLSx(t_A,t_x);
    }
    else {
      return xWMLSx(t_A,t_x,t_index);
    }
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Matrix<T> Interpolator<T>::xLinearSolvex(Matrix<T> t_A, Matrix<T> t_X,
                                           size_t t_index)
  {
    if (m_solvertype == SolverType::LS) {
      return xGELSx(t_A,t_X);
    }
    else if (m_solvertype == SolverType::MLS) {
      return xMLSx(t_A,t_X);
    }
    else {
      return xWMLSx(t_A,t_X,t_index);
    }
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Vector<T> Interpolator<T>::xLinearSolvex(Matrix<T> t_A, Vector<T> t_x,
                                           std::vector<T> t_point)
  {
    if (m_solvertype == SolverType::LS) {
      return xGELSx(t_A,t_x);
    }
    else if (m_solvertype == SolverType::MLS) {
      return xMLSx(t_A,t_x);
    }
    else {
      return xWMLSx(t_A,t_x,t_point);
    }
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Matrix<T> Interpolator<T>::xLinearSolvex(Matrix<T> t_A, Matrix<T> t_X,
                                           std::vector<T> t_point)
  {
    if (m_solvertype == SolverType::LS) {
      return xGELSx(t_A,t_X);
    }
    else if (m_solvertype == SolverType::MLS) {
      return xMLSx(t_A,t_X);
    }
    else {
      return xWMLSx(t_A,t_X,t_point);
    }
  }
  //----------------------------------------------------------------------------
  template<typename T>
	Vector<T> Interpolator<T>::xGELSx(Matrix<T> t_A, Vector<T> t_x)
	{
		if (m_lsdriver == LSDriver::xGELS) {
			return DGELS(t_A,t_x);
		}
		else if (m_lsdriver == LSDriver::xGELSY) {
			return DGELSY(t_A,t_x);
		}
		else if (m_lsdriver == LSDriver::xGELSD) {
			return DGELSD(t_A,t_x);
		}
		else {
			return DGELSS(t_A,t_x);
		}
	}
  //----------------------------------------------------------------------------
  template<typename T>
	Matrix<T> Interpolator<T>::xGELSx(Matrix<T> t_A, Matrix<T> t_X)
	{
		if (m_lsdriver == LSDriver::xGELS) {
			return DGELS(t_A,t_X);
		}
		else if (m_lsdriver == LSDriver::xGELSY) {
			return DGELSY(t_A,t_X);
		}
		else if (m_lsdriver == LSDriver::xGELSD) {
			return DGELSD(t_A,t_X);
		}
		else {
			return DGELSS(t_A,t_X);
		}
	}
  //----------------------------------------------------------------------------
  template<typename T>
	Vector<T> Interpolator<T>::xMLSx(Matrix<T> t_A, Vector<T> t_x)
	{
    //  Construct the transpose of A
    Matrix<T> A_T = t_A.transpose();
    //  Construct the product A*A_T
    Matrix<T> AA_T = t_A * A_T;
    //  Find the inverse of A*A_T
    Matrix<T> AA_T_inv = DGETRI(AA_T);
    //  Create the solution
    Vector<T> y = AA_T_inv * t_A * t_x;
    return y;
	}
  //----------------------------------------------------------------------------
  template<typename T>
	Matrix<T> Interpolator<T>::xMLSx(Matrix<T> t_A, Matrix<T> t_X)
	{
    //  Construct the transpose of A
    Matrix<T> A_T = t_A.transpose();
    //  Construct the product A*A_T
    Matrix<T> AA_T = t_A * A_T;
    //  Find the inverse of A*A_T
    Matrix<T> AA_T_inv = DGETRI(AA_T);
    //  Create the solution
    Matrix<T> Y = AA_T_inv * t_A * t_X;
    return Y;
	}
  //----------------------------------------------------------------------------
  template<typename T>
	Vector<T> Interpolator<T>::xWMLSx(Matrix<T> t_A, Vector<T> t_x,
                                    size_t t_index)
	{
    //  Construct the transpose of A
    Matrix<T> A_T = t_A.transpose();
    //  Construct the product A*A_T
    Matrix<T> AA_T = t_A * A_T;
    //  Find the inverse of A*A_T
    Matrix<T> AA_T_inv = DGETRI(AA_T);
    //  Create the solution
    Vector<T> y = AA_T_inv * t_A * t_x;
    return y;
	}
  //----------------------------------------------------------------------------
  template<typename T>
	Matrix<T> Interpolator<T>::xWMLSx(Matrix<T> t_A, Matrix<T> t_X,
                                    size_t t_index)
	{
    //  Construct the transpose of A
    Matrix<T> A_T = t_A.transpose();
    //  Construct the product A*A_T
    Matrix<T> AA_T = t_A * A_T;
    //  Find the inverse of A*A_T
    Matrix<T> AA_T_inv = DGETRI(AA_T);
    //  Create the solution
    Matrix<T> Y = AA_T_inv * t_A * t_X;
    return Y;
	}
  //----------------------------------------------------------------------------
  template<typename T>
	Vector<T> Interpolator<T>::xWMLSx(Matrix<T> t_A, Vector<T> t_x,
                                    std::vector<T> t_point)
	{
    //  Construct the transpose of A
    Matrix<T> A_T = t_A.transpose();
    //  Construct the product A*A_T
    Matrix<T> AA_T = t_A * A_T;
    //  Find the inverse of A*A_T
    Matrix<T> AA_T_inv = DGETRI(AA_T);
    //  Create the solution
    Vector<T> y = AA_T_inv * t_A * t_x;
    return y;
	}
  //----------------------------------------------------------------------------
  template<typename T>
	Matrix<T> Interpolator<T>::xWMLSx(Matrix<T> t_A, Matrix<T> t_X,
                                    std::vector<T> t_point)
	{
    //  Construct the transpose of A
    Matrix<T> A_T = t_A.transpose();
    //  Construct the product A*A_T
    Matrix<T> AA_T = t_A * A_T;
    //  Find the inverse of A*A_T
    Matrix<T> AA_T_inv = DGETRI(AA_T);
    //  Create the solution
    Matrix<T> Y = AA_T_inv * t_A * t_X;
    return Y;
	}
  //----------------------------------------------------------------------------
  template<typename T>
  Vector<T> Interpolator<T>::derivative(const size_t t_index,
                                        const size_t t_degree)
  {
   return Vector<T>();
  }
  //----------------------------------------------------------------------------
  template<typename T>
  T Interpolator<T>::derivative(const size_t t_index,
                                const size_t t_degree,
                                const size_t t_direction)
  {
   return 0;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Vector<T> Interpolator<T>::derivative(const std::vector<T>& point,
                                        const size_t t_degree)
  {
    return Vector<T>();
  }
  //----------------------------------------------------------------------------
  template<typename T>
  T Interpolator<T>::derivative(const std::vector<T>& point,
                                const size_t t_degree,
                                const size_t t_direction)
  {
    return 0;
  }
  //----------------------------------------------------------------------------

	// //----------------------------------------------------------------------------
	// //  nth-derivatives of scalar Field
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //  scalarDerivative       - approximate the derivative for an entire Field
	// //  Arguments:  Grid      - Grid<T> pointer
	// //              Field      - ScalarField<T> pointer
	// //              n          - order of derivative
	// //
	// //  Returns:    std::vector<std::vector<T>> of the derivatives
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<std::vector<T>>
	// Interpolator<T>::scalarDerivative(const std::shared_ptr<Grid<T>> Grid,
	// 												        	const std::shared_ptr<ScalarField<T>> Field,
	// 													        uint32_t n)
	// {
	// 	std::vector<std::vector<T>> result(Field->getN());
	// 	Monomial mono(Grid->getDim(),n);
	// 	std::vector<Vector<T>> vecs = xScalarDerivative(Grid,Field,n);
	// 	for (uint64_t i = 0; i < Field->getN(); i++)
	// 	{
	// 		result[i] = vecs[i].getVec();
	// 	}
	// 	return result;
	// }
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //  scalarDerivative       - approximate the derivative for an entire Field
	// //                           of order n in the direction dir
	// //  Arguments:  Grid      - Grid<T> pointer
	// //              Field      - ScalarField<T> pointer
	// //              dir        - direction of the derivative
	// //              n          - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the derivative along dir.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<T>
	// Interpolator<T>::scalarDerivative(const std::shared_ptr<Grid<T>> Grid,
	// 												        	const std::shared_ptr<ScalarField<T>> Field,
	// 													        uint32_t dir, uint32_t n)
	// {
	// 	std::vector<T> result(Field->getN());
	// 	std::vector<uint32_t> deriv(Field->getDim(),0);
	// 	deriv[dir] = n;
	// 	Monomial mono(Grid->getDim(),n);
	// 	uint32_t index = mono.getTaylorIndex(deriv);
	// 	std::vector<Vector<T>> vecs = xScalarDerivative(Grid,Field,n);
	// 	for(uint64_t i = 0; i < Field->getN(); i++)
	// 	{
	// 		result[i] = vecs[i](index);
	// 	}
	// 	return result;
	// }
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //  scalarDerivative       - approximate the derivative for an entire Field
	// //                           of order n in the direction dir
	// //  Arguments:  Grid      - Grid<T> pointer
	// //              Field      - ScalarField<T> pointer
	// //              deriv      - vector of ints denoting direction and order
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<T>
	// Interpolator<T>::scalarDerivative(const std::shared_ptr<Grid<T>> Grid,
	// 												        	const std::shared_ptr<ScalarField<T>> Field,
	// 													        std::vector<uint32_t> deriv)
	// {
	// 	std::vector<T> result(Field->getN());
	// 	uint32_t n = std::accumulate(deriv.begin(),deriv.end(),0);
	// 	Monomial mono(Grid->getDim(),n);
	// 	uint32_t index = mono.getTaylorIndex(deriv);
	// 	std::vector<Vector<T>> vecs = xScalarDerivative(Grid,Field,n);
	// 	for(uint64_t i = 0; i < Field->getN(); i++)
	// 	{
	// 		result[i] = vecs[i](index);
	// 	}
	// 	return result;
	// }
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //  scalarDerivativePoint  - approximate the derivative for a point
	// //                           of order n
	// //  Arguments:  Grid      - Grid<T> pointer
	// //              Field      - ScalarField<T> pointer
	// //              index      - index of the point
	// //              n          - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<T>
	// Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<Grid<T>> Grid,
	// 																	const std::shared_ptr<ScalarField<T>> Field,
	// 																	uint64_t index, uint32_t n)
	// {
	// 	std::vector<T> result(Field->getDim(),0.0);
	// 	Vector<T> coefficients = xScalarDerivativePoint(Grid,Field,index,n);
	// 	//  Grab the derivative determined by deriv
	// 	Monomial mono(Grid->getDim(),n);
	// 	//  Trim result to the first Field.getDim() elements
	// 	for (uint32_t j = 0; j < Field->getDim(); j++)
	// 	{
	// 		std::vector<uint32_t> deriv(Field->getDim(),0);
	// 		deriv[j] = n;
	// 		uint32_t l = mono.getTaylorIndex(deriv);
	// 		result[j] = coefficients(l);
	// 	}
	// 	return result;
	// }
	// //----------------------------------------------------------------------------
  //
  // //----------------------------------------------------------------------------
	// //  scalarDerivativePoint  - approximate the derivative for a point
	// //                           of order n
	// //  Arguments:  Grid      - Grid<T> pointer
	// //              Field      - ScalarField<T> pointer
	// //              point      - std::vector<T> of the point
	// //              n          - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<T>
	// Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<Grid<T>> Grid,
	// 																	const std::shared_ptr<ScalarField<T>> Field,
	// 																	std::vector<T> point, uint32_t n)
	// {
	// 	std::vector<T> result(Field->getDim(),0.0);
	// 	Vector<T> coefficients = xScalarDerivativePoint(Grid,Field,point,n);
	// 	//  Grab the derivative determined by deriv
	// 	Monomial mono(Grid->getDim(),n);
	// 	//  Trim result to the first Field.getDim() elements
	// 	for (uint32_t j = 0; j < Field->getDim(); j++)
	// 	{
	// 		std::vector<uint32_t> deriv(Field->getDim(),0);
	// 		deriv[j] = n;
	// 		uint32_t l = mono.getTaylorIndex(deriv);
	// 		result[j] = coefficients(l);
	// 	}
	// 	return result;
	// }
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //  scalarDerivativePoint  - approximate the derivative for a point
	// //                           of order n in direction dir
	// //  Arguments:  Grid      - Grid<T> pointer
	// //              Field      - ScalarField<T> pointer
	// //              index      - index of the point
	// //              dir        - direction of the derivative
	// //              n          - order of the derivative
	// //
	// //  Returns:    T          - the gradient in direction dir and order n.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// T	Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<Grid<T>> Grid,
	// 																	const std::shared_ptr<ScalarField<T>> Field,
	// 																	uint64_t index, uint32_t dir, uint32_t n)
	// {
	// 	Vector<T> coefficients = xScalarDerivativePoint(Grid,Field,index,n);
	// 	//  Grab the derivative determined by deriv
	// 	Monomial mono(Grid->getDim(),n);
	// 	std::vector<uint32_t> deriv(Field->getDim(),0);
	// 	deriv[dir] = n;
	// 	uint32_t l = mono.getTaylorIndex(deriv);
	// 	return coefficients(l);
	// }
	// //----------------------------------------------------------------------------
  //
  // //----------------------------------------------------------------------------
	// //  scalarDerivativePoint  - approximate the derivative for a point
	// //                           of order n in direction dir
	// //  Arguments:  Grid      - Grid<T> pointer
	// //              Field      - ScalarField<T> pointer
	// //              point      - std::vector<T> of the point
	// //              dir        - direction of the derivative
	// //              n          - order of the derivative
	// //
	// //  Returns:    T          - the gradient in direction dir and order n.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// T	Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<Grid<T>> Grid,
	// 																	const std::shared_ptr<ScalarField<T>> Field,
	// 																	std::vector<T> point, uint32_t dir,
  //                                   uint32_t n)
	// {
	// 	Vector<T> coefficients = xScalarDerivativePoint(Grid,Field,point,n);
	// 	//  Grab the derivative determined by deriv
	// 	Monomial mono(Grid->getDim(),n);
	// 	std::vector<uint32_t> deriv(Field->getDim(),0);
	// 	deriv[dir] = n;
	// 	uint32_t l = mono.getTaylorIndex(deriv);
	// 	return coefficients(l);
	// }
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //  scalarDerivativePoint  - approximate the derivative for a point
	// //                           of order n in direction dir
	// //  Arguments:  Grid      - Grid<T> pointer
	// //              Field      - ScalarField<T> pointer
	// //              index      - index of the point
	// //              deriv      - vector denoting the direction and order
	// //
	// //  Returns:    T          - the gradient in direction dir and order n.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// T	Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<Grid<T>> Grid,
	// 																	const std::shared_ptr<ScalarField<T>> Field,
	// 																	uint64_t index, std::vector<uint32_t> deriv)
	// {
	// 	uint32_t n = std::accumulate(deriv.begin(),deriv.end(),0);
	// 	Vector<T> coefficients = xScalarDerivativePoint(Grid,Field,index,n);
	// 	//  Grab the derivative determined by deriv
	// 	Monomial mono(Grid->getDim(),n);
	// 	uint32_t l = mono.getTaylorIndex(deriv);
	// 	return coefficients(l);
	// }
	// //----------------------------------------------------------------------------
  //
  // //----------------------------------------------------------------------------
	// //  scalarDerivativePoint  - approximate the derivative for a point
	// //                           of order n in direction dir
	// //  Arguments:  Grid      - Grid<T> pointer
	// //              Field      - ScalarField<T> pointer
	// //              point      - std::vector<T> of the point
  // //              deriv      - vector denoting the direction and order
	// //
	// //  Returns:    T          - the gradient in direction dir and order n.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// T	Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<Grid<T>> Grid,
	// 																	const std::shared_ptr<ScalarField<T>> Field,
	// 																	std::vector<T> point,
  //                                   std::vector<uint32_t> deriv)
	// {
	// 	uint32_t n = std::accumulate(deriv.begin(),deriv.end(),0);
	// 	Vector<T> coefficients = xScalarDerivativePoint(Grid,Field,point,n);
	// 	//  Grab the derivative determined by deriv
	// 	Monomial mono(Grid->getDim(),n);
	// 	uint32_t l = mono.getTaylorIndex(deriv);
	// 	return coefficients(l);
	// }
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //  Passing Field as a const reference
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //  scalarDerivative       - approximate the derivative for an entire Field
	// //  Arguments:  Grid      - Grid<T> pointer
	// //              Field      - const ScalarField<T>& reference
	// //              n          - order of derivative
	// //
	// //  Returns:    std::vector<std::vector<T>> of the derivatives
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<std::vector<T>>
	// Interpolator<T>::scalarDerivative(const std::shared_ptr<Grid<T>> Grid,
	// 												        	const ScalarField<T>& Field,
	// 													        uint32_t n)
	// {
	// 	std::vector<std::vector<T>> result(Field.getN());
	// 	Monomial mono(Grid->getDim(),n);
	// 	std::vector<Vector<T>> vecs = xScalarDerivative(Grid,Field,n);
	// 	for (uint64_t i = 0; i < Field.getN(); i++)
	// 	{
	// 		result[i] = vecs[i].getVec();
	// 	}
	// 	return result;
	// }
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //  scalarDerivative       - approximate the derivative for an entire Field
	// //                           of order n in the direction dir
	// //  Arguments:  Grid      - Grid<T> pointer
	// //              Field      - const ScalarField<T>& reference
	// //              dir        - direction of the derivative
	// //              n          - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the derivative along dir.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<T>
	// Interpolator<T>::scalarDerivative(const std::shared_ptr<Grid<T>> Grid,
	// 												        	const ScalarField<T>& Field,
	// 													        uint32_t dir, uint32_t n)
	// {
	// 	std::vector<T> result(Field.getN());
	// 	std::vector<uint32_t> deriv(Field.getDim(),0);
	// 	deriv[dir] = n;
	// 	Monomial mono(Grid->getDim(),n);
	// 	uint32_t index = mono.getTaylorIndex(deriv);
	// 	std::vector<Vector<T>> vecs = xScalarDerivative(Grid,Field,n);
	// 	for(uint64_t i = 0; i < Field.getN(); i++)
	// 	{
	// 		result[i] = vecs[i](index);
	// 	}
	// 	return result;
	// }
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //  scalarDerivative       - approximate the derivative for an entire Field
	// //                           of order n in the direction dir
	// //  Arguments:  Grid      - Grid<T> pointer
	// //              Field      - const ScalarField<T>& reference
	// //              deriv      - vector of ints denoting direction and order
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<T>
	// Interpolator<T>::scalarDerivative(const std::shared_ptr<Grid<T>> Grid,
	// 												        	const ScalarField<T>& Field,
	// 													        std::vector<uint32_t> deriv)
	// {
	// 	std::vector<T> result(Field.getN());
	// 	uint32_t n = std::accumulate(deriv.begin(),deriv.end(),0);
	// 	Monomial mono(Grid->getDim(),n);
	// 	uint32_t index = mono.getTaylorIndex(deriv);
	// 	std::vector<Vector<T>> vecs = xScalarDerivative(Grid,Field,n);
	// 	for(uint64_t i = 0; i < Field.getN(); i++)
	// 	{
	// 		result[i] = vecs[i](index);
	// 	}
	// 	return result;
	// }
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //  scalarDerivativePoint  - approximate the derivative for a point
	// //                           of order n
	// //  Arguments:  Grid      - Grid<T> pointer
	// //              Field      - const ScalarField<T>& reference
	// //              index      - index of the point
	// //              n          - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<T>
	// Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<Grid<T>> Grid,
	// 												        	const ScalarField<T>& Field,
	// 													        uint64_t index, uint32_t n)
	// {
	// 	std::vector<T> result(Field.getDim(),0.0);
	// 	Vector<T> coefficients = xScalarDerivativePoint(Grid,Field,index,n);
	// 	//  Grab the derivative determined by deriv
	// 	Monomial mono(Grid->getDim(),n);
	// 	//  Trim result to the first Field.getDim() elements
	// 	for (uint32_t j = 0; j < Field.getDim(); j++)
	// 	{
	// 		std::vector<uint32_t> deriv(Field.getDim(),0);
	// 		deriv[j] = n;
	// 		uint32_t l = mono.getTaylorIndex(deriv);
	// 		result[j] = coefficients(l);
	// 	}
	// 	return result;
	// }
	// //----------------------------------------------------------------------------
  //
  // //----------------------------------------------------------------------------
	// //  scalarDerivativePoint  - approximate the derivative for a point
	// //                           of order n
	// //  Arguments:  Grid      - Grid<T> pointer
	// //              Field      - const ScalarField<T>& reference
	// //              point      - std::vector<T> of the point
	// //              n          - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<T>
	// Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<Grid<T>> Grid,
	// 												        	const ScalarField<T>& Field,
	// 													        std::vector<T> point, uint32_t n)
	// {
	// 	std::vector<T> result(Field.getDim(),0.0);
	// 	Vector<T> coefficients = xScalarDerivativePoint(Grid,Field,point,n);
	// 	//  Grab the derivative determined by deriv
	// 	Monomial mono(Grid->getDim(),n);
	// 	//  Trim result to the first Field.getDim() elements
	// 	for (uint32_t j = 0; j < Field.getDim(); j++)
	// 	{
	// 		std::vector<uint32_t> deriv(Field.getDim(),0);
	// 		deriv[j] = n;
	// 		uint32_t l = mono.getTaylorIndex(deriv);
	// 		result[j] = coefficients(l);
	// 	}
	// 	return result;
	// }
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //  scalarDerivativePoint  - approximate the derivative for a point
	// //                           of order n in direction dir
	// //  Arguments:  Grid      - Grid<T> pointer
	// //              Field      - const ScalarField<T>& pointer
	// //              index      - index of the point
	// //              dir        - direction of the derivative
	// //              n          - order of the derivative
	// //
	// //  Returns:    T          - the gradient in direction dir and order n.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// T Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<Grid<T>> Grid,
	// 																	  const ScalarField<T>& Field,
	// 																	  uint64_t index, uint32_t dir, uint32_t n)
	// {
	// 	Vector<T> coefficients = xScalarDerivativePoint(Grid,Field,index,n);
	// 	//  Grab the derivative determined by deriv
	// 	Monomial mono(Grid->getDim(),n);
	// 	std::vector<uint32_t> deriv(Field.getDim(),0);
	// 	deriv[dir] = n;
	// 	uint32_t l = mono.getTaylorIndex(deriv);
	// 	return coefficients(l);
	// }
	// //----------------------------------------------------------------------------
  //
  // //----------------------------------------------------------------------------
	// //  scalarDerivativePoint  - approximate the derivative for a point
	// //                           of order n in direction dir
	// //  Arguments:  Grid      - Grid<T> pointer
	// //              Field      - const ScalarField<T>& pointer
	// //              point      - std::vector<T> of the point
	// //              dir        - direction of the derivative
	// //              n          - order of the derivative
	// //
	// //  Returns:    T          - the gradient in direction dir and order n.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// T Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<Grid<T>> Grid,
	// 																	  const ScalarField<T>& Field,
	// 																	  std::vector<T> point, uint32_t dir,
  //                                     uint32_t n)
	// {
	// 	Vector<T> coefficients = xScalarDerivativePoint(Grid,Field,point,n);
	// 	//  Grab the derivative determined by deriv
	// 	Monomial mono(Grid->getDim(),n);
	// 	std::vector<uint32_t> deriv(Field.getDim(),0);
	// 	deriv[dir] = n;
	// 	uint32_t l = mono.getTaylorIndex(deriv);
	// 	return coefficients(l);
	// }
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //  scalarDerivativePoint  - approximate the derivative for a point
	// //                           of order n in direction dir
	// //  Arguments:  Grid      - Grid<T> pointer
	// //              Field      - const ScalarField<T>& reference
	// //              index      - index of the point
	// //              deriv      - vector denoting the direction and order
	// //
	// //  Returns:    T          - the gradient in direction dir and order n.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// T Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<Grid<T>> Grid,
	// 																	  const ScalarField<T>& Field,
	// 																	  uint64_t index,
	// 																		std::vector<uint32_t> deriv)
	// {
	// 	uint32_t n = std::accumulate(deriv.begin(),deriv.end(),0);
	// 	Vector<T> coefficients = xScalarDerivativePoint(Grid,Field,index,n);
	// 	//  Grab the derivative determined by deriv
	// 	Monomial mono(Grid->getDim(),n);
	// 	uint32_t l = mono.getTaylorIndex(deriv);
	// 	return coefficients(l);
	// }
	// //----------------------------------------------------------------------------
  //
  // //----------------------------------------------------------------------------
	// //  scalarDerivativePoint  - approximate the derivative for a point
	// //                           of order n in direction dir
	// //  Arguments:  Grid      - Grid<T> pointer
	// //              Field      - const ScalarField<T>& reference
	// //              point      - std::vector<T> of the point
	// //              deriv      - vector denoting the direction and order
	// //
	// //  Returns:    T          - the gradient in direction dir and order n.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// T Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<Grid<T>> Grid,
	// 																	  const ScalarField<T>& Field,
	// 																	  std::vector<T> point,
	// 																		std::vector<uint32_t> deriv)
	// {
	// 	uint32_t n = std::accumulate(deriv.begin(),deriv.end(),0);
	// 	Vector<T> coefficients = xScalarDerivativePoint(Grid,Field,point,n);
	// 	//  Grab the derivative determined by deriv
	// 	Monomial mono(Grid->getDim(),n);
	// 	uint32_t l = mono.getTaylorIndex(deriv);
	// 	return coefficients(l);
	// }
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //  Driver routines for derivatives
	// //----------------------------------------------------------------------------
	// //----------------------------------------------------------------------------
	// //  Vanilla least squares
	// //----------------------------------------------------------------------------
	// //----------------------------------------------------------------------------
	// //  scalarDerivativeLSPoint - approximate the derivative for a point
	// //                             using the vanilla LS method
	// //  Arguments:  Grid        - Grid<T> pointer
	// //              Field        - ScalarField<T> pointer
	// //              index        - uint64_t
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// Vector<T>
	// Interpolator<T>::scalarDerivativeLSPoint(const std::shared_ptr<Grid<T>> Grid,
	// 												 const std::shared_ptr<ScalarField<T>> Field,
	// 												 uint64_t index,
	// 												 uint32_t n)
	// {
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(Grid->getDim(),n);
	// 	//	Query the nearest neighbors of Grid
	// 	Grid->queryNeighbors(_params.k);
	// 	//	Get the nearest neighbors for the point 'index'
	// 	std::vector<uint64_t> neighbors = Grid->getNeighbors(index);
	// 	//	Construct the taylor matrix B
	// 	Matrix<T> B = constructTaylorMatrix(Grid,neighbors,index,mono);
	// 	//	Construct the vector of corresponding Field values for neighboPointrs
	// 	Vector<T> Field_neighbors(_params.k,0.0);
	// 	for (uint32_t i = 0; i < _params.k; i++)
	// 	{
	// 		Field_neighbors(i) = (*Field)(neighbors[i]);
	// 	}
	// 	//	Complete the vanilla least squares to get the coefficients
	// 	Vector<T> coefficients = xGELSx(B,Field_neighbors);
	// 	if (B.getFlag() == -1)
	// 	{
	// 		m_log->ERROR(B.getInfo());
	// 	}
	// 	return coefficients;
	// }
	// //----------------------------------------------------------------------------
	// //  scalarDerivativeLSPoint - approximate the derivative for a point
	// //                             using the vanilla LS method
	// //  Arguments:  Grid        - Grid<T> pointer
	// //              Field        - ScalarField<T> pointer
	// //              point        - std::vector<T> of the point
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// Vector<T>
	// Interpolator<T>::scalarDerivativeLSPoint(const std::shared_ptr<Grid<T>> Grid,
	// 												 const std::shared_ptr<ScalarField<T>> Field,
	// 												 std::vector<T> point,
	// 												 uint32_t n)
	// {
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(Grid->getDim(),n);
	// 	//	Get the nearest neighbors for the point 'index'
	// 	std::vector<uint64_t> neighbors = Grid->queryNeighbors(point,_params.k);
	// 	//	Construct the taylor matrix B
	// 	Matrix<T> B = constructTaylorMatrix(Grid,neighbors,point,n);
	// 	//	Construct the vector of corresponding Field values for neighboPointrs
	// 	Vector<T> Field_neighbors(_params.k,0.0);
	// 	for (uint32_t i = 0; i < _params.k; i++)
	// 	{
	// 		Field_neighbors(i) = (*Field)(neighbors[i]);
	// 	}
	// 	//	Complete the vanilla least squares to get the coefficients
	// 	Vector<T> coefficients = xGELSx(B,Field_neighbors);
	// 	if (B.getFlag() == -1)
	// 	{
	// 		m_log->ERROR(B.getInfo());
	// 	}
	// 	return coefficients;
	// }
	// //----------------------------------------------------------------------------
	// //  scalarDerivativeLS      - approximate the derivative for an entire
	// //                             Field using the vanilla LS method
	// //  Arguments:  Grid        - Grid<T> pointer
	// //              Field        - ScalarField<T> pointer
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<Vector<T>>
	// Interpolator<T>::scalarDerivativeLS(const std::shared_ptr<Grid<T>> Grid,
	// 										const std::shared_ptr<ScalarField<T>> Field,
	// 										uint32_t n)
	// {
	// 	std::vector<Vector<T>> result(Field->getN());
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(Grid->getDim(),n);
	// 	//	Query the nearest neighbors of Grid
	// 	Grid->queryNeighbors(_params.k);
	// 	for (uint64_t i = 0; i < Field->getN(); i++)
	// 	{
	// 		//	Get the nearest neighbors for the point 'index'
	// 		std::vector<uint64_t> neighbors = Grid->getNeighbors(i);
	// 		//	Construct the taylor matrix B
	// 		Matrix<T> B = constructTaylorMatrix(Grid,neighbors,i,mono);
	// 		//	Construct the vector of corresponding Field values for neighbors
	// 		Vector<T> Field_neighbors(_params.k,0.0);
	// 		for (uint32_t i = 0; i < _params.k; i++)
	// 		{
	// 			Field_neighbors(i) = (*Field)(neighbors[i]);
	// 		}
	// 		//	Complete the vanilla least squares to get the coefficients
	// 		Vector<T> coefficients = xGELSx(B,Field_neighbors);
	// 		if (B.getFlag() == -1)
	// 		{
	// 			m_log->ERROR(B.getInfo());
	// 		}
	// 		result[i] = coefficients;
	// 	}
	// 	return result;
	// }
	// //----------------------------------------------------------------------------
	// //----------------------------------------------------------------------------
	// //  Passing Field as a const reference
	// //----------------------------------------------------------------------------
	// //----------------------------------------------------------------------------
	// //  scalarDerivativeLSPoint - approximate the derivative for a point
	// //                             using the vanilla LS method
	// //  Arguments:  Grid        - Grid<T> pointer
	// //              Field        - ScalarField<T> pointer
	// //              index        - uint64_t
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// Vector<T>
	// Interpolator<T>::scalarDerivativeLSPoint(const std::shared_ptr<Grid<T>> Grid,
	// 												 const ScalarField<T>& Field,
	// 												 uint64_t index,
	// 												 uint32_t n)
	// {
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(Grid->getDim(),n);
	// 	//	Query the nearest neighbors of Grid
	// 	Grid->queryNeighbors(_params.k);
	// 	//	Get the nearest neighbors for the point 'index'
	// 	std::vector<uint64_t> neighbors = Grid->getNeighbors(index);
	// 	//	Construct the taylor matrix B
	// 	Matrix<T> B = constructTaylorMatrix(Grid,neighbors,index,mono);
	// 	//	Construct the vector of corresponding Field values for neighbors
	// 	Vector<T> Field_neighbors(_params.k,0.0);
	// 	for (uint32_t i = 0; i < _params.k; i++)
	// 	{
	// 		Field_neighbors(i) = Field(neighbors[i]);
	// 	}
	// 	//	Complete the vanilla least squares to get the coefficients
	// 	Vector<T> coefficients = xGELSx(B,Field_neighbors);
	// 	if (B.getFlag() == -1)
	// 	{
	// 		m_log->ERROR(B.getInfo());
	// 	}
	// 	return coefficients;
	// }
	// //----------------------------------------------------------------------------
	// //  scalarDerivativeLSPoint - approximate the derivative for a point
	// //                             using the vanilla LS method
	// //  Arguments:  Grid        - Grid<T> pointer
	// //              Field        - const ScalarField<T>& reference
	// //              point        - std::vector<T> of the point
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// Vector<T>
	// Interpolator<T>::scalarDerivativeLSPoint(const std::shared_ptr<Grid<T>> Grid,
	// 												 const ScalarField<T>& Field,
	// 												 std::vector<T> point,
	// 												 uint32_t n)
	// {
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(Grid->getDim(),n);
	// 	//	Get the nearest neighbors for the point 'index'
	// 	std::vector<uint64_t> neighbors = Grid->queryNeighbors(point,_params.k);
	// 	//	Construct the taylor matrix B
	// 	Matrix<T> B = constructTaylorMatrix(Grid,neighbors,point,n);
	// 	//	Construct the vector of corresponding Field values for neighboPointrs
	// 	Vector<T> Field_neighbors(_params.k,0.0);
	// 	for (uint32_t i = 0; i < _params.k; i++)
	// 	{
	// 		Field_neighbors(i) = Field(neighbors[i]);
	// 	}
	// 	//	Complete the vanilla least squares to get the coefficients
	// 	Vector<T> coefficients = xGELSx(B,Field_neighbors);
	// 	if (B.getFlag() == -1)
	// 	{
	// 		m_log->ERROR(B.getInfo());
	// 	}
	// 	return coefficients;
	// }
	// //----------------------------------------------------------------------------
	// //  scalarDerivativeLS      - approximate the derivative for an entire
	// //                             Field using the vanilla LS method
	// //  Arguments:  Grid        - Grid<T> pointer
	// //              Field        - ScalarField<T> pointer
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<Vector<T>>
	// Interpolator<T>::scalarDerivativeLS(const std::shared_ptr<Grid<T>> Grid,
	// 										const ScalarField<T>& Field,
	// 										uint32_t n)
	// {
	// 	std::vector<Vector<T>> result(Field.getN());
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(Grid->getDim(),n);
	// 	//	Query the nearest neighbors of Grid
	// 	Grid->queryNeighbors(_params.k);
	// 	for (uint64_t i = 0; i < Field.getN(); i++)
	// 	{
	// 		//	Get the nearest neighbors for the point 'index'
	// 		std::vector<uint64_t> neighbors = Grid->getNeighbors(i);
	// 		//	Construct the taylor matrix B
	// 		Matrix<T> B = constructTaylorMatrix(Grid,neighbors,i,mono);
	// 		//	Construct the vector of corresponding Field values for neighbors
	// 		Vector<T> Field_neighbors(_params.k,0.0);
	// 		for (uint32_t i = 0; i < _params.k; i++)
	// 		{
	// 			Field_neighbors(i) = Field(neighbors[i]);
	// 		}
	// 		//	Complete the vanilla least squares to get the coefficients
	// 		Vector<T> coefficients = xGELSx(B,Field_neighbors);
	// 		if (B.getFlag() == -1)
	// 		{
	// 			m_log->ERROR(B.getInfo());
	// 		}
	// 		result[i] = coefficients;
	// 	}
	// 	return result;
	// }
	// //----------------------------------------------------------------------------
	// //----------------------------------------------------------------------------
	// //  Moving least squares
	// //----------------------------------------------------------------------------
	// //----------------------------------------------------------------------------
	// //  scalarDerivativeMLSPoint - approximate the derivative for a point
	// //                             using the MLS method
	// //  Arguments:  Grid        - Grid<T> pointer
	// //              Field        - ScalarField<T> pointer
	// //              index        - uint64_t
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// Vector<T>
	// Interpolator<T>::scalarDerivativeMLSPoint(const std::shared_ptr<Grid<T>> Grid,
	// 												 const std::shared_ptr<ScalarField<T>> Field,
	// 												 uint64_t index,
	// 												 uint32_t n)
	// {
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(Grid->getDim(),n);
	// 	//	Query the nearest neighbors of Grid
	// 	Grid->queryNeighbors(_params.k);
	// 	//	Get the nearest neighbors for the point 'index'
	// 	std::vector<uint64_t> neighbors = Grid->getNeighbors(index);
	// 	//	Construct the taylor matrix B
	// 	Matrix<T> B = constructTaylorMatrix(Grid,neighbors,index,mono);
	// 	//	Construct the vector of corresponding Field values for neighboPointrs
	// 	Vector<T> Field_neighbors(_params.k,0.0);
	// 	for (uint32_t i = 0; i < _params.k; i++)
	// 	{
	// 		Field_neighbors(i) = (*Field)(neighbors[i]);
	// 	}
	// 	//	Complete the vanilla least squares to get the coefficients
	// 	//	compute the transpose of B
	// 	Matrix<T> B_T = B.transpose();
	// 	//	compute the product (B^T*B)
	// 	Matrix<T> B_TB = B_T * B;
	// 	//	compute the inverse of the product
	// 	Matrix<T> B_TB_inv = DGETRI(B_TB);
	// 	//	finally take the product
	// 	Vector<T> coefficients = B_TB_inv * B_T * Field_neighbors;
	// 	if (B.getFlag() == -1)
	// 	{
	// 		m_log->ERROR(B.getInfo());
	// 		m_log->ERROR(B_TB.getInfo());
	// 	}
	// 	return coefficients;
	// }
	// //----------------------------------------------------------------------------
	// //  scalarDerivativeMLSPoint - approximate the derivative for a point
	// //                             using the MLS method
	// //  Arguments:  Grid        - Grid<T> pointer
	// //              Field        - ScalarField<T> pointer
	// //              point        - coordinates of the point of interest
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// Vector<T>
	// Interpolator<T>::scalarDerivativeMLSPoint(const std::shared_ptr<Grid<T>> Grid,
	// 												 const std::shared_ptr<ScalarField<T>> Field,
	// 												 std::vector<T> point,
	// 												 uint32_t n)
	// {
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(Grid->getDim(),n);
	// 	//	Get the nearest neighbors for the point 'index'
	// 	std::vector<uint64_t> neighbors = Grid->queryNeighbors(point,_params.k);
	// 	//	Construct the taylor matrix B
	// 	Matrix<T> B = constructTaylorMatrix(Grid,neighbors,point,n);
	// 	//	Construct the vector of corresponding Field values for neighboPointrs
	// 	Vector<T> Field_neighbors(_params.k,0.0);
	// 	for (uint32_t i = 0; i < _params.k; i++)
	// 	{
	// 		Field_neighbors(i) = (*Field)(neighbors[i]);
	// 	}
	// 	//	Complete the vanilla least squares to get the coefficients
	// 	//	compute the transpose of B
	// 	Matrix<T> B_T = B.transpose();
	// 	//	compute the product (B^T*B)
	// 	Matrix<T> B_TB = B_T * B;
	// 	//	compute the inverse of the product
	// 	Matrix<T> B_TB_inv = DGETRI(B_TB);
	// 	//	finally take the product
	// 	Vector<T> coefficients = B_TB_inv * B_T * Field_neighbors;
	// 	if (B.getFlag() == -1)
	// 	{
	// 		m_log->ERROR(B.getInfo());
	// 		m_log->ERROR(B_TB.getInfo());
	// 	}
	// 	return coefficients;
	// }
	// //----------------------------------------------------------------------------
	// //  scalarDerivativeMLS      - approximate the derivative for an entire
	// //                             Field using the MLS method
	// //  Arguments:  Grid        - Grid<T> pointer
	// //              Field        - ScalarField<T> pointer
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<Vector<T>>
	// Interpolator<T>::scalarDerivativeMLS(const std::shared_ptr<Grid<T>> Grid,
	// 										const std::shared_ptr<ScalarField<T>> Field,
	// 										uint32_t n)
	// {
	// 	std::vector<Vector<T>> result(Field->getN());
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(Grid->getDim(),n);
	// 	//	Query the nearest neighbors of Grid
	// 	Grid->queryNeighbors(_params.k);
	// 	for (uint64_t i = 0; i < Field->getN(); i++)
	// 	{
	// 		//	Get the nearest neighbors for the point 'index'
	// 		std::vector<uint64_t> neighbors = Grid->getNeighbors(i);
	// 		//	Construct the taylor matrix B
	// 		Matrix<T> B = constructTaylorMatrix(Grid,neighbors,i,mono);
	// 		//	Construct the vector of corresponding Field values for neighbors
	// 		Vector<T> Field_neighbors(_params.k,0.0);
	// 		for (uint32_t i = 0; i < _params.k; i++)
	// 		{
	// 			Field_neighbors(i) = (*Field)(neighbors[i]);
	// 		}
	// 		//	Complete the vanilla least squares to get the coefficients
	// 		//	Complete the vanilla least squares to get the coefficients
	// 		//	compute the transpose of B
	// 		Matrix<T> B_T = B.transpose();
	// 		//	compute the product (B^T*B)
	// 		Matrix<T> B_TB = B_T * B;
	// 		//	compute the inverse of the product
	// 		Matrix<T> B_TB_inv = DGETRI(B_TB);
	// 		//	finally take the product
	// 		Vector<T> coefficients = B_TB_inv * B_T * Field_neighbors;
	// 		if (B.getFlag() == -1)
	// 		{
	// 			m_log->ERROR(B.getInfo());
	// 			m_log->ERROR(B_TB.getInfo());
	// 		}
	// 		result[i] = coefficients;
	// 	}
	// 	return result;
	// }
	// //----------------------------------------------------------------------------
	// //----------------------------------------------------------------------------
	// //  Passing Field as const reference
	// //----------------------------------------------------------------------------
	// //----------------------------------------------------------------------------
	// //  scalarDerivativeMLSPoint - approximate the derivative for a point
	// //                             using the MLS method
	// //  Arguments:  Grid        - Grid<T> pointer
	// //              Field        - const ScalarField<T>& reference
	// //              index        - uint64_t
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// Vector<T>
	// Interpolator<T>::scalarDerivativeMLSPoint(const std::shared_ptr<Grid<T>> Grid,
	// 												 const ScalarField<T>& Field,
	// 												 uint64_t index,
	// 												 uint32_t n)
	// {
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(Grid->getDim(),n);
	// 	//	Query the nearest neighbors of Grid
	// 	Grid->queryNeighbors(_params.k);
	// 	//	Get the nearest neighbors for the point 'index'
	// 	std::vector<uint64_t> neighbors = Grid->getNeighbors(index);
	// 	//	Construct the taylor matrix B
	// 	Matrix<T> B = constructTaylorMatrix(Grid,neighbors,index,mono);
	// 	//	Construct the vector of corresponding Field values for neighboPointrs
	// 	Vector<T> Field_neighbors(_params.k,0.0);
	// 	for (uint32_t i = 0; i < _params.k; i++)
	// 	{
	// 		Field_neighbors(i) = Field(neighbors[i]);
	// 	}
	// 	//	Complete the vanilla least squares to get the coefficients
	// 	//	compute the transpose of B
	// 	Matrix<T> B_T = B.transpose();
	// 	//	compute the product (B^T*B)
	// 	Matrix<T> B_TB = B_T * B;
	// 	//	compute the inverse of the product
	// 	Matrix<T> B_TB_inv = DGETRI(B_TB);
	// 	//	finally take the product
	// 	Vector<T> coefficients = B_TB_inv * B_T * Field_neighbors;
	// 	if (B.getFlag() == -1)
	// 	{
	// 		m_log->ERROR(B.getInfo());
	// 		m_log->ERROR(B_TB.getInfo());
	// 	}
	// 	return coefficients;
	// }
	// //----------------------------------------------------------------------------
	// //  scalarDerivativeMLSPoint - approximate the derivative for a point
	// //                             using the MLS method
	// //  Arguments:  Grid        - Grid<T> pointer
	// //              Field        - const ScalarField<T>& reference
	// //              point        - coordinates of the point of interest
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// Vector<T>
	// Interpolator<T>::scalarDerivativeMLSPoint(const std::shared_ptr<Grid<T>> Grid,
	// 												 const ScalarField<T>& Field,
	// 												 std::vector<T> point,
	// 												 uint32_t n)
	// {
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(Grid->getDim(),n);
	// 	//	Get the nearest neighbors for the point 'index'
	// 	std::vector<uint64_t> neighbors = Grid->queryNeighbors(point,_params.k);
	// 	//	Construct the taylor matrix B
	// 	Matrix<T> B = constructTaylorMatrix(Grid,neighbors,point,n);
	// 	//	Construct the vector of corresponding Field values for neighboPointrs
	// 	Vector<T> Field_neighbors(_params.k,0.0);
	// 	for (uint32_t i = 0; i < _params.k; i++)
	// 	{
	// 		Field_neighbors(i) = Field(neighbors[i]);
	// 	}
	// 	//	Complete the vanilla least squares to get the coefficients
	// 	//	compute the transpose of B
	// 	Matrix<T> B_T = B.transpose();
	// 	//	compute the product (B^T*B)
	// 	Matrix<T> B_TB = B_T * B;
	// 	//	compute the inverse of the product
	// 	Matrix<T> B_TB_inv = DGETRI(B_TB);
	// 	//	finally take the product
	// 	Vector<T> coefficients = B_TB_inv * B_T * Field_neighbors;
	// 	if (B.getFlag() == -1)
	// 	{
	// 		m_log->ERROR(B.getInfo());
	// 		m_log->ERROR(B_TB.getInfo());
	// 	}
	// 	return coefficients;
	// }
	// //----------------------------------------------------------------------------
	// //  scalarDerivativeLS      - approximate the derivative for an entire
	// //                             Field using the vanilla LS method
	// //  Arguments:  Grid        - Grid<T> pointer
	// //              Field        - const ScalarField<T>& reference
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<Vector<T>>
	// Interpolator<T>::scalarDerivativeMLS(const std::shared_ptr<Grid<T>> Grid,
	// 										const ScalarField<T>& Field,
	// 										uint32_t n)
	// {
	// 	std::vector<Vector<T>> result(Field.getN());
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(Grid->getDim(),n);
	// 	//	Query the nearest neighbors of Grid
	// 	Grid->queryNeighbors(_params.k);
	// 	for (uint64_t i = 0; i < Field.getN(); i++)
	// 	{
	// 		//	Get the nearest neighbors for the point 'index'
	// 		std::vector<uint64_t> neighbors = Grid->getNeighbors(i);
	// 		//	Construct the taylor matrix B
	// 		Matrix<T> B = constructTaylorMatrix(Grid,neighbors,i,mono);
	// 		//	Construct the vector of corresponding Field values for neighbors
	// 		Vector<T> Field_neighbors(_params.k,0.0);
	// 		for (uint32_t i = 0; i < _params.k; i++)
	// 		{
	// 			Field_neighbors(i) = Field(neighbors[i]);
	// 		}
	// 		//	Complete the vanilla least squares to get the coefficients
	// 		//	Complete the vanilla least squares to get the coefficients
	// 		//	compute the transpose of B
	// 		Matrix<T> B_T = B.transpose();
	// 		//	compute the product (B^T*B)
	// 		Matrix<T> B_TB = B_T * B;
	// 		//	compute the inverse of the product
	// 		Matrix<T> B_TB_inv = DGETRI(B_TB);
	// 		//	finally take the product
	// 		Vector<T> coefficients = B_TB_inv * B_T * Field_neighbors;
	// 		if (B.getFlag() == -1)
	// 		{
	// 			m_log->ERROR(B.getInfo());
	// 			m_log->ERROR(B_TB.getInfo());
	// 		}
	// 		result[i] = coefficients;
	// 	}
	// 	return result;
	// }
	// //----------------------------------------------------------------------------
	// //----------------------------------------------------------------------------
	// //  Weighted Moving least squares
	// //----------------------------------------------------------------------------
	// //----------------------------------------------------------------------------
	// //  scalarDerivativeWMLSPoint - approximate the derivative for a point
	// //                             using the WMLS method
	// //  Arguments:  Grid        - Grid<T> pointer
	// //              Field        - ScalarField<T> pointer
	// //              index        - uint64_t
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// Vector<T>
	// Interpolator<T>::scalarDerivativeWMLSPoint(const std::shared_ptr<Grid<T>> Grid,
	// 												 const std::shared_ptr<ScalarField<T>> Field,
	// 												 uint64_t index,
	// 												 uint32_t n)
	// {
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(Grid->getDim(),n);
	// 	//	Query the nearest neighbors of Grid
	// 	Grid->queryNeighbors(_params.k);
	// 	//	Get the nearest neighbors for the point 'index'
	// 	std::vector<uint64_t> neighbors = Grid->getNeighbors(index);
	// 	//	Construct the taylor matrix B
	// 	Matrix<T> B = constructTaylorMatrix(Grid,neighbors,index,mono);
	// 	//	Construct the weight matrix W
	// 	Matrix<T> W = constructWeightMatrix(Grid,neighbors,index);
	// 	//	Construct the vector of corresponding Field values for neighboPointrs
	// 	Vector<T> Field_neighbors(_params.k,0.0);
	// 	for (uint32_t i = 0; i < _params.k; i++)
	// 	{
	// 		Field_neighbors(i) = (*Field)(neighbors[i]);
	// 	}
	// 	//	Complete the vanilla least squares to get the coefficients
	// 	//	compute the transpose of B
	// 	Matrix<T> B_T = B.transpose();
	// 	//	compute the product (B^T*B)
	// 	Matrix<T> B_TWB = B_T * W * B;
	// 	//	compute the inverse of the product
	// 	Matrix<T> B_TWB_inv = DGETRI(B_TWB);
	// 	//	finally take the product
	// 	Vector<T> coefficients = B_TWB_inv * B_T * W * Field_neighbors;
	// 	if (B.getFlag() == -1)
	// 	{
	// 		m_log->ERROR(B.getInfo());
	// 		m_log->ERROR(B_TWB.getInfo());
	// 	}
	// 	return coefficients;
	// }
	// //----------------------------------------------------------------------------
	// //  scalarDerivativeWMLSPoint - approximate the derivative for a point
	// //                             using the WMLS method
	// //  Arguments:  Grid        - Grid<T> pointer
	// //              Field        - ScalarField<T> pointer
	// //              point        - std::vector<T> of the point
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// Vector<T>
	// Interpolator<T>::scalarDerivativeWMLSPoint(const std::shared_ptr<Grid<T>> Grid,
	// 												 const std::shared_ptr<ScalarField<T>> Field,
	// 												 std::vector<T> point,
	// 												 uint32_t n)
	// {
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(Grid->getDim(),n);
	// 	//	Get the nearest neighbors for the point 'index'
	// 	std::vector<uint64_t> neighbors = Grid->queryNeighbors(point,_params.k);
	// 	//	Construct the taylor matrix B
	// 	Matrix<T> B = constructTaylorMatrix(Grid,neighbors,point,n);
	// 	//	Construct the weight matrix W
	// 	Matrix<T> W = constructWeightMatrix(Grid,neighbors,point);
	// 	//	Construct the vector of corresponding Field values for neighboPointrs
	// 	Vector<T> Field_neighbors(_params.k,0.0);
	// 	for (uint32_t i = 0; i < _params.k; i++)
	// 	{
	// 		Field_neighbors(i) = (*Field)(neighbors[i]);
	// 	}
	// 	//	Complete the vanilla least squares to get the coefficients
	// 	//	compute the transpose of B
	// 	Matrix<T> B_T = B.transpose();
	// 	//	compute the product (B^T*B)
	// 	Matrix<T> B_TWB = B_T * W * B;
	// 	//	compute the inverse of the product
	// 	Matrix<T> B_TWB_inv = DGETRI(B_TWB);
	// 	//	finally take the product
	// 	Vector<T> coefficients = B_TWB_inv * B_T * W * Field_neighbors;
	// 	if (B.getFlag() == -1)
	// 	{
	// 		m_log->ERROR(B.getInfo());
	// 		m_log->ERROR(B_TWB.getInfo());
	// 	}
	// 	return coefficients;
	// }
	// //----------------------------------------------------------------------------
	// //  scalarDerivativeWMLS      - approximate the derivative for an entire
	// //                             Field using the weighted MLS method
	// //  Arguments:  Grid        - Grid<T> pointer
	// //              Field        - ScalarField<T> pointer
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<Vector<T>>
	// Interpolator<T>::scalarDerivativeWMLS(const std::shared_ptr<Grid<T>> Grid,
	// 										const std::shared_ptr<ScalarField<T>> Field,
	// 										uint32_t n)
	// {
	// 	std::vector<Vector<T>> result(Field->getN());
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(Grid->getDim(),n);
	// 	//	Query the nearest neighbors of Grid
	// 	Grid->queryNeighbors(_params.k);
	// 	for (uint64_t i = 0; i < Field->getN(); i++)
	// 	{
	// 		//	Get the nearest neighbors for the point 'index'
	// 		std::vector<uint64_t> neighbors = Grid->getNeighbors(i);
	// 		//	Construct the taylor matrix B
	// 		Matrix<T> B = constructTaylorMatrix(Grid,neighbors,i,mono);
	// 		//	Construct the weight matrix W
	// 		Matrix<T> W = constructWeightMatrix(Grid,neighbors,i);
	// 		//	Construct the vector of corresponding Field values for neighbors
	// 		Vector<T> Field_neighbors(_params.k,0.0);
	// 		for (uint32_t i = 0; i < _params.k; i++)
	// 		{
	// 			Field_neighbors(i) = (*Field)(neighbors[i]);
	// 		}
	// 		//	compute the transpose of B
	// 		Matrix<T> B_T = B.transpose();
	// 		//	compute the product (B^T*B)
	// 		Matrix<T> B_TWB = B_T * W * B;
	// 		//	compute the inverse of the product
	// 		Matrix<T> B_TWB_inv = DGETRI(B_TWB);
	// 		//	finally take the product
	// 		Vector<T> coefficients = B_TWB_inv * B_T * W * Field_neighbors;
	// 		if (B.getFlag() == -1)
	// 		{
	// 			m_log->ERROR(B.getInfo());
	// 			m_log->ERROR(B_TWB.getInfo());
	// 		}
	// 		result[i] = coefficients;
	// 	}
	// 	return result;
	// }
	// //----------------------------------------------------------------------------
	// //----------------------------------------------------------------------------
	// //  Passing Field as const reference
	// //----------------------------------------------------------------------------
	// //----------------------------------------------------------------------------
	// //  scalarDerivativeWMLSPoint - approximate the derivative for a point
	// //                             using the WMLS method
	// //  Arguments:  Grid        - Grid<T> pointer
	// //              Field        - const ScalarField<T>& reference
	// //              index        - uint64_t
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// Vector<T>
	// Interpolator<T>::scalarDerivativeWMLSPoint(const std::shared_ptr<Grid<T>> Grid,
	// 												 const ScalarField<T>& Field,
	// 												 uint64_t index,
	// 												 uint32_t n)
	// {
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(Grid->getDim(),n);
	// 	//	Query the nearest neighbors of Grid
	// 	Grid->queryNeighbors(_params.k);
	// 	//	Get the nearest neighbors for the point 'index'
	// 	std::vector<uint64_t> neighbors = Grid->getNeighbors(index);
	// 	//	Construct the taylor matrix B
	// 	Matrix<T> B = constructTaylorMatrix(Grid,neighbors,index,mono);
	// 	//	Construct the weight matrix W
	// 	Matrix<T> W = constructWeightMatrix(Grid,neighbors,index);
	// 	//	Construct the vector of corresponding Field values for neighboPointrs
	// 	Vector<T> Field_neighbors(_params.k,0.0);
	// 	for (uint32_t i = 0; i < _params.k; i++)
	// 	{
	// 		Field_neighbors(i) = Field(neighbors[i]);
	// 	}
	// 	//	compute the transpose of B
	// 	Matrix<T> B_T = B.transpose();
	// 	//	compute the product (B^T*B)
	// 	Matrix<T> B_TWB = B_T * W * B;
	// 	//	compute the inverse of the product
	// 	Matrix<T> B_TWB_inv = DGETRI(B_TWB);
	// 	//	finally take the product
	// 	Vector<T> coefficients = B_TWB_inv * B_T * W * Field_neighbors;
	// 	if (B.getFlag() == -1)
	// 	{
	// 		m_log->ERROR(B.getInfo());
	// 		m_log->ERROR(B_TWB.getInfo());
	// 	}
	// 	return coefficients;
	// }
	// //----------------------------------------------------------------------------
	// //  scalarDerivativeWMLSPoint - approximate the derivative for a point
	// //                             using the WMLS method
	// //  Arguments:  Grid        - Grid<T> pointer
	// //              Field        - const ScalarField<T>& reference
	// //              point        - std::vector<T> of the point
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// Vector<T>
	// Interpolator<T>::scalarDerivativeWMLSPoint(const std::shared_ptr<Grid<T>> Grid,
	// 												 const ScalarField<T>& Field,
	// 												 std::vector<T> point,
	// 												 uint32_t n)
	// {
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(Grid->getDim(),n);
	// 	//	Get the nearest neighbors for the point 'index'
	// 	std::vector<uint64_t> neighbors = Grid->queryNeighbors(point,_params.k);
	// 	//	Construct the taylor matrix B
	// 	Matrix<T> B = constructTaylorMatrix(Grid,neighbors,point,n);
	// 	//	Construct the weight matrix W
	// 	Matrix<T> W = constructWeightMatrix(Grid,neighbors,point);
	// 	//	Construct the vector of corresponding Field values for neighboPointrs
	// 	Vector<T> Field_neighbors(_params.k,0.0);
	// 	for (uint32_t i = 0; i < _params.k; i++)
	// 	{
	// 		Field_neighbors(i) = Field(neighbors[i]);
	// 	}
	// 	//	Complete the vanilla least squares to get the coefficients
	// 	//	compute the transpose of B
	// 	Matrix<T> B_T = B.transpose();
	// 	//	compute the product (B^T*B)
	// 	Matrix<T> B_TWB = B_T * W * B;
	// 	//	compute the inverse of the product
	// 	Matrix<T> B_TWB_inv = DGETRI(B_TWB);
	// 	//	finally take the product
	// 	Vector<T> coefficients = B_TWB_inv * B_T * W * Field_neighbors;
	// 	if (B.getFlag() == -1)
	// 	{
	// 		m_log->ERROR(B.getInfo());
	// 		m_log->ERROR(B_TWB.getInfo());
	// 	}
	// 	return coefficients;
	// }
	// //----------------------------------------------------------------------------
	// //  scalarDerivativeWMLS      - approximate the derivative for an entire
	// //                             Field using the weighted MLS method
	// //  Arguments:  Grid        - Grid<T> pointer
	// //              Field        - const ScalarField<T>& reference
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<Vector<T>>
	// Interpolator<T>::scalarDerivativeWMLS(const std::shared_ptr<Grid<T>> Grid,
	// 										const ScalarField<T>& Field,
	// 										uint32_t n)
	// {
	// 	std::vector<Vector<T>> result(Field.getN());
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(Grid->getDim(),n);
	// 	//	Query the nearest neighbors of Grid
	// 	Grid->queryNeighbors(_params.k);
	// 	for (uint64_t i = 0; i < Field.getN(); i++)
	// 	{
	// 		//	Get the nearest neighbors for the point 'index'
	// 		std::vector<uint64_t> neighbors = Grid->getNeighbors(i);
	// 		//	Construct the taylor matrix B
	// 		Matrix<T> B = constructTaylorMatrix(Grid,neighbors,i,mono);
	// 		//	Construct the weight matrix W
	// 		Matrix<T> W = constructWeightMatrix(Grid,neighbors,i);
	// 		//	Construct the vector of corres
	// 		Vector<T> Field_neighbors(_params.k,0.0);
	// 		for (uint32_t i = 0; i < _params.k; i++)
	// 		{
	// 			Field_neighbors(i) = Field(neighbors[i]);
	// 		}
	// 		//	Complete the vanilla least squares to get the coefficients
	// 		//	compute the transpose of B
	// 		Matrix<T> B_T = B.transpose();
	// 		//	compute the product (B^T*B)
	// 		Matrix<T> B_TWB = B_T * W * B;
	// 		//	compute the inverse of the product
	// 		Matrix<T> B_TWB_inv = DGETRI(B_TWB);
	// 		//	finally take the product
	// 		Vector<T> coefficients = B_TWB_inv * B_T * W * Field_neighbors;
	// 		if (B.getFlag() == -1)
	// 		{
	// 			m_log->ERROR(B.getInfo());
	// 			m_log->ERROR(B_TWB.getInfo());
	// 		}
	// 		result[i] = coefficients;
	// 	}
	// 	return result;
	// }
	// //----------------------------------------------------------------------------

	// //----------------------------------------------------------------------------
	// //  Helper functions
	// //----------------------------------------------------------------------------
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<Vector<T>>
	// Interpolator<T>::xScalarDerivative(const std::shared_ptr<Grid<T>> Grid,
	// 																 	 const std::shared_ptr<ScalarField<T>> Field,
	// 																 	 uint32_t n)
	// {
	// 	if (_type == InterpolatorType::MLS)
	// 	{
	// 		return scalarDerivativeMLS(Grid,Field,n);
	// 	}
	// 	else if (_type == InterpolatorType::WMLS)
	// 	{
	// 		return scalarDerivativeWMLS(Grid,Field,n);
	// 	}
	// 	else
	// 	{
	// 		return scalarDerivativeLS(Grid,Field,n);
	// 	}
	// }
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<Vector<T>>
	// Interpolator<T>::xScalarDerivative(const std::shared_ptr<Grid<T>> Grid,
	// 																 	 const ScalarField<T>& Field,
	// 																 	 uint32_t n)
	// {
	// 	if (_type == InterpolatorType::MLS)
	// 	{
	// 		return scalarDerivativeMLS(Grid,Field,n);
	// 	}
	// 	else if (_type == InterpolatorType::WMLS)
	// 	{
	// 		return scalarDerivativeWMLS(Grid,Field,n);
	// 	}
	// 	else
	// 	{
	// 		return scalarDerivativeLS(Grid,Field,n);
	// 	}
	// }
	// //----------------------------------------------------------------------------
	// template<typename T>
	// Vector<T>
	// Interpolator<T>::xScalarDerivativePoint(const std::shared_ptr<Grid<T>> Grid,
	// 																 const std::shared_ptr<ScalarField<T>> Field,
	// 																 uint64_t index, uint32_t n)
	// {
	// 	if (_type == InterpolatorType::MLS)
	// 	{
	// 		return scalarDerivativeMLSPoint(Grid,Field,index,n);
	// 	}
	// 	else if (_type == InterpolatorType::WMLS)
	// 	{
	// 		return scalarDerivativeWMLSPoint(Grid,Field,index,n);
	// 	}
	// 	else
	// 	{
	// 		return scalarDerivativeLSPoint(Grid,Field,index,n);
	// 	}
	// }
	// //----------------------------------------------------------------------------
	// template<typename T>
	// Vector<T>
	// Interpolator<T>::xScalarDerivativePoint(const std::shared_ptr<Grid<T>> Grid,
	// 																 const ScalarField<T>& Field,
	// 										             uint64_t index, uint32_t n)
	// {
	// 	if (_type == InterpolatorType::MLS)
	// 	{
	// 		return scalarDerivativeMLSPoint(Grid,Field,index,n);
	// 	}
	// 	else if (_type == InterpolatorType::WMLS)
	// 	{
	// 		return scalarDerivativeWMLSPoint(Grid,Field,index,n);
	// 	}
	// 	else
	// 	{
	// 		return scalarDerivativeLSPoint(Grid,Field,index,n);
	// 	}
	// }
	// //----------------------------------------------------------------------------
	// template<typename T>
	// Vector<T>
	// Interpolator<T>::xScalarDerivativePoint(const std::shared_ptr<Grid<T>> Grid,
	// 																 const std::shared_ptr<ScalarField<T>> Field,
	// 																 std::vector<T> point, uint32_t n)
	// {
	// 	if (_type == InterpolatorType::MLS)
	// 	{
	// 		return scalarDerivativeMLSPoint(Grid,Field,point,n);
	// 	}
	// 	else if (_type == InterpolatorType::WMLS)
	// 	{
	// 		return scalarDerivativeWMLSPoint(Grid,Field,point,n);
	// 	}
	// 	else
	// 	{
	// 		return scalarDerivativeLSPoint(Grid,Field,point,n);
	// 	}
	// }
	// //----------------------------------------------------------------------------
	// template<typename T>
	// Vector<T>
	// Interpolator<T>::xScalarDerivativePoint(const std::shared_ptr<Grid<T>> Grid,
	// 																 const ScalarField<T>& Field,
	// 																 std::vector<T> point, uint32_t n)
	// {
	// 	if (_type == InterpolatorType::MLS)
	// 	{
	// 		return scalarDerivativeMLSPoint(Grid,Field,point,n);
	// 	}
	// 	else if (_type == InterpolatorType::WMLS)
	// 	{
	// 		return scalarDerivativeWMLSPoint(Grid,Field,point,n);
	// 	}
	// 	else
	// 	{
	// 		return scalarDerivativeLSPoint(Grid,Field,point,n);
	// 	}
	// }
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //	Derivatives of vector Fields
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<std::vector<T>>
	// Interpolator<T>::vectorDerivative(const std::shared_ptr<Grid<T>> Grid,
	// 								                  const VectorField<T>& Field,
	// 								                  uint32_t dir, uint32_t n)
	// {
	// 	std::vector<std::vector<T>> result(Field.getN());
	// 	Monomial mono(Grid->getDim(),n);
	// 	Grid->queryNeighbors(_params.k);
	// 	for (uint64_t p = 0; p < Field.getN(); p++)
	// 	{
	// 		std::vector<uint64_t> neighbors = Grid->getNeighbors(p);
	// 		Matrix<T> B = constructTaylorMatrix(Grid,neighbors,p,mono);
	// 		std::vector<T> Field_neighbors_x(_params.k);
	// 		std::vector<T> Field_neighbors_y(_params.k);
	// 		std::vector<T> Field_neighbors_z(_params.k);
	// 		for (uint32_t i = 0; i < _params.k; i++)
	// 		{
	// 			Field_neighbors_x[i] = Field(neighbors[i],0);
	// 			Field_neighbors_y[i] = Field(neighbors[i],1);
	// 			Field_neighbors_z[i] = Field(neighbors[i],2);
	// 		}
	// 		Vector<T> Field_vals_x(Field_neighbors_x);
	// 		Vector<T> answer_x = xGELSx(B,Field_vals_x);
	// 		if (B.getFlag() == -1)
	// 		{
	// 			m_log->ERROR(B.getInfo());
	// 		}
	// 		Vector<T> Field_vals_y(Field_neighbors_y);
	// 		Vector<T> answer_y = xGELSx(B,Field_vals_y);
	// 		if (B.getFlag() == -1)
	// 		{
	// 			m_log->ERROR(B.getInfo());
	// 		}
	// 		Vector<T> Field_vals_z(Field_neighbors_z);
	// 		Vector<T> answer_z = xGELSx(B,Field_vals_z);
	// 		if (B.getFlag() == -1)
	// 		{
	// 			m_log->ERROR(B.getInfo());
	// 		}
	// 		//	Get the index of the monomial expansion corresponding to
	// 		//	the nth-derivative in the 'dir'-direction
	// 		std::vector<uint32_t> deriv(Field.getDim(),0);
	// 		deriv[dir] = n;
	// 		uint32_t index = mono.getTaylorIndex(deriv);
	// 		result[p][0] = answer_x(index);
	// 		result[p][1] = answer_y(index);
	// 		result[p][2] = answer_z(index);
	// 	}
	// 	return result;
	// }
	// //----------------------------------------------------------------------------
  //
  // //----------------------------------------------------------------------------
  // //  Construct the Taylor expansion matrix
  // //----------------------------------------------------------------------------
  // template<typename T>
  // Matrix<T>
	// Interpolator<T>::constructTaylorMatrix(				 const std::shared_ptr<Grid<T>> Grid,
  //   const std::vector<uint64_t> neighbors, uint64_t index, uint64_t order)
  // {
  //   Monomial mono(Grid->getDim(),order);
  //   std::vector<std::vector<double>> B;
  //   for (auto i = 0; i < neighbors.size(); i++)
  //   {
  //     auto id = neighbors[i];
  //     std::vector<double>
	// 		temp = mono.taylorMonomialExpansion(Grid->getPoint(index),
  //                                         Grid->getPoint(id));
  //     B.push_back(temp);
  //   }
  //   return Matrix<T>("B",B);
  // }
  // //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
  // //  Construct the Taylor expansion matrix about a point
  // //------------------------------------				 ----------------------------------------
  // template<typename T>
  // Matrix<T>
	// Interpolator<T>::constructTaylorMatrix(const std::shared_ptr<Grid<T>> Grid,
  //   const std::vector<uint64_t> neighbors, std::vector<T> point, uint64_t order)
  // {
  //   Monomial mono(Grid->getDim(),order);
  //   std::vector<std::vector<double>> B;
  //   for (uint64_t i = 0; i < neighbors.size(); i++)
  //   {
  //     uint64_t id = neighbors[i];
  //     std::vector<double>
	// 		temp = mono.taylorMonomialExpansion(point,
  //                                         Grid->getPoint(id));
  //     B.push_back(temp);
  //   }
  //   return Matrix<T>("B",B);
  // }
  // //----------------------------------------------------------------------------
  //
  // //----------------------------------------------------------------------------
  // //  Construct the Taylor expansion matrix
  // //  with mono
  // //----------------------------------------------------------------------------
  // template<typename T>
  // Matrix<T>
	// Interpolator<T>::constructTaylorMatrix(const std::shared_ptr<Grid<T>> Grid,
  //   const std::vector<uint64_t> neighbors, uint64_t index, Monomial& mono)
  // {
  //   std::vector<std::vector<double>> B;
  //   for (uint64_t i = 0; i < neighbors.size(); i++)
  //   {
  //     uint64_t id = neighbors[i];
  //     std::vector<double>
	// 		temp = mono.taylorMonomialExpansion(Grid->getPoint(index),
  //                                         Grid->getPoint(id));
  //     B.push_back(temp);
  //   }
  //   return Matrix<T>("B",B);
  // }
  // //----------------------------------------------------------------------------

	// //----------------------------------------------------------------------------
	// //	Construct weight matrix
	// //----------------------------------------------------------------------------
	// template<typename T>
	// Matrix<T>
	// Interpolator<T>::constructWeightMatrix(const std::shared_ptr<Grid<T>> Grid,
  //   const std::vector<uint64_t> neighbors, uint64_t index)
	// {
	//   //	For a simple gaussian weight, find the distances from index
	// 	//	to each point in neighbors.
	// 	std::vector<double> distances = Grid->getDistances()[index];
	// 	Matrix<T> W(neighbors.size());
	// 	for (uint32_t i = 0; i < neighbors.size(); i++)
	// 	{
	// 		W(i,i) = exp(-.5*distances[i]);
	// 	}
	// 	return W;
	// }
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //	Construct weight matrix
	// //----------------------------------------------------------------------------
	// template<typename T>
	// Matrix<T>
	// Interpolator<T>::constructWeightMatrix(const std::shared_ptr<Grid<T>> Grid,
  //   const std::vector<uint64_t> neighbors, std::vector<T> point)
	// {
	//   //	For a simple gaussian weight, find the distances from point
	// 	//	to each point in neighbors.
	// 	std::vector<double> distances = Grid->queryDistances(point,_params.k);
	// 	Matrix<T> W(neighbors.size());
	// 	for (uint32_t i = 0; i < neighbors.size(); i++)
	// 	{
	// 		W(i,i) = exp(-.5*distances[i]);
	// 	}
	// 	return W;
	// }
	// //----------------------------------------------------------------------------


}
