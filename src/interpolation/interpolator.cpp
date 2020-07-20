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
  //  map for a string to enum of LSDriver type
  std::map<std::string, LSDriver> LSDriverMap =
  {
    { "xGELS",  LSDriver::xGELS },
    { "xGELSY", LSDriver::xGELSY },
    { "xGELSD", LSDriver::xGELSD },
    { "xGELSS", LSDriver::xGELSS }
  };
  //  map for a int to enum of LSDriver type
  std::map<int, LSDriver> LSDriverMapInt =
  {
    { 0, LSDriver::xGELS },
    { 1, LSDriver::xGELSY },
    { 2, LSDriver::xGELSD },
    { 3, LSDriver::xGELSS }
  };
  //  map for a string to enum of LSDriver type
  std::map<LSDriver, std::string> LSDriverNameMap =
  {
    { LSDriver::xGELS,  "xGELS" },
    { LSDriver::xGELSY, "xGELSY" },
    { LSDriver::xGELSD, "xGELSD" },
    { LSDriver::xGELSS, "xGELSS" }
  };
  //----------------------------------------------------------------------------
	template<typename T>
  Interpolator<T>::Interpolator() : m_name("default")
  {
    m_lsdriver = LSDriver::xGELS;
    m_grid = std::make_shared<Grid<T>>();
		m_log = std::make_shared<Log>();
		m_log->init("ET:Interpolator:default", ".logs/interpolator_default.txt");
		m_log->TRACE("Interpolator 'default' created at location "
		            + getMem(*this));
  }
	//----------------------------------------------------------------------------
  //  Destructor
  //----------------------------------------------------------------------------
  template<typename T>
  Interpolator<T>::~Interpolator()
  {
		m_log->TRACE("Interpolator '" + m_name
								+ "' destroyed at location " + getMem(*this));
	}
  //----------------------------------------------------------------------------
	//	Constructor with shared Grid
	//----------------------------------------------------------------------------
  template<typename T>
  Interpolator<T>::Interpolator(std::shared_ptr<Grid<T>> t_grid)
  : m_name("default"), m_grid(t_grid)
  {
    m_lsdriver = LSDriver::xGELS;
		m_log = std::make_shared<Log>();
    m_log->init("ET:Interpolator:default", ".logs/interpolator_default.txt");
		m_log->TRACE("Interpolator 'default' created at location "
		            + getMem(*this));
  }
	//----------------------------------------------------------------------------
	//	Constructor with shared logger
	//----------------------------------------------------------------------------
	template<typename T>
  Interpolator<T>::Interpolator(std::shared_ptr<Log> t_log)
  : m_name("default")
  {
    m_lsdriver = LSDriver::xGELS;
    m_grid = std::make_shared<Grid<T>>();
		m_log = t_log;
		m_log->TRACE("Interpolator 'default' created at location "
		            + getMem(*this));
		m_log->INFO("Logger passed to Interpolator 'default'");
  }
  //----------------------------------------------------------------------------
	//	Constructor with shared Grid and logger
	//----------------------------------------------------------------------------
  template<typename T>
  Interpolator<T>::Interpolator(std::shared_ptr<Grid<T>> t_grid,
                                std::shared_ptr<Log> t_log)
  : m_name("default"), m_grid(t_grid), m_log(t_log)
  {
    m_lsdriver = LSDriver::xGELS;
    m_log->TRACE("Interpolator 'default' created at location "
                + getMem(*this));
    m_log->INFO("Logger passed to Interpolator 'default'");
  }
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
  //  Getters and Setters
  //----------------------------------------------------------------------------
  template<typename T>
  std::string Interpolator<T>::getName() const
  {
    return m_name;
  }
  template<typename T>
  std::shared_ptr<Grid<T>> Interpolator<T>::getGrid() const
  {
    return m_grid;
  }
  template<typename T>
	std::shared_ptr<Log> Interpolator<T>::getLogger() const
	{
		return m_log;
	}
  template<typename T>
  LSDriver Interpolator<T>::getLSDriver() const
  {
    return m_lsdriver;
  }
  template<typename T>
  int Interpolator<T>::getFlag() const
  {
    return m_flag;
  }
  template<typename T>
  std::string Interpolator<T>::getInfo() const
  {
    return m_info;
  }
  template<typename T>
  void Interpolator<T>::setName(std::string t_name)
  {
    m_name = t_name;
  }

  template<typename T>
  void Interpolator<T>::setGrid(std::shared_ptr<Grid<T>> t_grid)
  {
    m_grid = t_grid;
  }
  template<typename T>
  void Interpolator<T>::setLogger(std::shared_ptr<Log> t_log)
  {
    m_log = t_log;
  }
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
  template<typename T>
  void Interpolator<T>::setFlag(int t_flag)
  {
    m_flag = t_flag;
  }
  template<typename T>
  void Interpolator<T>::setInfo(std::string t_info)
  {
    m_info = t_info;
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
  template<typename T>
  Vector<T> Interpolator<T>::derivative(const size_t t_index,
                                        const size_t t_degree)
  {
   return Vector<T>();
  }
  template<typename T>
  T Interpolator<T>::derivative(const size_t t_index,
                                const size_t t_degree,
                                const size_t t_direction)
  {
   return 0;
  }
  template<typename T>
  Vector<T> Interpolator<T>::derivative(const std::vector<T>& point,
                                        const size_t t_degree)
  {
    return Vector<T>();
  }
  template<typename T>
  T Interpolator<T>::derivative(const std::vector<T>& point,
                                const size_t t_degree,
                                const size_t t_direction)
  {
    return 0;
  }

	// //----------------------------------------------------------------------------
	// //  nth-derivatives of scalar field
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //  scalarDerivative       - approximate the derivative for an entire field
	// //  Arguments:  grid      - Grid<T> pointer
	// //              field      - ScalarField<T> pointer
	// //              n          - order of derivative
	// //
	// //  Returns:    std::vector<std::vector<T>> of the derivatives
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<std::vector<T>>
	// Interpolator<T>::scalarDerivative(const std::shared_ptr<Grid<T>> grid,
	// 												        	const std::shared_ptr<ScalarField<T>> field,
	// 													        uint32_t n)
	// {
	// 	std::vector<std::vector<T>> result(field->getN());
	// 	Monomial mono(grid->getDim(),n);
	// 	std::vector<Vector<T>> vecs = xScalarDerivative(grid,field,n);
	// 	for (uint64_t i = 0; i < field->getN(); i++)
	// 	{
	// 		result[i] = vecs[i].getVec();
	// 	}
	// 	return result;
	// }
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //  scalarDerivative       - approximate the derivative for an entire field
	// //                           of order n in the direction dir
	// //  Arguments:  grid      - Grid<T> pointer
	// //              field      - ScalarField<T> pointer
	// //              dir        - direction of the derivative
	// //              n          - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the derivative along dir.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<T>
	// Interpolator<T>::scalarDerivative(const std::shared_ptr<Grid<T>> grid,
	// 												        	const std::shared_ptr<ScalarField<T>> field,
	// 													        uint32_t dir, uint32_t n)
	// {
	// 	std::vector<T> result(field->getN());
	// 	std::vector<uint32_t> deriv(field->getDim(),0);
	// 	deriv[dir] = n;
	// 	Monomial mono(grid->getDim(),n);
	// 	uint32_t index = mono.getTaylorIndex(deriv);
	// 	std::vector<Vector<T>> vecs = xScalarDerivative(grid,field,n);
	// 	for(uint64_t i = 0; i < field->getN(); i++)
	// 	{
	// 		result[i] = vecs[i](index);
	// 	}
	// 	return result;
	// }
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //  scalarDerivative       - approximate the derivative for an entire field
	// //                           of order n in the direction dir
	// //  Arguments:  grid      - Grid<T> pointer
	// //              field      - ScalarField<T> pointer
	// //              deriv      - vector of ints denoting direction and order
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<T>
	// Interpolator<T>::scalarDerivative(const std::shared_ptr<Grid<T>> grid,
	// 												        	const std::shared_ptr<ScalarField<T>> field,
	// 													        std::vector<uint32_t> deriv)
	// {
	// 	std::vector<T> result(field->getN());
	// 	uint32_t n = std::accumulate(deriv.begin(),deriv.end(),0);
	// 	Monomial mono(grid->getDim(),n);
	// 	uint32_t index = mono.getTaylorIndex(deriv);
	// 	std::vector<Vector<T>> vecs = xScalarDerivative(grid,field,n);
	// 	for(uint64_t i = 0; i < field->getN(); i++)
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
	// //  Arguments:  grid      - Grid<T> pointer
	// //              field      - ScalarField<T> pointer
	// //              index      - index of the point
	// //              n          - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<T>
	// Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<Grid<T>> grid,
	// 																	const std::shared_ptr<ScalarField<T>> field,
	// 																	uint64_t index, uint32_t n)
	// {
	// 	std::vector<T> result(field->getDim(),0.0);
	// 	Vector<T> coefficients = xScalarDerivativePoint(grid,field,index,n);
	// 	//  Grab the derivative determined by deriv
	// 	Monomial mono(grid->getDim(),n);
	// 	//  Trim result to the first field.getDim() elements
	// 	for (uint32_t j = 0; j < field->getDim(); j++)
	// 	{
	// 		std::vector<uint32_t> deriv(field->getDim(),0);
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
	// //  Arguments:  grid      - Grid<T> pointer
	// //              field      - ScalarField<T> pointer
	// //              point      - std::vector<T> of the point
	// //              n          - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<T>
	// Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<Grid<T>> grid,
	// 																	const std::shared_ptr<ScalarField<T>> field,
	// 																	std::vector<T> point, uint32_t n)
	// {
	// 	std::vector<T> result(field->getDim(),0.0);
	// 	Vector<T> coefficients = xScalarDerivativePoint(grid,field,point,n);
	// 	//  Grab the derivative determined by deriv
	// 	Monomial mono(grid->getDim(),n);
	// 	//  Trim result to the first field.getDim() elements
	// 	for (uint32_t j = 0; j < field->getDim(); j++)
	// 	{
	// 		std::vector<uint32_t> deriv(field->getDim(),0);
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
	// //  Arguments:  grid      - Grid<T> pointer
	// //              field      - ScalarField<T> pointer
	// //              index      - index of the point
	// //              dir        - direction of the derivative
	// //              n          - order of the derivative
	// //
	// //  Returns:    T          - the gradient in direction dir and order n.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// T	Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<Grid<T>> grid,
	// 																	const std::shared_ptr<ScalarField<T>> field,
	// 																	uint64_t index, uint32_t dir, uint32_t n)
	// {
	// 	Vector<T> coefficients = xScalarDerivativePoint(grid,field,index,n);
	// 	//  Grab the derivative determined by deriv
	// 	Monomial mono(grid->getDim(),n);
	// 	std::vector<uint32_t> deriv(field->getDim(),0);
	// 	deriv[dir] = n;
	// 	uint32_t l = mono.getTaylorIndex(deriv);
	// 	return coefficients(l);
	// }
	// //----------------------------------------------------------------------------
  //
  // //----------------------------------------------------------------------------
	// //  scalarDerivativePoint  - approximate the derivative for a point
	// //                           of order n in direction dir
	// //  Arguments:  grid      - Grid<T> pointer
	// //              field      - ScalarField<T> pointer
	// //              point      - std::vector<T> of the point
	// //              dir        - direction of the derivative
	// //              n          - order of the derivative
	// //
	// //  Returns:    T          - the gradient in direction dir and order n.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// T	Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<Grid<T>> grid,
	// 																	const std::shared_ptr<ScalarField<T>> field,
	// 																	std::vector<T> point, uint32_t dir,
  //                                   uint32_t n)
	// {
	// 	Vector<T> coefficients = xScalarDerivativePoint(grid,field,point,n);
	// 	//  Grab the derivative determined by deriv
	// 	Monomial mono(grid->getDim(),n);
	// 	std::vector<uint32_t> deriv(field->getDim(),0);
	// 	deriv[dir] = n;
	// 	uint32_t l = mono.getTaylorIndex(deriv);
	// 	return coefficients(l);
	// }
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //  scalarDerivativePoint  - approximate the derivative for a point
	// //                           of order n in direction dir
	// //  Arguments:  grid      - Grid<T> pointer
	// //              field      - ScalarField<T> pointer
	// //              index      - index of the point
	// //              deriv      - vector denoting the direction and order
	// //
	// //  Returns:    T          - the gradient in direction dir and order n.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// T	Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<Grid<T>> grid,
	// 																	const std::shared_ptr<ScalarField<T>> field,
	// 																	uint64_t index, std::vector<uint32_t> deriv)
	// {
	// 	uint32_t n = std::accumulate(deriv.begin(),deriv.end(),0);
	// 	Vector<T> coefficients = xScalarDerivativePoint(grid,field,index,n);
	// 	//  Grab the derivative determined by deriv
	// 	Monomial mono(grid->getDim(),n);
	// 	uint32_t l = mono.getTaylorIndex(deriv);
	// 	return coefficients(l);
	// }
	// //----------------------------------------------------------------------------
  //
  // //----------------------------------------------------------------------------
	// //  scalarDerivativePoint  - approximate the derivative for a point
	// //                           of order n in direction dir
	// //  Arguments:  grid      - Grid<T> pointer
	// //              field      - ScalarField<T> pointer
	// //              point      - std::vector<T> of the point
  // //              deriv      - vector denoting the direction and order
	// //
	// //  Returns:    T          - the gradient in direction dir and order n.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// T	Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<Grid<T>> grid,
	// 																	const std::shared_ptr<ScalarField<T>> field,
	// 																	std::vector<T> point,
  //                                   std::vector<uint32_t> deriv)
	// {
	// 	uint32_t n = std::accumulate(deriv.begin(),deriv.end(),0);
	// 	Vector<T> coefficients = xScalarDerivativePoint(grid,field,point,n);
	// 	//  Grab the derivative determined by deriv
	// 	Monomial mono(grid->getDim(),n);
	// 	uint32_t l = mono.getTaylorIndex(deriv);
	// 	return coefficients(l);
	// }
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //  Passing field as a const reference
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //  scalarDerivative       - approximate the derivative for an entire field
	// //  Arguments:  grid      - Grid<T> pointer
	// //              field      - const ScalarField<T>& reference
	// //              n          - order of derivative
	// //
	// //  Returns:    std::vector<std::vector<T>> of the derivatives
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<std::vector<T>>
	// Interpolator<T>::scalarDerivative(const std::shared_ptr<Grid<T>> grid,
	// 												        	const ScalarField<T>& field,
	// 													        uint32_t n)
	// {
	// 	std::vector<std::vector<T>> result(field.getN());
	// 	Monomial mono(grid->getDim(),n);
	// 	std::vector<Vector<T>> vecs = xScalarDerivative(grid,field,n);
	// 	for (uint64_t i = 0; i < field.getN(); i++)
	// 	{
	// 		result[i] = vecs[i].getVec();
	// 	}
	// 	return result;
	// }
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //  scalarDerivative       - approximate the derivative for an entire field
	// //                           of order n in the direction dir
	// //  Arguments:  grid      - Grid<T> pointer
	// //              field      - const ScalarField<T>& reference
	// //              dir        - direction of the derivative
	// //              n          - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the derivative along dir.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<T>
	// Interpolator<T>::scalarDerivative(const std::shared_ptr<Grid<T>> grid,
	// 												        	const ScalarField<T>& field,
	// 													        uint32_t dir, uint32_t n)
	// {
	// 	std::vector<T> result(field.getN());
	// 	std::vector<uint32_t> deriv(field.getDim(),0);
	// 	deriv[dir] = n;
	// 	Monomial mono(grid->getDim(),n);
	// 	uint32_t index = mono.getTaylorIndex(deriv);
	// 	std::vector<Vector<T>> vecs = xScalarDerivative(grid,field,n);
	// 	for(uint64_t i = 0; i < field.getN(); i++)
	// 	{
	// 		result[i] = vecs[i](index);
	// 	}
	// 	return result;
	// }
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //  scalarDerivative       - approximate the derivative for an entire field
	// //                           of order n in the direction dir
	// //  Arguments:  grid      - Grid<T> pointer
	// //              field      - const ScalarField<T>& reference
	// //              deriv      - vector of ints denoting direction and order
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<T>
	// Interpolator<T>::scalarDerivative(const std::shared_ptr<Grid<T>> grid,
	// 												        	const ScalarField<T>& field,
	// 													        std::vector<uint32_t> deriv)
	// {
	// 	std::vector<T> result(field.getN());
	// 	uint32_t n = std::accumulate(deriv.begin(),deriv.end(),0);
	// 	Monomial mono(grid->getDim(),n);
	// 	uint32_t index = mono.getTaylorIndex(deriv);
	// 	std::vector<Vector<T>> vecs = xScalarDerivative(grid,field,n);
	// 	for(uint64_t i = 0; i < field.getN(); i++)
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
	// //  Arguments:  grid      - Grid<T> pointer
	// //              field      - const ScalarField<T>& reference
	// //              index      - index of the point
	// //              n          - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<T>
	// Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<Grid<T>> grid,
	// 												        	const ScalarField<T>& field,
	// 													        uint64_t index, uint32_t n)
	// {
	// 	std::vector<T> result(field.getDim(),0.0);
	// 	Vector<T> coefficients = xScalarDerivativePoint(grid,field,index,n);
	// 	//  Grab the derivative determined by deriv
	// 	Monomial mono(grid->getDim(),n);
	// 	//  Trim result to the first field.getDim() elements
	// 	for (uint32_t j = 0; j < field.getDim(); j++)
	// 	{
	// 		std::vector<uint32_t> deriv(field.getDim(),0);
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
	// //  Arguments:  grid      - Grid<T> pointer
	// //              field      - const ScalarField<T>& reference
	// //              point      - std::vector<T> of the point
	// //              n          - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<T>
	// Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<Grid<T>> grid,
	// 												        	const ScalarField<T>& field,
	// 													        std::vector<T> point, uint32_t n)
	// {
	// 	std::vector<T> result(field.getDim(),0.0);
	// 	Vector<T> coefficients = xScalarDerivativePoint(grid,field,point,n);
	// 	//  Grab the derivative determined by deriv
	// 	Monomial mono(grid->getDim(),n);
	// 	//  Trim result to the first field.getDim() elements
	// 	for (uint32_t j = 0; j < field.getDim(); j++)
	// 	{
	// 		std::vector<uint32_t> deriv(field.getDim(),0);
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
	// //  Arguments:  grid      - Grid<T> pointer
	// //              field      - const ScalarField<T>& pointer
	// //              index      - index of the point
	// //              dir        - direction of the derivative
	// //              n          - order of the derivative
	// //
	// //  Returns:    T          - the gradient in direction dir and order n.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// T Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<Grid<T>> grid,
	// 																	  const ScalarField<T>& field,
	// 																	  uint64_t index, uint32_t dir, uint32_t n)
	// {
	// 	Vector<T> coefficients = xScalarDerivativePoint(grid,field,index,n);
	// 	//  Grab the derivative determined by deriv
	// 	Monomial mono(grid->getDim(),n);
	// 	std::vector<uint32_t> deriv(field.getDim(),0);
	// 	deriv[dir] = n;
	// 	uint32_t l = mono.getTaylorIndex(deriv);
	// 	return coefficients(l);
	// }
	// //----------------------------------------------------------------------------
  //
  // //----------------------------------------------------------------------------
	// //  scalarDerivativePoint  - approximate the derivative for a point
	// //                           of order n in direction dir
	// //  Arguments:  grid      - Grid<T> pointer
	// //              field      - const ScalarField<T>& pointer
	// //              point      - std::vector<T> of the point
	// //              dir        - direction of the derivative
	// //              n          - order of the derivative
	// //
	// //  Returns:    T          - the gradient in direction dir and order n.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// T Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<Grid<T>> grid,
	// 																	  const ScalarField<T>& field,
	// 																	  std::vector<T> point, uint32_t dir,
  //                                     uint32_t n)
	// {
	// 	Vector<T> coefficients = xScalarDerivativePoint(grid,field,point,n);
	// 	//  Grab the derivative determined by deriv
	// 	Monomial mono(grid->getDim(),n);
	// 	std::vector<uint32_t> deriv(field.getDim(),0);
	// 	deriv[dir] = n;
	// 	uint32_t l = mono.getTaylorIndex(deriv);
	// 	return coefficients(l);
	// }
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //  scalarDerivativePoint  - approximate the derivative for a point
	// //                           of order n in direction dir
	// //  Arguments:  grid      - Grid<T> pointer
	// //              field      - const ScalarField<T>& reference
	// //              index      - index of the point
	// //              deriv      - vector denoting the direction and order
	// //
	// //  Returns:    T          - the gradient in direction dir and order n.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// T Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<Grid<T>> grid,
	// 																	  const ScalarField<T>& field,
	// 																	  uint64_t index,
	// 																		std::vector<uint32_t> deriv)
	// {
	// 	uint32_t n = std::accumulate(deriv.begin(),deriv.end(),0);
	// 	Vector<T> coefficients = xScalarDerivativePoint(grid,field,index,n);
	// 	//  Grab the derivative determined by deriv
	// 	Monomial mono(grid->getDim(),n);
	// 	uint32_t l = mono.getTaylorIndex(deriv);
	// 	return coefficients(l);
	// }
	// //----------------------------------------------------------------------------
  //
  // //----------------------------------------------------------------------------
	// //  scalarDerivativePoint  - approximate the derivative for a point
	// //                           of order n in direction dir
	// //  Arguments:  grid      - Grid<T> pointer
	// //              field      - const ScalarField<T>& reference
	// //              point      - std::vector<T> of the point
	// //              deriv      - vector denoting the direction and order
	// //
	// //  Returns:    T          - the gradient in direction dir and order n.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// T Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<Grid<T>> grid,
	// 																	  const ScalarField<T>& field,
	// 																	  std::vector<T> point,
	// 																		std::vector<uint32_t> deriv)
	// {
	// 	uint32_t n = std::accumulate(deriv.begin(),deriv.end(),0);
	// 	Vector<T> coefficients = xScalarDerivativePoint(grid,field,point,n);
	// 	//  Grab the derivative determined by deriv
	// 	Monomial mono(grid->getDim(),n);
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
	// //  Arguments:  grid        - Grid<T> pointer
	// //              field        - ScalarField<T> pointer
	// //              index        - uint64_t
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// Vector<T>
	// Interpolator<T>::scalarDerivativeLSPoint(const std::shared_ptr<Grid<T>> grid,
	// 												 const std::shared_ptr<ScalarField<T>> field,
	// 												 uint64_t index,
	// 												 uint32_t n)
	// {
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(grid->getDim(),n);
	// 	//	Query the nearest neighbors of grid
	// 	grid->queryNeighbors(_params.k);
	// 	//	Get the nearest neighbors for the point 'index'
	// 	std::vector<uint64_t> neighbors = grid->getNeighbors(index);
	// 	//	Construct the taylor matrix B
	// 	Matrix<T> B = constructTaylorMatrix(grid,neighbors,index,mono);
	// 	//	Construct the vector of corresponding field values for neighboPointrs
	// 	Vector<T> field_neighbors(_params.k,0.0);
	// 	for (uint32_t i = 0; i < _params.k; i++)
	// 	{
	// 		field_neighbors(i) = (*field)(neighbors[i]);
	// 	}
	// 	//	Complete the vanilla least squares to get the coefficients
	// 	Vector<T> coefficients = xGELSx(B,field_neighbors);
	// 	if (B.getFlag() == -1)
	// 	{
	// 		m_log->ERROR(B.getInfo());
	// 	}
	// 	return coefficients;
	// }
	// //----------------------------------------------------------------------------
	// //  scalarDerivativeLSPoint - approximate the derivative for a point
	// //                             using the vanilla LS method
	// //  Arguments:  grid        - Grid<T> pointer
	// //              field        - ScalarField<T> pointer
	// //              point        - std::vector<T> of the point
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// Vector<T>
	// Interpolator<T>::scalarDerivativeLSPoint(const std::shared_ptr<Grid<T>> grid,
	// 												 const std::shared_ptr<ScalarField<T>> field,
	// 												 std::vector<T> point,
	// 												 uint32_t n)
	// {
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(grid->getDim(),n);
	// 	//	Get the nearest neighbors for the point 'index'
	// 	std::vector<uint64_t> neighbors = grid->queryNeighbors(point,_params.k);
	// 	//	Construct the taylor matrix B
	// 	Matrix<T> B = constructTaylorMatrix(grid,neighbors,point,n);
	// 	//	Construct the vector of corresponding field values for neighboPointrs
	// 	Vector<T> field_neighbors(_params.k,0.0);
	// 	for (uint32_t i = 0; i < _params.k; i++)
	// 	{
	// 		field_neighbors(i) = (*field)(neighbors[i]);
	// 	}
	// 	//	Complete the vanilla least squares to get the coefficients
	// 	Vector<T> coefficients = xGELSx(B,field_neighbors);
	// 	if (B.getFlag() == -1)
	// 	{
	// 		m_log->ERROR(B.getInfo());
	// 	}
	// 	return coefficients;
	// }
	// //----------------------------------------------------------------------------
	// //  scalarDerivativeLS      - approximate the derivative for an entire
	// //                             field using the vanilla LS method
	// //  Arguments:  grid        - Grid<T> pointer
	// //              field        - ScalarField<T> pointer
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<Vector<T>>
	// Interpolator<T>::scalarDerivativeLS(const std::shared_ptr<Grid<T>> grid,
	// 										const std::shared_ptr<ScalarField<T>> field,
	// 										uint32_t n)
	// {
	// 	std::vector<Vector<T>> result(field->getN());
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(grid->getDim(),n);
	// 	//	Query the nearest neighbors of grid
	// 	grid->queryNeighbors(_params.k);
	// 	for (uint64_t i = 0; i < field->getN(); i++)
	// 	{
	// 		//	Get the nearest neighbors for the point 'index'
	// 		std::vector<uint64_t> neighbors = grid->getNeighbors(i);
	// 		//	Construct the taylor matrix B
	// 		Matrix<T> B = constructTaylorMatrix(grid,neighbors,i,mono);
	// 		//	Construct the vector of corresponding field values for neighbors
	// 		Vector<T> field_neighbors(_params.k,0.0);
	// 		for (uint32_t i = 0; i < _params.k; i++)
	// 		{
	// 			field_neighbors(i) = (*field)(neighbors[i]);
	// 		}
	// 		//	Complete the vanilla least squares to get the coefficients
	// 		Vector<T> coefficients = xGELSx(B,field_neighbors);
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
	// //  Passing field as a const reference
	// //----------------------------------------------------------------------------
	// //----------------------------------------------------------------------------
	// //  scalarDerivativeLSPoint - approximate the derivative for a point
	// //                             using the vanilla LS method
	// //  Arguments:  grid        - Grid<T> pointer
	// //              field        - ScalarField<T> pointer
	// //              index        - uint64_t
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// Vector<T>
	// Interpolator<T>::scalarDerivativeLSPoint(const std::shared_ptr<Grid<T>> grid,
	// 												 const ScalarField<T>& field,
	// 												 uint64_t index,
	// 												 uint32_t n)
	// {
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(grid->getDim(),n);
	// 	//	Query the nearest neighbors of grid
	// 	grid->queryNeighbors(_params.k);
	// 	//	Get the nearest neighbors for the point 'index'
	// 	std::vector<uint64_t> neighbors = grid->getNeighbors(index);
	// 	//	Construct the taylor matrix B
	// 	Matrix<T> B = constructTaylorMatrix(grid,neighbors,index,mono);
	// 	//	Construct the vector of corresponding field values for neighbors
	// 	Vector<T> field_neighbors(_params.k,0.0);
	// 	for (uint32_t i = 0; i < _params.k; i++)
	// 	{
	// 		field_neighbors(i) = field(neighbors[i]);
	// 	}
	// 	//	Complete the vanilla least squares to get the coefficients
	// 	Vector<T> coefficients = xGELSx(B,field_neighbors);
	// 	if (B.getFlag() == -1)
	// 	{
	// 		m_log->ERROR(B.getInfo());
	// 	}
	// 	return coefficients;
	// }
	// //----------------------------------------------------------------------------
	// //  scalarDerivativeLSPoint - approximate the derivative for a point
	// //                             using the vanilla LS method
	// //  Arguments:  grid        - Grid<T> pointer
	// //              field        - const ScalarField<T>& reference
	// //              point        - std::vector<T> of the point
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// Vector<T>
	// Interpolator<T>::scalarDerivativeLSPoint(const std::shared_ptr<Grid<T>> grid,
	// 												 const ScalarField<T>& field,
	// 												 std::vector<T> point,
	// 												 uint32_t n)
	// {
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(grid->getDim(),n);
	// 	//	Get the nearest neighbors for the point 'index'
	// 	std::vector<uint64_t> neighbors = grid->queryNeighbors(point,_params.k);
	// 	//	Construct the taylor matrix B
	// 	Matrix<T> B = constructTaylorMatrix(grid,neighbors,point,n);
	// 	//	Construct the vector of corresponding field values for neighboPointrs
	// 	Vector<T> field_neighbors(_params.k,0.0);
	// 	for (uint32_t i = 0; i < _params.k; i++)
	// 	{
	// 		field_neighbors(i) = field(neighbors[i]);
	// 	}
	// 	//	Complete the vanilla least squares to get the coefficients
	// 	Vector<T> coefficients = xGELSx(B,field_neighbors);
	// 	if (B.getFlag() == -1)
	// 	{
	// 		m_log->ERROR(B.getInfo());
	// 	}
	// 	return coefficients;
	// }
	// //----------------------------------------------------------------------------
	// //  scalarDerivativeLS      - approximate the derivative for an entire
	// //                             field using the vanilla LS method
	// //  Arguments:  grid        - Grid<T> pointer
	// //              field        - ScalarField<T> pointer
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<Vector<T>>
	// Interpolator<T>::scalarDerivativeLS(const std::shared_ptr<Grid<T>> grid,
	// 										const ScalarField<T>& field,
	// 										uint32_t n)
	// {
	// 	std::vector<Vector<T>> result(field.getN());
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(grid->getDim(),n);
	// 	//	Query the nearest neighbors of grid
	// 	grid->queryNeighbors(_params.k);
	// 	for (uint64_t i = 0; i < field.getN(); i++)
	// 	{
	// 		//	Get the nearest neighbors for the point 'index'
	// 		std::vector<uint64_t> neighbors = grid->getNeighbors(i);
	// 		//	Construct the taylor matrix B
	// 		Matrix<T> B = constructTaylorMatrix(grid,neighbors,i,mono);
	// 		//	Construct the vector of corresponding field values for neighbors
	// 		Vector<T> field_neighbors(_params.k,0.0);
	// 		for (uint32_t i = 0; i < _params.k; i++)
	// 		{
	// 			field_neighbors(i) = field(neighbors[i]);
	// 		}
	// 		//	Complete the vanilla least squares to get the coefficients
	// 		Vector<T> coefficients = xGELSx(B,field_neighbors);
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
	// //  Arguments:  grid        - Grid<T> pointer
	// //              field        - ScalarField<T> pointer
	// //              index        - uint64_t
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// Vector<T>
	// Interpolator<T>::scalarDerivativeMLSPoint(const std::shared_ptr<Grid<T>> grid,
	// 												 const std::shared_ptr<ScalarField<T>> field,
	// 												 uint64_t index,
	// 												 uint32_t n)
	// {
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(grid->getDim(),n);
	// 	//	Query the nearest neighbors of grid
	// 	grid->queryNeighbors(_params.k);
	// 	//	Get the nearest neighbors for the point 'index'
	// 	std::vector<uint64_t> neighbors = grid->getNeighbors(index);
	// 	//	Construct the taylor matrix B
	// 	Matrix<T> B = constructTaylorMatrix(grid,neighbors,index,mono);
	// 	//	Construct the vector of corresponding field values for neighboPointrs
	// 	Vector<T> field_neighbors(_params.k,0.0);
	// 	for (uint32_t i = 0; i < _params.k; i++)
	// 	{
	// 		field_neighbors(i) = (*field)(neighbors[i]);
	// 	}
	// 	//	Complete the vanilla least squares to get the coefficients
	// 	//	compute the transpose of B
	// 	Matrix<T> B_T = B.transpose();
	// 	//	compute the product (B^T*B)
	// 	Matrix<T> B_TB = B_T * B;
	// 	//	compute the inverse of the product
	// 	Matrix<T> B_TB_inv = DGETRI(B_TB);
	// 	//	finally take the product
	// 	Vector<T> coefficients = B_TB_inv * B_T * field_neighbors;
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
	// //  Arguments:  grid        - Grid<T> pointer
	// //              field        - ScalarField<T> pointer
	// //              point        - coordinates of the point of interest
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// Vector<T>
	// Interpolator<T>::scalarDerivativeMLSPoint(const std::shared_ptr<Grid<T>> grid,
	// 												 const std::shared_ptr<ScalarField<T>> field,
	// 												 std::vector<T> point,
	// 												 uint32_t n)
	// {
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(grid->getDim(),n);
	// 	//	Get the nearest neighbors for the point 'index'
	// 	std::vector<uint64_t> neighbors = grid->queryNeighbors(point,_params.k);
	// 	//	Construct the taylor matrix B
	// 	Matrix<T> B = constructTaylorMatrix(grid,neighbors,point,n);
	// 	//	Construct the vector of corresponding field values for neighboPointrs
	// 	Vector<T> field_neighbors(_params.k,0.0);
	// 	for (uint32_t i = 0; i < _params.k; i++)
	// 	{
	// 		field_neighbors(i) = (*field)(neighbors[i]);
	// 	}
	// 	//	Complete the vanilla least squares to get the coefficients
	// 	//	compute the transpose of B
	// 	Matrix<T> B_T = B.transpose();
	// 	//	compute the product (B^T*B)
	// 	Matrix<T> B_TB = B_T * B;
	// 	//	compute the inverse of the product
	// 	Matrix<T> B_TB_inv = DGETRI(B_TB);
	// 	//	finally take the product
	// 	Vector<T> coefficients = B_TB_inv * B_T * field_neighbors;
	// 	if (B.getFlag() == -1)
	// 	{
	// 		m_log->ERROR(B.getInfo());
	// 		m_log->ERROR(B_TB.getInfo());
	// 	}
	// 	return coefficients;
	// }
	// //----------------------------------------------------------------------------
	// //  scalarDerivativeMLS      - approximate the derivative for an entire
	// //                             field using the MLS method
	// //  Arguments:  grid        - Grid<T> pointer
	// //              field        - ScalarField<T> pointer
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<Vector<T>>
	// Interpolator<T>::scalarDerivativeMLS(const std::shared_ptr<Grid<T>> grid,
	// 										const std::shared_ptr<ScalarField<T>> field,
	// 										uint32_t n)
	// {
	// 	std::vector<Vector<T>> result(field->getN());
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(grid->getDim(),n);
	// 	//	Query the nearest neighbors of grid
	// 	grid->queryNeighbors(_params.k);
	// 	for (uint64_t i = 0; i < field->getN(); i++)
	// 	{
	// 		//	Get the nearest neighbors for the point 'index'
	// 		std::vector<uint64_t> neighbors = grid->getNeighbors(i);
	// 		//	Construct the taylor matrix B
	// 		Matrix<T> B = constructTaylorMatrix(grid,neighbors,i,mono);
	// 		//	Construct the vector of corresponding field values for neighbors
	// 		Vector<T> field_neighbors(_params.k,0.0);
	// 		for (uint32_t i = 0; i < _params.k; i++)
	// 		{
	// 			field_neighbors(i) = (*field)(neighbors[i]);
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
	// 		Vector<T> coefficients = B_TB_inv * B_T * field_neighbors;
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
	// //  Passing field as const reference
	// //----------------------------------------------------------------------------
	// //----------------------------------------------------------------------------
	// //  scalarDerivativeMLSPoint - approximate the derivative for a point
	// //                             using the MLS method
	// //  Arguments:  grid        - Grid<T> pointer
	// //              field        - const ScalarField<T>& reference
	// //              index        - uint64_t
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// Vector<T>
	// Interpolator<T>::scalarDerivativeMLSPoint(const std::shared_ptr<Grid<T>> grid,
	// 												 const ScalarField<T>& field,
	// 												 uint64_t index,
	// 												 uint32_t n)
	// {
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(grid->getDim(),n);
	// 	//	Query the nearest neighbors of grid
	// 	grid->queryNeighbors(_params.k);
	// 	//	Get the nearest neighbors for the point 'index'
	// 	std::vector<uint64_t> neighbors = grid->getNeighbors(index);
	// 	//	Construct the taylor matrix B
	// 	Matrix<T> B = constructTaylorMatrix(grid,neighbors,index,mono);
	// 	//	Construct the vector of corresponding field values for neighboPointrs
	// 	Vector<T> field_neighbors(_params.k,0.0);
	// 	for (uint32_t i = 0; i < _params.k; i++)
	// 	{
	// 		field_neighbors(i) = field(neighbors[i]);
	// 	}
	// 	//	Complete the vanilla least squares to get the coefficients
	// 	//	compute the transpose of B
	// 	Matrix<T> B_T = B.transpose();
	// 	//	compute the product (B^T*B)
	// 	Matrix<T> B_TB = B_T * B;
	// 	//	compute the inverse of the product
	// 	Matrix<T> B_TB_inv = DGETRI(B_TB);
	// 	//	finally take the product
	// 	Vector<T> coefficients = B_TB_inv * B_T * field_neighbors;
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
	// //  Arguments:  grid        - Grid<T> pointer
	// //              field        - const ScalarField<T>& reference
	// //              point        - coordinates of the point of interest
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// Vector<T>
	// Interpolator<T>::scalarDerivativeMLSPoint(const std::shared_ptr<Grid<T>> grid,
	// 												 const ScalarField<T>& field,
	// 												 std::vector<T> point,
	// 												 uint32_t n)
	// {
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(grid->getDim(),n);
	// 	//	Get the nearest neighbors for the point 'index'
	// 	std::vector<uint64_t> neighbors = grid->queryNeighbors(point,_params.k);
	// 	//	Construct the taylor matrix B
	// 	Matrix<T> B = constructTaylorMatrix(grid,neighbors,point,n);
	// 	//	Construct the vector of corresponding field values for neighboPointrs
	// 	Vector<T> field_neighbors(_params.k,0.0);
	// 	for (uint32_t i = 0; i < _params.k; i++)
	// 	{
	// 		field_neighbors(i) = field(neighbors[i]);
	// 	}
	// 	//	Complete the vanilla least squares to get the coefficients
	// 	//	compute the transpose of B
	// 	Matrix<T> B_T = B.transpose();
	// 	//	compute the product (B^T*B)
	// 	Matrix<T> B_TB = B_T * B;
	// 	//	compute the inverse of the product
	// 	Matrix<T> B_TB_inv = DGETRI(B_TB);
	// 	//	finally take the product
	// 	Vector<T> coefficients = B_TB_inv * B_T * field_neighbors;
	// 	if (B.getFlag() == -1)
	// 	{
	// 		m_log->ERROR(B.getInfo());
	// 		m_log->ERROR(B_TB.getInfo());
	// 	}
	// 	return coefficients;
	// }
	// //----------------------------------------------------------------------------
	// //  scalarDerivativeLS      - approximate the derivative for an entire
	// //                             field using the vanilla LS method
	// //  Arguments:  grid        - Grid<T> pointer
	// //              field        - const ScalarField<T>& reference
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<Vector<T>>
	// Interpolator<T>::scalarDerivativeMLS(const std::shared_ptr<Grid<T>> grid,
	// 										const ScalarField<T>& field,
	// 										uint32_t n)
	// {
	// 	std::vector<Vector<T>> result(field.getN());
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(grid->getDim(),n);
	// 	//	Query the nearest neighbors of grid
	// 	grid->queryNeighbors(_params.k);
	// 	for (uint64_t i = 0; i < field.getN(); i++)
	// 	{
	// 		//	Get the nearest neighbors for the point 'index'
	// 		std::vector<uint64_t> neighbors = grid->getNeighbors(i);
	// 		//	Construct the taylor matrix B
	// 		Matrix<T> B = constructTaylorMatrix(grid,neighbors,i,mono);
	// 		//	Construct the vector of corresponding field values for neighbors
	// 		Vector<T> field_neighbors(_params.k,0.0);
	// 		for (uint32_t i = 0; i < _params.k; i++)
	// 		{
	// 			field_neighbors(i) = field(neighbors[i]);
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
	// 		Vector<T> coefficients = B_TB_inv * B_T * field_neighbors;
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
	// //  Arguments:  grid        - Grid<T> pointer
	// //              field        - ScalarField<T> pointer
	// //              index        - uint64_t
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// Vector<T>
	// Interpolator<T>::scalarDerivativeWMLSPoint(const std::shared_ptr<Grid<T>> grid,
	// 												 const std::shared_ptr<ScalarField<T>> field,
	// 												 uint64_t index,
	// 												 uint32_t n)
	// {
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(grid->getDim(),n);
	// 	//	Query the nearest neighbors of grid
	// 	grid->queryNeighbors(_params.k);
	// 	//	Get the nearest neighbors for the point 'index'
	// 	std::vector<uint64_t> neighbors = grid->getNeighbors(index);
	// 	//	Construct the taylor matrix B
	// 	Matrix<T> B = constructTaylorMatrix(grid,neighbors,index,mono);
	// 	//	Construct the weight matrix W
	// 	Matrix<T> W = constructWeightMatrix(grid,neighbors,index);
	// 	//	Construct the vector of corresponding field values for neighboPointrs
	// 	Vector<T> field_neighbors(_params.k,0.0);
	// 	for (uint32_t i = 0; i < _params.k; i++)
	// 	{
	// 		field_neighbors(i) = (*field)(neighbors[i]);
	// 	}
	// 	//	Complete the vanilla least squares to get the coefficients
	// 	//	compute the transpose of B
	// 	Matrix<T> B_T = B.transpose();
	// 	//	compute the product (B^T*B)
	// 	Matrix<T> B_TWB = B_T * W * B;
	// 	//	compute the inverse of the product
	// 	Matrix<T> B_TWB_inv = DGETRI(B_TWB);
	// 	//	finally take the product
	// 	Vector<T> coefficients = B_TWB_inv * B_T * W * field_neighbors;
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
	// //  Arguments:  grid        - Grid<T> pointer
	// //              field        - ScalarField<T> pointer
	// //              point        - std::vector<T> of the point
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// Vector<T>
	// Interpolator<T>::scalarDerivativeWMLSPoint(const std::shared_ptr<Grid<T>> grid,
	// 												 const std::shared_ptr<ScalarField<T>> field,
	// 												 std::vector<T> point,
	// 												 uint32_t n)
	// {
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(grid->getDim(),n);
	// 	//	Get the nearest neighbors for the point 'index'
	// 	std::vector<uint64_t> neighbors = grid->queryNeighbors(point,_params.k);
	// 	//	Construct the taylor matrix B
	// 	Matrix<T> B = constructTaylorMatrix(grid,neighbors,point,n);
	// 	//	Construct the weight matrix W
	// 	Matrix<T> W = constructWeightMatrix(grid,neighbors,point);
	// 	//	Construct the vector of corresponding field values for neighboPointrs
	// 	Vector<T> field_neighbors(_params.k,0.0);
	// 	for (uint32_t i = 0; i < _params.k; i++)
	// 	{
	// 		field_neighbors(i) = (*field)(neighbors[i]);
	// 	}
	// 	//	Complete the vanilla least squares to get the coefficients
	// 	//	compute the transpose of B
	// 	Matrix<T> B_T = B.transpose();
	// 	//	compute the product (B^T*B)
	// 	Matrix<T> B_TWB = B_T * W * B;
	// 	//	compute the inverse of the product
	// 	Matrix<T> B_TWB_inv = DGETRI(B_TWB);
	// 	//	finally take the product
	// 	Vector<T> coefficients = B_TWB_inv * B_T * W * field_neighbors;
	// 	if (B.getFlag() == -1)
	// 	{
	// 		m_log->ERROR(B.getInfo());
	// 		m_log->ERROR(B_TWB.getInfo());
	// 	}
	// 	return coefficients;
	// }
	// //----------------------------------------------------------------------------
	// //  scalarDerivativeWMLS      - approximate the derivative for an entire
	// //                             field using the weighted MLS method
	// //  Arguments:  grid        - Grid<T> pointer
	// //              field        - ScalarField<T> pointer
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<Vector<T>>
	// Interpolator<T>::scalarDerivativeWMLS(const std::shared_ptr<Grid<T>> grid,
	// 										const std::shared_ptr<ScalarField<T>> field,
	// 										uint32_t n)
	// {
	// 	std::vector<Vector<T>> result(field->getN());
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(grid->getDim(),n);
	// 	//	Query the nearest neighbors of grid
	// 	grid->queryNeighbors(_params.k);
	// 	for (uint64_t i = 0; i < field->getN(); i++)
	// 	{
	// 		//	Get the nearest neighbors for the point 'index'
	// 		std::vector<uint64_t> neighbors = grid->getNeighbors(i);
	// 		//	Construct the taylor matrix B
	// 		Matrix<T> B = constructTaylorMatrix(grid,neighbors,i,mono);
	// 		//	Construct the weight matrix W
	// 		Matrix<T> W = constructWeightMatrix(grid,neighbors,i);
	// 		//	Construct the vector of corresponding field values for neighbors
	// 		Vector<T> field_neighbors(_params.k,0.0);
	// 		for (uint32_t i = 0; i < _params.k; i++)
	// 		{
	// 			field_neighbors(i) = (*field)(neighbors[i]);
	// 		}
	// 		//	compute the transpose of B
	// 		Matrix<T> B_T = B.transpose();
	// 		//	compute the product (B^T*B)
	// 		Matrix<T> B_TWB = B_T * W * B;
	// 		//	compute the inverse of the product
	// 		Matrix<T> B_TWB_inv = DGETRI(B_TWB);
	// 		//	finally take the product
	// 		Vector<T> coefficients = B_TWB_inv * B_T * W * field_neighbors;
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
	// //  Passing field as const reference
	// //----------------------------------------------------------------------------
	// //----------------------------------------------------------------------------
	// //  scalarDerivativeWMLSPoint - approximate the derivative for a point
	// //                             using the WMLS method
	// //  Arguments:  grid        - Grid<T> pointer
	// //              field        - const ScalarField<T>& reference
	// //              index        - uint64_t
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// Vector<T>
	// Interpolator<T>::scalarDerivativeWMLSPoint(const std::shared_ptr<Grid<T>> grid,
	// 												 const ScalarField<T>& field,
	// 												 uint64_t index,
	// 												 uint32_t n)
	// {
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(grid->getDim(),n);
	// 	//	Query the nearest neighbors of grid
	// 	grid->queryNeighbors(_params.k);
	// 	//	Get the nearest neighbors for the point 'index'
	// 	std::vector<uint64_t> neighbors = grid->getNeighbors(index);
	// 	//	Construct the taylor matrix B
	// 	Matrix<T> B = constructTaylorMatrix(grid,neighbors,index,mono);
	// 	//	Construct the weight matrix W
	// 	Matrix<T> W = constructWeightMatrix(grid,neighbors,index);
	// 	//	Construct the vector of corresponding field values for neighboPointrs
	// 	Vector<T> field_neighbors(_params.k,0.0);
	// 	for (uint32_t i = 0; i < _params.k; i++)
	// 	{
	// 		field_neighbors(i) = field(neighbors[i]);
	// 	}
	// 	//	compute the transpose of B
	// 	Matrix<T> B_T = B.transpose();
	// 	//	compute the product (B^T*B)
	// 	Matrix<T> B_TWB = B_T * W * B;
	// 	//	compute the inverse of the product
	// 	Matrix<T> B_TWB_inv = DGETRI(B_TWB);
	// 	//	finally take the product
	// 	Vector<T> coefficients = B_TWB_inv * B_T * W * field_neighbors;
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
	// //  Arguments:  grid        - Grid<T> pointer
	// //              field        - const ScalarField<T>& reference
	// //              point        - std::vector<T> of the point
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// Vector<T>
	// Interpolator<T>::scalarDerivativeWMLSPoint(const std::shared_ptr<Grid<T>> grid,
	// 												 const ScalarField<T>& field,
	// 												 std::vector<T> point,
	// 												 uint32_t n)
	// {
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(grid->getDim(),n);
	// 	//	Get the nearest neighbors for the point 'index'
	// 	std::vector<uint64_t> neighbors = grid->queryNeighbors(point,_params.k);
	// 	//	Construct the taylor matrix B
	// 	Matrix<T> B = constructTaylorMatrix(grid,neighbors,point,n);
	// 	//	Construct the weight matrix W
	// 	Matrix<T> W = constructWeightMatrix(grid,neighbors,point);
	// 	//	Construct the vector of corresponding field values for neighboPointrs
	// 	Vector<T> field_neighbors(_params.k,0.0);
	// 	for (uint32_t i = 0; i < _params.k; i++)
	// 	{
	// 		field_neighbors(i) = field(neighbors[i]);
	// 	}
	// 	//	Complete the vanilla least squares to get the coefficients
	// 	//	compute the transpose of B
	// 	Matrix<T> B_T = B.transpose();
	// 	//	compute the product (B^T*B)
	// 	Matrix<T> B_TWB = B_T * W * B;
	// 	//	compute the inverse of the product
	// 	Matrix<T> B_TWB_inv = DGETRI(B_TWB);
	// 	//	finally take the product
	// 	Vector<T> coefficients = B_TWB_inv * B_T * W * field_neighbors;
	// 	if (B.getFlag() == -1)
	// 	{
	// 		m_log->ERROR(B.getInfo());
	// 		m_log->ERROR(B_TWB.getInfo());
	// 	}
	// 	return coefficients;
	// }
	// //----------------------------------------------------------------------------
	// //  scalarDerivativeWMLS      - approximate the derivative for an entire
	// //                             field using the weighted MLS method
	// //  Arguments:  grid        - Grid<T> pointer
	// //              field        - const ScalarField<T>& reference
	// //              n            - order of the derivative
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<Vector<T>>
	// Interpolator<T>::scalarDerivativeWMLS(const std::shared_ptr<Grid<T>> grid,
	// 										const ScalarField<T>& field,
	// 										uint32_t n)
	// {
	// 	std::vector<Vector<T>> result(field.getN());
	// 	//	Generate a monomial up to order n for the vanilla method
	// 	Monomial mono(grid->getDim(),n);
	// 	//	Query the nearest neighbors of grid
	// 	grid->queryNeighbors(_params.k);
	// 	for (uint64_t i = 0; i < field.getN(); i++)
	// 	{
	// 		//	Get the nearest neighbors for the point 'index'
	// 		std::vector<uint64_t> neighbors = grid->getNeighbors(i);
	// 		//	Construct the taylor matrix B
	// 		Matrix<T> B = constructTaylorMatrix(grid,neighbors,i,mono);
	// 		//	Construct the weight matrix W
	// 		Matrix<T> W = constructWeightMatrix(grid,neighbors,i);
	// 		//	Construct the vector of corres
	// 		Vector<T> field_neighbors(_params.k,0.0);
	// 		for (uint32_t i = 0; i < _params.k; i++)
	// 		{
	// 			field_neighbors(i) = field(neighbors[i]);
	// 		}
	// 		//	Complete the vanilla least squares to get the coefficients
	// 		//	compute the transpose of B
	// 		Matrix<T> B_T = B.transpose();
	// 		//	compute the product (B^T*B)
	// 		Matrix<T> B_TWB = B_T * W * B;
	// 		//	compute the inverse of the product
	// 		Matrix<T> B_TWB_inv = DGETRI(B_TWB);
	// 		//	finally take the product
	// 		Vector<T> coefficients = B_TWB_inv * B_T * W * field_neighbors;
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
	// Interpolator<T>::xScalarDerivative(const std::shared_ptr<Grid<T>> grid,
	// 																 	 const std::shared_ptr<ScalarField<T>> field,
	// 																 	 uint32_t n)
	// {
	// 	if (_type == InterpolatorType::MLS)
	// 	{
	// 		return scalarDerivativeMLS(grid,field,n);
	// 	}
	// 	else if (_type == InterpolatorType::WMLS)
	// 	{
	// 		return scalarDerivativeWMLS(grid,field,n);
	// 	}
	// 	else
	// 	{
	// 		return scalarDerivativeLS(grid,field,n);
	// 	}
	// }
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<Vector<T>>
	// Interpolator<T>::xScalarDerivative(const std::shared_ptr<Grid<T>> grid,
	// 																 	 const ScalarField<T>& field,
	// 																 	 uint32_t n)
	// {
	// 	if (_type == InterpolatorType::MLS)
	// 	{
	// 		return scalarDerivativeMLS(grid,field,n);
	// 	}
	// 	else if (_type == InterpolatorType::WMLS)
	// 	{
	// 		return scalarDerivativeWMLS(grid,field,n);
	// 	}
	// 	else
	// 	{
	// 		return scalarDerivativeLS(grid,field,n);
	// 	}
	// }
	// //----------------------------------------------------------------------------
	// template<typename T>
	// Vector<T>
	// Interpolator<T>::xScalarDerivativePoint(const std::shared_ptr<Grid<T>> grid,
	// 																 const std::shared_ptr<ScalarField<T>> field,
	// 																 uint64_t index, uint32_t n)
	// {
	// 	if (_type == InterpolatorType::MLS)
	// 	{
	// 		return scalarDerivativeMLSPoint(grid,field,index,n);
	// 	}
	// 	else if (_type == InterpolatorType::WMLS)
	// 	{
	// 		return scalarDerivativeWMLSPoint(grid,field,index,n);
	// 	}
	// 	else
	// 	{
	// 		return scalarDerivativeLSPoint(grid,field,index,n);
	// 	}
	// }
	// //----------------------------------------------------------------------------
	// template<typename T>
	// Vector<T>
	// Interpolator<T>::xScalarDerivativePoint(const std::shared_ptr<Grid<T>> grid,
	// 																 const ScalarField<T>& field,
	// 										             uint64_t index, uint32_t n)
	// {
	// 	if (_type == InterpolatorType::MLS)
	// 	{
	// 		return scalarDerivativeMLSPoint(grid,field,index,n);
	// 	}
	// 	else if (_type == InterpolatorType::WMLS)
	// 	{
	// 		return scalarDerivativeWMLSPoint(grid,field,index,n);
	// 	}
	// 	else
	// 	{
	// 		return scalarDerivativeLSPoint(grid,field,index,n);
	// 	}
	// }
	// //----------------------------------------------------------------------------
	// template<typename T>
	// Vector<T>
	// Interpolator<T>::xScalarDerivativePoint(const std::shared_ptr<Grid<T>> grid,
	// 																 const std::shared_ptr<ScalarField<T>> field,
	// 																 std::vector<T> point, uint32_t n)
	// {
	// 	if (_type == InterpolatorType::MLS)
	// 	{
	// 		return scalarDerivativeMLSPoint(grid,field,point,n);
	// 	}
	// 	else if (_type == InterpolatorType::WMLS)
	// 	{
	// 		return scalarDerivativeWMLSPoint(grid,field,point,n);
	// 	}
	// 	else
	// 	{
	// 		return scalarDerivativeLSPoint(grid,field,point,n);
	// 	}
	// }
	// //----------------------------------------------------------------------------
	// template<typename T>
	// Vector<T>
	// Interpolator<T>::xScalarDerivativePoint(const std::shared_ptr<Grid<T>> grid,
	// 																 const ScalarField<T>& field,
	// 																 std::vector<T> point, uint32_t n)
	// {
	// 	if (_type == InterpolatorType::MLS)
	// 	{
	// 		return scalarDerivativeMLSPoint(grid,field,point,n);
	// 	}
	// 	else if (_type == InterpolatorType::WMLS)
	// 	{
	// 		return scalarDerivativeWMLSPoint(grid,field,point,n);
	// 	}
	// 	else
	// 	{
	// 		return scalarDerivativeLSPoint(grid,field,point,n);
	// 	}
	// }
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //	Derivatives of vector fields
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<std::vector<T>>
	// Interpolator<T>::vectorDerivative(const std::shared_ptr<Grid<T>> grid,
	// 								                  const VectorField<T>& field,
	// 								                  uint32_t dir, uint32_t n)
	// {
	// 	std::vector<std::vector<T>> result(field.getN());
	// 	Monomial mono(grid->getDim(),n);
	// 	grid->queryNeighbors(_params.k);
	// 	for (uint64_t p = 0; p < field.getN(); p++)
	// 	{
	// 		std::vector<uint64_t> neighbors = grid->getNeighbors(p);
	// 		Matrix<T> B = constructTaylorMatrix(grid,neighbors,p,mono);
	// 		std::vector<T> field_neighbors_x(_params.k);
	// 		std::vector<T> field_neighbors_y(_params.k);
	// 		std::vector<T> field_neighbors_z(_params.k);
	// 		for (uint32_t i = 0; i < _params.k; i++)
	// 		{
	// 			field_neighbors_x[i] = field(neighbors[i],0);
	// 			field_neighbors_y[i] = field(neighbors[i],1);
	// 			field_neighbors_z[i] = field(neighbors[i],2);
	// 		}
	// 		Vector<T> field_vals_x(field_neighbors_x);
	// 		Vector<T> answer_x = xGELSx(B,field_vals_x);
	// 		if (B.getFlag() == -1)
	// 		{
	// 			m_log->ERROR(B.getInfo());
	// 		}
	// 		Vector<T> field_vals_y(field_neighbors_y);
	// 		Vector<T> answer_y = xGELSx(B,field_vals_y);
	// 		if (B.getFlag() == -1)
	// 		{
	// 			m_log->ERROR(B.getInfo());
	// 		}
	// 		Vector<T> field_vals_z(field_neighbors_z);
	// 		Vector<T> answer_z = xGELSx(B,field_vals_z);
	// 		if (B.getFlag() == -1)
	// 		{
	// 			m_log->ERROR(B.getInfo());
	// 		}
	// 		//	Get the index of the monomial expansion corresponding to
	// 		//	the nth-derivative in the 'dir'-direction
	// 		std::vector<uint32_t> deriv(field.getDim(),0);
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
	// Interpolator<T>::constructTaylorMatrix(				 const std::shared_ptr<Grid<T>> grid,
  //   const std::vector<uint64_t> neighbors, uint64_t index, uint64_t order)
  // {
  //   Monomial mono(grid->getDim(),order);
  //   std::vector<std::vector<double>> B;
  //   for (auto i = 0; i < neighbors.size(); i++)
  //   {
  //     auto id = neighbors[i];
  //     std::vector<double>
	// 		temp = mono.taylorMonomialExpansion(grid->getPoint(index),
  //                                         grid->getPoint(id));
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
	// Interpolator<T>::constructTaylorMatrix(const std::shared_ptr<Grid<T>> grid,
  //   const std::vector<uint64_t> neighbors, std::vector<T> point, uint64_t order)
  // {
  //   Monomial mono(grid->getDim(),order);
  //   std::vector<std::vector<double>> B;
  //   for (uint64_t i = 0; i < neighbors.size(); i++)
  //   {
  //     uint64_t id = neighbors[i];
  //     std::vector<double>
	// 		temp = mono.taylorMonomialExpansion(point,
  //                                         grid->getPoint(id));
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
	// Interpolator<T>::constructTaylorMatrix(const std::shared_ptr<Grid<T>> grid,
  //   const std::vector<uint64_t> neighbors, uint64_t index, Monomial& mono)
  // {
  //   std::vector<std::vector<double>> B;
  //   for (uint64_t i = 0; i < neighbors.size(); i++)
  //   {
  //     uint64_t id = neighbors[i];
  //     std::vector<double>
	// 		temp = mono.taylorMonomialExpansion(grid->getPoint(index),
  //                                         grid->getPoint(id));
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
	// Interpolator<T>::constructWeightMatrix(const std::shared_ptr<Grid<T>> grid,
  //   const std::vector<uint64_t> neighbors, uint64_t index)
	// {
	//   //	For a simple gaussian weight, find the distances from index
	// 	//	to each point in neighbors.
	// 	std::vector<double> distances = grid->getDistances()[index];
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
	// Interpolator<T>::constructWeightMatrix(const std::shared_ptr<Grid<T>> grid,
  //   const std::vector<uint64_t> neighbors, std::vector<T> point)
	// {
	//   //	For a simple gaussian weight, find the distances from point
	// 	//	to each point in neighbors.
	// 	std::vector<double> distances = grid->queryDistances(point,_params.k);
	// 	Matrix<T> W(neighbors.size());
	// 	for (uint32_t i = 0; i < neighbors.size(); i++)
	// 	{
	// 		W(i,i) = exp(-.5*distances[i]);
	// 	}
	// 	return W;
	// }
	// //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Various functions
  //----------------------------------------------------------------------------
  template<typename T>
  const std::string Interpolator<T>::summary()
  {
		std::string sum = "---------------------------------------------------";
		sum += "\n<ET::Interpolator<double";
		sum += "> object at " + getMem(this) + ">";
		sum += "\n---------------------------------------------------";
    sum += "\n       name:  '" + m_name + "'";
    sum += "\n   LSDriver:  '" + LSDriverNameMap[m_lsdriver] + "'";
		sum += "\n---------------------------------------------------";
		sum += "\nLogger at: " + getMem(*getLogger()) + ",";
		sum += "\n   ref at: " + getMem(m_log);
		sum += "\n++++++++++++++++++++++++++++++++++++++++++++++++++++";
		return sum;
  }
  //----------------------------------------------------------------------------
}
