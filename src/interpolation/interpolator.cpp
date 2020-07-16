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
  //----------------------------------------------------------------------------
  //  Constructors
  //----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
  //  Default constructor
  //    sets name = "default"
  //----------------------------------------------------------------------------
	template<typename T>
  Interpolator<T>::Interpolator() : m_name("default")
  {
    m_lsdriver = LSDriver::xGELS;
    m_ugrid = std::make_shared<UGrid<T>>();
		m_log = std::make_shared<Log>();
		m_log->init("ET:Interpolator:default", ".logs/approx_default.txt");
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
	//	Constructor with shared UGrid
	//----------------------------------------------------------------------------
  template<typename T>
  Interpolator<T>::Interpolator(std::shared_ptr<UGrid<T>> t_ugrid)
  : m_name("default"), m_ugrid(t_ugrid)
  {
    m_lsdriver = LSDriver::xGELS;
		m_log = log;
		m_log->TRACE("Interpolator 'default' created at location "
		            + getMem(*this));
		m_log->INFO("Logger passed to Interpolator 'default'");
  }
	//----------------------------------------------------------------------------
	//	Constructor with shared logger
	//----------------------------------------------------------------------------
	template<typename T>
  Interpolator<T>::Interpolator(std::shared_ptr<Log> t_log)
  : m_name("default")
  {
    m_lsdriver = LSDriver::xGELS;
    m_ugrid = std::make_shared<UGrid<T>>();
		m_log = t_log;
		m_log->TRACE("Interpolator 'default' created at location "
		            + getMem(*this));
		m_log->INFO("Logger passed to Interpolator 'default'");
  }
  //----------------------------------------------------------------------------
	//	Constructor with shared UGrid and logger
	//----------------------------------------------------------------------------
  template<typename T>
  Interpolator<T>::Interpolator(std::shared_ptr<UGrid<T>> t_ugrid,
                                std::shared_ptr<Log> t_log)
  : m_name("default"), m_ugrid(t_ugrid)
  {
    m_lsdriver = LSDriver::xGELS;
    m_log = t_log;
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
  std::shared_ptr<UGrid<T>> Interpolator<T>::getUGrid() const
  {
    return m_ugrid;
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
  void Interpolator<T>::setUGrid(std::shared_ptr<UGrid<T>> t_ugrid)
  {
    m_ugrid = t_ugrid;
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

	//----------------------------------------------------------------------------
	//  Gradient functions
	//----------------------------------------------------------------------------

	// //----------------------------------------------------------------------------
	// //  Scalar fields
	// //----------------------------------------------------------------------------
	// //----------------------------------------------------------------------------
	// //  scalarGradientPoint - approximate the gradient at a point
	// //  Arguments:  ugrid   - UGrid<T> pointer
	// //              field   - ScalarField<T> pointer
	// //              index   - index of the point
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
  // template<typename T>
  // std::vector<T>
	// Interpolator<T>::scalarGradientPoint(const std::shared_ptr<UGrid<T>> ugrid,
  //   const std::shared_ptr<ScalarField<T>> field, uint64_t index)
  // {
  //   if (_type == InterpolatorType::LS) {
  //     return scalarGradientLSPoint(ugrid, field, index);
  //   }
  //   else {
  //     return scalarGradientLSPoint(ugrid, field, index);
  //   }
  // }
  // //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //  scalarGradientLSPoint - approximate the gradient at a point
	// //                           using the vanilla LS method
	// //  Arguments:  ugrid      - UGrid<T> pointer
	// //              field      - ScalarField<T> pointer
	// //              index      - index of the point
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
  // template<typename T>
  // std::vector<T>
	// Interpolator<T>::scalarGradientLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
  //   const std::shared_ptr<ScalarField<T>> field, uint64_t index)
  // {
  //   std::vector<T> result(field->getDim());
	// 	Monomial mono(ugrid->getDim(),_params.n);
  //   //  First, find the nearest neighbors associated to the point specified by
  //   //  index.
  //   ugrid->queryNeighbors(_params.k);
  //   std::vector<uint64_t> neighbors = ugrid->getNeighbors(index);
  //   //  Construct the matrix associated with the ugrid spacing
  //   Matrix<T> B = constructTaylorMatrix(ugrid,neighbors,index,mono);
  //   //  Construct the vector of field values associated to each point
  //   std::vector<T> field_neighbors(_params.k);
  //   for (uint32_t i = 0; i < _params.k; i++) {
  //     field_neighbors[i] = (*field)(neighbors[i]);
  //   }
  //   Vector<T> field_vals(field_neighbors);
  //   Vector<T> answer = xGELSx(B,field_vals);
	// 	if (B.getFlag() == -1) {
	// 		m_log->ERROR(B.getInfo());
	// 	}
  //   return answer.getVec();
  // }
  // //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //  scalarGradient      - approximate the gradient for an entire field
	// //  Arguments:  ugrid   - UGrid<T> pointer
	// //              field   - ScalarField<T> pointer
	// //
	// //  Returns:    std::vector<std::vector<T>> of the gradient.
	// //----------------------------------------------------------------------------
  // template<typename T>
  // std::vector<std::vector<T>>
	// Interpolator<T>::scalarGradient(const std::shared_ptr<UGrid<T>> ugrid,
  //                                 const std::shared_ptr<ScalarField<T>> field)
  // {
  //   if (_type == InterpolatorType::LS) {
  //     return scalarGradientLS(ugrid, field);
  //   }
  //   else {
  //     return scalarGradientLS(ugrid, field);
  //   }
  // }
  // //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //  scalarGradientLS      - approximate the gradient for an entire field
	// //                           using the vanilla LS method
	// //  Arguments:  ugrid      - UGrid<T> pointer
	// //              field      - ScalarField<T> pointer
	// //
	// //  Returns:    std::vector<std::vector<T>> of the gradient.
	// //----------------------------------------------------------------------------
  // template<typename T>
  // std::vector<std::vector<T>>
	// Interpolator<T>::scalarGradientLS(const std::shared_ptr<UGrid<T>> ugrid,
  //   													const std::shared_ptr<ScalarField<T>> field)
  // {
  //   std::vector<std::vector<T>> result(field->getN());
  //   Monomial mono(ugrid->getDim(),_params.n);
	// 	ugrid->queryNeighbors(_params.k);
  //   for (uint64_t i = 0; i < field->getN(); i++) {
  //     std::vector<uint64_t> neighbors = ugrid->getNeighbors(i);
  //     Matrix<T> B = constructTaylorMatrix(ugrid,neighbors,i,mono);
  //     std::vector<T> field_neighbors(_params.k);
  //
  //     for (uint32_t i = 0; i < _params.k; i++) {
  //       field_neighbors[i] = (*field)(neighbors[i]);
  //     }
  //     Vector<T> field_vals(field_neighbors);
  //     Vector<T> answer = xGELSx(B,field_vals);
	// 		if (B.getFlag() == -1) {
	// 			m_log->ERROR(B.getInfo());
	// 		}
  //     //  Trim result to the first field->getDim() elements
	// 		std::vector<T> v = answer.getVec();
  //     std::vector<T> u = {v.begin()+1,v.begin()+field->getDim()+1};
  //     result[i] = u;
  //   }
	// 	return result;
  // }
  // //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
  // //  Passing field as a const reference
  // //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //  scalarGradientPoint - approximate the gradient at a point
	// //  Arguments:  ugrid   - UGrid<T> pointer
	// //              field   - const ScalarField<T>& reference
	// //              index   - index of the point
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
  // template<typename T>
  // std::vector<T>
	// Interpolator<T>::scalarGradientPoint(const std::shared_ptr<UGrid<T>> ugrid,
  //   const ScalarField<T>& field, uint64_t index)
  // {
  //   if (_type == InterpolatorType::LS)
  //   {
  //     return scalarGradientLSPoint(ugrid, field, index);
  //   }
  //   else
  //   {
  //     return scalarGradientLSPoint(ugrid, field, index);
  //   }
  // }
  // //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //  scalarGradientLSPoint - approximate the gradient at a point
	// //                           using the vanilla LS method
	// //  Arguments:  ugrid      - UGrid<T> pointer
	// //              field      - const ScalarField<T>& reference
	// //              index      - index of the point
	// //
	// //  Returns:    std::vector<T> of the gradient.
	// //----------------------------------------------------------------------------
  // template<typename T>
  // std::vector<T>
	// Interpolator<T>::scalarGradientLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
  //   const ScalarField<T>& field, uint64_t index)
  // {
  //   std::vector<T> result(field.getDim());
	// 	Monomial mono(ugrid->getDim(),_params.n);
  //   //  First, find the nearest neighbors associated to the point specified by
  //   //  index.
  //   ugrid->queryNeighbors(_params.k);
  //   std::vector<uint64_t> neighbors = ugrid->getNeighbors(index);
  //   //  Construct the matrix associated with the ugrid spacing
  //   Matrix<T> B = constructTaylorMatrix(ugrid,neighbors,index,mono);
  //   //  Construct the vector of field values associated to each point
  //   std::vector<T> field_neighbors(_params.k);
  //   for (uint32_t i = 0; i < _params.k; i++)
  //   {
  //     field_neighbors[i] = field(neighbors[i]);
  //   }
  //   Vector<T> field_vals(field_neighbors);
  //   Vector<T> answer = xGELSx(B,field_vals);
	// 	if (B.getFlag() == -1)
	// 	{
	// 		m_log->ERROR(B.getInfo());
	// 	}
  //   return answer.getVec();
  // }
  // //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //  scalarGradient      - approximate the gradient for an entire field
	// //  Arguments:  ugrid   - UGrid<T> pointer
	// //              field   - const ScalarField<T>& reference
	// //              index   - index of the point
	// //
	// //  Returns:    std::vector<std::vector<T>> of the gradient.
	// //----------------------------------------------------------------------------
  // template<typename T>
  // std::vector<std::vector<T>>
	// Interpolator<T>::scalarGradient(const std::shared_ptr<UGrid<T>> ugrid,
  //                                 const ScalarField<T>& field)
  // {
  //   if (_type == InterpolatorType::LS)
  //   {
  //     return scalarGradientLS(ugrid, field);
  //   }
  //   else
  //   {
  //     return scalarGradientLS(ugrid, field);
  //   }
  // }
  // //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //  scalarGradientLS      - approximate the gradient for an entire field
	// //                           using the vanilla LS method
	// //  Arguments:  ugrid      - UGrid<T> pointer
	// //              field      - const ScalarField<T>& reference
	// //              index      - index of the point
	// //
	// //  Returns:    std::vector<std::vector<T>> of the gradient.
	// //----------------------------------------------------------------------------
  // template<typename T>
  // std::vector<std::vector<T>>
	// Interpolator<T>::scalarGradientLS(const std::shared_ptr<UGrid<T>> ugrid,
  //   													const ScalarField<T>& field)
  // {
  //   std::vector<std::vector<T>> result(field.getN());
  //   Monomial mono(ugrid->getDim(),_params.n);
	// 	ugrid->queryNeighbors(_params.k);
  //   for (uint64_t i = 0; i < field.getN(); i++)
  //   {
  //     std::vector<uint64_t> neighbors = ugrid->getNeighbors(i);
  //     Matrix<T> B = constructTaylorMatrix(ugrid,neighbors,i,mono);
  //     std::vector<T> field_neighbors(_params.k);
  //     for (uint32_t i = 0; i < _params.k; i++)
  //     {
  //       field_neighbors[i] = field(neighbors[i]);
  //     }
  //     Vector<T> field_vals(field_neighbors);
  //     Vector<T> answer = xGELSx(B,field_vals);
	// 		if (B.getFlag() == -1)
	// 		{
	// 			m_log->ERROR(B.getInfo());
	// 		}
  //     //  Trim result to the first field->getDim() elements
	// 		std::vector<T> v = answer.getVec();
  //     std::vector<T> u = {v.begin()+1,v.begin()+field.getDim()+1};
  //     result[i] = u;
  //   }
	// 	return result;
  // }
  // //----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  nth-derivatives of scalar field
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  scalarDerivative       - approximate the derivative for an entire field
	//  Arguments:  ugrid      - UGrid<T> pointer
	//              field      - ScalarField<T> pointer
	//              n          - order of derivative
	//
	//  Returns:    std::vector<std::vector<T>> of the derivatives
	//----------------------------------------------------------------------------
	template<typename T>
	std::vector<std::vector<T>>
	Interpolator<T>::scalarDerivative(const std::shared_ptr<UGrid<T>> ugrid,
													        	const std::shared_ptr<ScalarField<T>> field,
														        uint32_t n)
	{
		std::vector<std::vector<T>> result(field->getN());
		Monomial mono(ugrid->getDim(),n);
		std::vector<Vector<T>> vecs = xScalarDerivative(ugrid,field,n);
		for (uint64_t i = 0; i < field->getN(); i++)
		{
			result[i] = vecs[i].getVec();
		}
		return result;
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  scalarDerivative       - approximate the derivative for an entire field
	//                           of order n in the direction dir
	//  Arguments:  ugrid      - UGrid<T> pointer
	//              field      - ScalarField<T> pointer
	//              dir        - direction of the derivative
	//              n          - order of the derivative
	//
	//  Returns:    std::vector<T> of the derivative along dir.
	//----------------------------------------------------------------------------
	template<typename T>
	std::vector<T>
	Interpolator<T>::scalarDerivative(const std::shared_ptr<UGrid<T>> ugrid,
													        	const std::shared_ptr<ScalarField<T>> field,
														        uint32_t dir, uint32_t n)
	{
		std::vector<T> result(field->getN());
		std::vector<uint32_t> deriv(field->getDim(),0);
		deriv[dir] = n;
		Monomial mono(ugrid->getDim(),n);
		uint32_t index = mono.getTaylorIndex(deriv);
		std::vector<Vector<T>> vecs = xScalarDerivative(ugrid,field,n);
		for(uint64_t i = 0; i < field->getN(); i++)
		{
			result[i] = vecs[i](index);
		}
		return result;
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  scalarDerivative       - approximate the derivative for an entire field
	//                           of order n in the direction dir
	//  Arguments:  ugrid      - UGrid<T> pointer
	//              field      - ScalarField<T> pointer
	//              deriv      - vector of ints denoting direction and order
	//
	//  Returns:    std::vector<T> of the gradient.
	//----------------------------------------------------------------------------
	template<typename T>
	std::vector<T>
	Interpolator<T>::scalarDerivative(const std::shared_ptr<UGrid<T>> ugrid,
													        	const std::shared_ptr<ScalarField<T>> field,
														        std::vector<uint32_t> deriv)
	{
		std::vector<T> result(field->getN());
		uint32_t n = std::accumulate(deriv.begin(),deriv.end(),0);
		Monomial mono(ugrid->getDim(),n);
		uint32_t index = mono.getTaylorIndex(deriv);
		std::vector<Vector<T>> vecs = xScalarDerivative(ugrid,field,n);
		for(uint64_t i = 0; i < field->getN(); i++)
		{
			result[i] = vecs[i](index);
		}
		return result;
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  scalarDerivativePoint  - approximate the derivative for a point
	//                           of order n
	//  Arguments:  ugrid      - UGrid<T> pointer
	//              field      - ScalarField<T> pointer
	//              index      - index of the point
	//              n          - order of the derivative
	//
	//  Returns:    std::vector<T> of the gradient.
	//----------------------------------------------------------------------------
	template<typename T>
	std::vector<T>
	Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
																		const std::shared_ptr<ScalarField<T>> field,
																		uint64_t index, uint32_t n)
	{
		std::vector<T> result(field->getDim(),0.0);
		Vector<T> coefficients = xScalarDerivativePoint(ugrid,field,index,n);
		//  Grab the derivative determined by deriv
		Monomial mono(ugrid->getDim(),n);
		//  Trim result to the first field.getDim() elements
		for (uint32_t j = 0; j < field->getDim(); j++)
		{
			std::vector<uint32_t> deriv(field->getDim(),0);
			deriv[j] = n;
			uint32_t l = mono.getTaylorIndex(deriv);
			result[j] = coefficients(l);
		}
		return result;
	}
	//----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
	//  scalarDerivativePoint  - approximate the derivative for a point
	//                           of order n
	//  Arguments:  ugrid      - UGrid<T> pointer
	//              field      - ScalarField<T> pointer
	//              point      - std::vector<T> of the point
	//              n          - order of the derivative
	//
	//  Returns:    std::vector<T> of the gradient.
	//----------------------------------------------------------------------------
	template<typename T>
	std::vector<T>
	Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
																		const std::shared_ptr<ScalarField<T>> field,
																		std::vector<T> point, uint32_t n)
	{
		std::vector<T> result(field->getDim(),0.0);
		Vector<T> coefficients = xScalarDerivativePoint(ugrid,field,point,n);
		//  Grab the derivative determined by deriv
		Monomial mono(ugrid->getDim(),n);
		//  Trim result to the first field.getDim() elements
		for (uint32_t j = 0; j < field->getDim(); j++)
		{
			std::vector<uint32_t> deriv(field->getDim(),0);
			deriv[j] = n;
			uint32_t l = mono.getTaylorIndex(deriv);
			result[j] = coefficients(l);
		}
		return result;
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  scalarDerivativePoint  - approximate the derivative for a point
	//                           of order n in direction dir
	//  Arguments:  ugrid      - UGrid<T> pointer
	//              field      - ScalarField<T> pointer
	//              index      - index of the point
	//              dir        - direction of the derivative
	//              n          - order of the derivative
	//
	//  Returns:    T          - the gradient in direction dir and order n.
	//----------------------------------------------------------------------------
	template<typename T>
	T	Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
																		const std::shared_ptr<ScalarField<T>> field,
																		uint64_t index, uint32_t dir, uint32_t n)
	{
		Vector<T> coefficients = xScalarDerivativePoint(ugrid,field,index,n);
		//  Grab the derivative determined by deriv
		Monomial mono(ugrid->getDim(),n);
		std::vector<uint32_t> deriv(field->getDim(),0);
		deriv[dir] = n;
		uint32_t l = mono.getTaylorIndex(deriv);
		return coefficients(l);
	}
	//----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
	//  scalarDerivativePoint  - approximate the derivative for a point
	//                           of order n in direction dir
	//  Arguments:  ugrid      - UGrid<T> pointer
	//              field      - ScalarField<T> pointer
	//              point      - std::vector<T> of the point
	//              dir        - direction of the derivative
	//              n          - order of the derivative
	//
	//  Returns:    T          - the gradient in direction dir and order n.
	//----------------------------------------------------------------------------
	template<typename T>
	T	Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
																		const std::shared_ptr<ScalarField<T>> field,
																		std::vector<T> point, uint32_t dir,
                                    uint32_t n)
	{
		Vector<T> coefficients = xScalarDerivativePoint(ugrid,field,point,n);
		//  Grab the derivative determined by deriv
		Monomial mono(ugrid->getDim(),n);
		std::vector<uint32_t> deriv(field->getDim(),0);
		deriv[dir] = n;
		uint32_t l = mono.getTaylorIndex(deriv);
		return coefficients(l);
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  scalarDerivativePoint  - approximate the derivative for a point
	//                           of order n in direction dir
	//  Arguments:  ugrid      - UGrid<T> pointer
	//              field      - ScalarField<T> pointer
	//              index      - index of the point
	//              deriv      - vector denoting the direction and order
	//
	//  Returns:    T          - the gradient in direction dir and order n.
	//----------------------------------------------------------------------------
	template<typename T>
	T	Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
																		const std::shared_ptr<ScalarField<T>> field,
																		uint64_t index, std::vector<uint32_t> deriv)
	{
		uint32_t n = std::accumulate(deriv.begin(),deriv.end(),0);
		Vector<T> coefficients = xScalarDerivativePoint(ugrid,field,index,n);
		//  Grab the derivative determined by deriv
		Monomial mono(ugrid->getDim(),n);
		uint32_t l = mono.getTaylorIndex(deriv);
		return coefficients(l);
	}
	//----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
	//  scalarDerivativePoint  - approximate the derivative for a point
	//                           of order n in direction dir
	//  Arguments:  ugrid      - UGrid<T> pointer
	//              field      - ScalarField<T> pointer
	//              point      - std::vector<T> of the point
  //              deriv      - vector denoting the direction and order
	//
	//  Returns:    T          - the gradient in direction dir and order n.
	//----------------------------------------------------------------------------
	template<typename T>
	T	Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
																		const std::shared_ptr<ScalarField<T>> field,
																		std::vector<T> point,
                                    std::vector<uint32_t> deriv)
	{
		uint32_t n = std::accumulate(deriv.begin(),deriv.end(),0);
		Vector<T> coefficients = xScalarDerivativePoint(ugrid,field,point,n);
		//  Grab the derivative determined by deriv
		Monomial mono(ugrid->getDim(),n);
		uint32_t l = mono.getTaylorIndex(deriv);
		return coefficients(l);
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  Passing field as a const reference
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  scalarDerivative       - approximate the derivative for an entire field
	//  Arguments:  ugrid      - UGrid<T> pointer
	//              field      - const ScalarField<T>& reference
	//              n          - order of derivative
	//
	//  Returns:    std::vector<std::vector<T>> of the derivatives
	//----------------------------------------------------------------------------
	template<typename T>
	std::vector<std::vector<T>>
	Interpolator<T>::scalarDerivative(const std::shared_ptr<UGrid<T>> ugrid,
													        	const ScalarField<T>& field,
														        uint32_t n)
	{
		std::vector<std::vector<T>> result(field.getN());
		Monomial mono(ugrid->getDim(),n);
		std::vector<Vector<T>> vecs = xScalarDerivative(ugrid,field,n);
		for (uint64_t i = 0; i < field.getN(); i++)
		{
			result[i] = vecs[i].getVec();
		}
		return result;
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  scalarDerivative       - approximate the derivative for an entire field
	//                           of order n in the direction dir
	//  Arguments:  ugrid      - UGrid<T> pointer
	//              field      - const ScalarField<T>& reference
	//              dir        - direction of the derivative
	//              n          - order of the derivative
	//
	//  Returns:    std::vector<T> of the derivative along dir.
	//----------------------------------------------------------------------------
	template<typename T>
	std::vector<T>
	Interpolator<T>::scalarDerivative(const std::shared_ptr<UGrid<T>> ugrid,
													        	const ScalarField<T>& field,
														        uint32_t dir, uint32_t n)
	{
		std::vector<T> result(field.getN());
		std::vector<uint32_t> deriv(field.getDim(),0);
		deriv[dir] = n;
		Monomial mono(ugrid->getDim(),n);
		uint32_t index = mono.getTaylorIndex(deriv);
		std::vector<Vector<T>> vecs = xScalarDerivative(ugrid,field,n);
		for(uint64_t i = 0; i < field.getN(); i++)
		{
			result[i] = vecs[i](index);
		}
		return result;
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  scalarDerivative       - approximate the derivative for an entire field
	//                           of order n in the direction dir
	//  Arguments:  ugrid      - UGrid<T> pointer
	//              field      - const ScalarField<T>& reference
	//              deriv      - vector of ints denoting direction and order
	//
	//  Returns:    std::vector<T> of the gradient.
	//----------------------------------------------------------------------------
	template<typename T>
	std::vector<T>
	Interpolator<T>::scalarDerivative(const std::shared_ptr<UGrid<T>> ugrid,
													        	const ScalarField<T>& field,
														        std::vector<uint32_t> deriv)
	{
		std::vector<T> result(field.getN());
		uint32_t n = std::accumulate(deriv.begin(),deriv.end(),0);
		Monomial mono(ugrid->getDim(),n);
		uint32_t index = mono.getTaylorIndex(deriv);
		std::vector<Vector<T>> vecs = xScalarDerivative(ugrid,field,n);
		for(uint64_t i = 0; i < field.getN(); i++)
		{
			result[i] = vecs[i](index);
		}
		return result;
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  scalarDerivativePoint  - approximate the derivative for a point
	//                           of order n
	//  Arguments:  ugrid      - UGrid<T> pointer
	//              field      - const ScalarField<T>& reference
	//              index      - index of the point
	//              n          - order of the derivative
	//
	//  Returns:    std::vector<T> of the gradient.
	//----------------------------------------------------------------------------
	template<typename T>
	std::vector<T>
	Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
													        	const ScalarField<T>& field,
														        uint64_t index, uint32_t n)
	{
		std::vector<T> result(field.getDim(),0.0);
		Vector<T> coefficients = xScalarDerivativePoint(ugrid,field,index,n);
		//  Grab the derivative determined by deriv
		Monomial mono(ugrid->getDim(),n);
		//  Trim result to the first field.getDim() elements
		for (uint32_t j = 0; j < field.getDim(); j++)
		{
			std::vector<uint32_t> deriv(field.getDim(),0);
			deriv[j] = n;
			uint32_t l = mono.getTaylorIndex(deriv);
			result[j] = coefficients(l);
		}
		return result;
	}
	//----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
	//  scalarDerivativePoint  - approximate the derivative for a point
	//                           of order n
	//  Arguments:  ugrid      - UGrid<T> pointer
	//              field      - const ScalarField<T>& reference
	//              point      - std::vector<T> of the point
	//              n          - order of the derivative
	//
	//  Returns:    std::vector<T> of the gradient.
	//----------------------------------------------------------------------------
	template<typename T>
	std::vector<T>
	Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
													        	const ScalarField<T>& field,
														        std::vector<T> point, uint32_t n)
	{
		std::vector<T> result(field.getDim(),0.0);
		Vector<T> coefficients = xScalarDerivativePoint(ugrid,field,point,n);
		//  Grab the derivative determined by deriv
		Monomial mono(ugrid->getDim(),n);
		//  Trim result to the first field.getDim() elements
		for (uint32_t j = 0; j < field.getDim(); j++)
		{
			std::vector<uint32_t> deriv(field.getDim(),0);
			deriv[j] = n;
			uint32_t l = mono.getTaylorIndex(deriv);
			result[j] = coefficients(l);
		}
		return result;
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  scalarDerivativePoint  - approximate the derivative for a point
	//                           of order n in direction dir
	//  Arguments:  ugrid      - UGrid<T> pointer
	//              field      - const ScalarField<T>& pointer
	//              index      - index of the point
	//              dir        - direction of the derivative
	//              n          - order of the derivative
	//
	//  Returns:    T          - the gradient in direction dir and order n.
	//----------------------------------------------------------------------------
	template<typename T>
	T Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
																		  const ScalarField<T>& field,
																		  uint64_t index, uint32_t dir, uint32_t n)
	{
		Vector<T> coefficients = xScalarDerivativePoint(ugrid,field,index,n);
		//  Grab the derivative determined by deriv
		Monomial mono(ugrid->getDim(),n);
		std::vector<uint32_t> deriv(field.getDim(),0);
		deriv[dir] = n;
		uint32_t l = mono.getTaylorIndex(deriv);
		return coefficients(l);
	}
	//----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
	//  scalarDerivativePoint  - approximate the derivative for a point
	//                           of order n in direction dir
	//  Arguments:  ugrid      - UGrid<T> pointer
	//              field      - const ScalarField<T>& pointer
	//              point      - std::vector<T> of the point
	//              dir        - direction of the derivative
	//              n          - order of the derivative
	//
	//  Returns:    T          - the gradient in direction dir and order n.
	//----------------------------------------------------------------------------
	template<typename T>
	T Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
																		  const ScalarField<T>& field,
																		  std::vector<T> point, uint32_t dir,
                                      uint32_t n)
	{
		Vector<T> coefficients = xScalarDerivativePoint(ugrid,field,point,n);
		//  Grab the derivative determined by deriv
		Monomial mono(ugrid->getDim(),n);
		std::vector<uint32_t> deriv(field.getDim(),0);
		deriv[dir] = n;
		uint32_t l = mono.getTaylorIndex(deriv);
		return coefficients(l);
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  scalarDerivativePoint  - approximate the derivative for a point
	//                           of order n in direction dir
	//  Arguments:  ugrid      - UGrid<T> pointer
	//              field      - const ScalarField<T>& reference
	//              index      - index of the point
	//              deriv      - vector denoting the direction and order
	//
	//  Returns:    T          - the gradient in direction dir and order n.
	//----------------------------------------------------------------------------
	template<typename T>
	T Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
																		  const ScalarField<T>& field,
																		  uint64_t index,
																			std::vector<uint32_t> deriv)
	{
		uint32_t n = std::accumulate(deriv.begin(),deriv.end(),0);
		Vector<T> coefficients = xScalarDerivativePoint(ugrid,field,index,n);
		//  Grab the derivative determined by deriv
		Monomial mono(ugrid->getDim(),n);
		uint32_t l = mono.getTaylorIndex(deriv);
		return coefficients(l);
	}
	//----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
	//  scalarDerivativePoint  - approximate the derivative for a point
	//                           of order n in direction dir
	//  Arguments:  ugrid      - UGrid<T> pointer
	//              field      - const ScalarField<T>& reference
	//              point      - std::vector<T> of the point
	//              deriv      - vector denoting the direction and order
	//
	//  Returns:    T          - the gradient in direction dir and order n.
	//----------------------------------------------------------------------------
	template<typename T>
	T Interpolator<T>::scalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
																		  const ScalarField<T>& field,
																		  std::vector<T> point,
																			std::vector<uint32_t> deriv)
	{
		uint32_t n = std::accumulate(deriv.begin(),deriv.end(),0);
		Vector<T> coefficients = xScalarDerivativePoint(ugrid,field,point,n);
		//  Grab the derivative determined by deriv
		Monomial mono(ugrid->getDim(),n);
		uint32_t l = mono.getTaylorIndex(deriv);
		return coefficients(l);
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  Driver routines for derivatives
	//----------------------------------------------------------------------------
	//----------------------------------------------------------------------------
	//  Vanilla least squares
	//----------------------------------------------------------------------------
	//----------------------------------------------------------------------------
	//  scalarDerivativeLSPoint - approximate the derivative for a point
	//                             using the vanilla LS method
	//  Arguments:  ugrid        - UGrid<T> pointer
	//              field        - ScalarField<T> pointer
	//              index        - uint64_t
	//              n            - order of the derivative
	//
	//  Returns:    std::vector<T> of the gradient.
	//----------------------------------------------------------------------------
	template<typename T>
	Vector<T>
	Interpolator<T>::scalarDerivativeLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
													 const std::shared_ptr<ScalarField<T>> field,
													 uint64_t index,
													 uint32_t n)
	{
		//	Generate a monomial up to order n for the vanilla method
		Monomial mono(ugrid->getDim(),n);
		//	Query the nearest neighbors of ugrid
		ugrid->queryNeighbors(_params.k);
		//	Get the nearest neighbors for the point 'index'
		std::vector<uint64_t> neighbors = ugrid->getNeighbors(index);
		//	Construct the taylor matrix B
		Matrix<T> B = constructTaylorMatrix(ugrid,neighbors,index,mono);
		//	Construct the vector of corresponding field values for neighboPointrs
		Vector<T> field_neighbors(_params.k,0.0);
		for (uint32_t i = 0; i < _params.k; i++)
		{
			field_neighbors(i) = (*field)(neighbors[i]);
		}
		//	Complete the vanilla least squares to get the coefficients
		Vector<T> coefficients = xGELSx(B,field_neighbors);
		if (B.getFlag() == -1)
		{
			m_log->ERROR(B.getInfo());
		}
		return coefficients;
	}
	//----------------------------------------------------------------------------
	//  scalarDerivativeLSPoint - approximate the derivative for a point
	//                             using the vanilla LS method
	//  Arguments:  ugrid        - UGrid<T> pointer
	//              field        - ScalarField<T> pointer
	//              point        - std::vector<T> of the point
	//              n            - order of the derivative
	//
	//  Returns:    std::vector<T> of the gradient.
	//----------------------------------------------------------------------------
	template<typename T>
	Vector<T>
	Interpolator<T>::scalarDerivativeLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
													 const std::shared_ptr<ScalarField<T>> field,
													 std::vector<T> point,
													 uint32_t n)
	{
		//	Generate a monomial up to order n for the vanilla method
		Monomial mono(ugrid->getDim(),n);
		//	Get the nearest neighbors for the point 'index'
		std::vector<uint64_t> neighbors = ugrid->queryNeighbors(point,_params.k);
		//	Construct the taylor matrix B
		Matrix<T> B = constructTaylorMatrix(ugrid,neighbors,point,n);
		//	Construct the vector of corresponding field values for neighboPointrs
		Vector<T> field_neighbors(_params.k,0.0);
		for (uint32_t i = 0; i < _params.k; i++)
		{
			field_neighbors(i) = (*field)(neighbors[i]);
		}
		//	Complete the vanilla least squares to get the coefficients
		Vector<T> coefficients = xGELSx(B,field_neighbors);
		if (B.getFlag() == -1)
		{
			m_log->ERROR(B.getInfo());
		}
		return coefficients;
	}
	//----------------------------------------------------------------------------
	//  scalarDerivativeLS      - approximate the derivative for an entire
	//                             field using the vanilla LS method
	//  Arguments:  ugrid        - UGrid<T> pointer
	//              field        - ScalarField<T> pointer
	//              n            - order of the derivative
	//
	//  Returns:    std::vector<T> of the gradient.
	//----------------------------------------------------------------------------
	template<typename T>
	std::vector<Vector<T>>
	Interpolator<T>::scalarDerivativeLS(const std::shared_ptr<UGrid<T>> ugrid,
											const std::shared_ptr<ScalarField<T>> field,
											uint32_t n)
	{
		std::vector<Vector<T>> result(field->getN());
		//	Generate a monomial up to order n for the vanilla method
		Monomial mono(ugrid->getDim(),n);
		//	Query the nearest neighbors of ugrid
		ugrid->queryNeighbors(_params.k);
		for (uint64_t i = 0; i < field->getN(); i++)
		{
			//	Get the nearest neighbors for the point 'index'
			std::vector<uint64_t> neighbors = ugrid->getNeighbors(i);
			//	Construct the taylor matrix B
			Matrix<T> B = constructTaylorMatrix(ugrid,neighbors,i,mono);
			//	Construct the vector of corresponding field values for neighbors
			Vector<T> field_neighbors(_params.k,0.0);
			for (uint32_t i = 0; i < _params.k; i++)
			{
				field_neighbors(i) = (*field)(neighbors[i]);
			}
			//	Complete the vanilla least squares to get the coefficients
			Vector<T> coefficients = xGELSx(B,field_neighbors);
			if (B.getFlag() == -1)
			{
				m_log->ERROR(B.getInfo());
			}
			result[i] = coefficients;
		}
		return result;
	}
	//----------------------------------------------------------------------------
	//----------------------------------------------------------------------------
	//  Passing field as a const reference
	//----------------------------------------------------------------------------
	//----------------------------------------------------------------------------
	//  scalarDerivativeLSPoint - approximate the derivative for a point
	//                             using the vanilla LS method
	//  Arguments:  ugrid        - UGrid<T> pointer
	//              field        - ScalarField<T> pointer
	//              index        - uint64_t
	//              n            - order of the derivative
	//
	//  Returns:    std::vector<T> of the gradient.
	//----------------------------------------------------------------------------
	template<typename T>
	Vector<T>
	Interpolator<T>::scalarDerivativeLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
													 const ScalarField<T>& field,
													 uint64_t index,
													 uint32_t n)
	{
		//	Generate a monomial up to order n for the vanilla method
		Monomial mono(ugrid->getDim(),n);
		//	Query the nearest neighbors of ugrid
		ugrid->queryNeighbors(_params.k);
		//	Get the nearest neighbors for the point 'index'
		std::vector<uint64_t> neighbors = ugrid->getNeighbors(index);
		//	Construct the taylor matrix B
		Matrix<T> B = constructTaylorMatrix(ugrid,neighbors,index,mono);
		//	Construct the vector of corresponding field values for neighbors
		Vector<T> field_neighbors(_params.k,0.0);
		for (uint32_t i = 0; i < _params.k; i++)
		{
			field_neighbors(i) = field(neighbors[i]);
		}
		//	Complete the vanilla least squares to get the coefficients
		Vector<T> coefficients = xGELSx(B,field_neighbors);
		if (B.getFlag() == -1)
		{
			m_log->ERROR(B.getInfo());
		}
		return coefficients;
	}
	//----------------------------------------------------------------------------
	//  scalarDerivativeLSPoint - approximate the derivative for a point
	//                             using the vanilla LS method
	//  Arguments:  ugrid        - UGrid<T> pointer
	//              field        - const ScalarField<T>& reference
	//              point        - std::vector<T> of the point
	//              n            - order of the derivative
	//
	//  Returns:    std::vector<T> of the gradient.
	//----------------------------------------------------------------------------
	template<typename T>
	Vector<T>
	Interpolator<T>::scalarDerivativeLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
													 const ScalarField<T>& field,
													 std::vector<T> point,
													 uint32_t n)
	{
		//	Generate a monomial up to order n for the vanilla method
		Monomial mono(ugrid->getDim(),n);
		//	Get the nearest neighbors for the point 'index'
		std::vector<uint64_t> neighbors = ugrid->queryNeighbors(point,_params.k);
		//	Construct the taylor matrix B
		Matrix<T> B = constructTaylorMatrix(ugrid,neighbors,point,n);
		//	Construct the vector of corresponding field values for neighboPointrs
		Vector<T> field_neighbors(_params.k,0.0);
		for (uint32_t i = 0; i < _params.k; i++)
		{
			field_neighbors(i) = field(neighbors[i]);
		}
		//	Complete the vanilla least squares to get the coefficients
		Vector<T> coefficients = xGELSx(B,field_neighbors);
		if (B.getFlag() == -1)
		{
			m_log->ERROR(B.getInfo());
		}
		return coefficients;
	}
	//----------------------------------------------------------------------------
	//  scalarDerivativeLS      - approximate the derivative for an entire
	//                             field using the vanilla LS method
	//  Arguments:  ugrid        - UGrid<T> pointer
	//              field        - ScalarField<T> pointer
	//              n            - order of the derivative
	//
	//  Returns:    std::vector<T> of the gradient.
	//----------------------------------------------------------------------------
	template<typename T>
	std::vector<Vector<T>>
	Interpolator<T>::scalarDerivativeLS(const std::shared_ptr<UGrid<T>> ugrid,
											const ScalarField<T>& field,
											uint32_t n)
	{
		std::vector<Vector<T>> result(field.getN());
		//	Generate a monomial up to order n for the vanilla method
		Monomial mono(ugrid->getDim(),n);
		//	Query the nearest neighbors of ugrid
		ugrid->queryNeighbors(_params.k);
		for (uint64_t i = 0; i < field.getN(); i++)
		{
			//	Get the nearest neighbors for the point 'index'
			std::vector<uint64_t> neighbors = ugrid->getNeighbors(i);
			//	Construct the taylor matrix B
			Matrix<T> B = constructTaylorMatrix(ugrid,neighbors,i,mono);
			//	Construct the vector of corresponding field values for neighbors
			Vector<T> field_neighbors(_params.k,0.0);
			for (uint32_t i = 0; i < _params.k; i++)
			{
				field_neighbors(i) = field(neighbors[i]);
			}
			//	Complete the vanilla least squares to get the coefficients
			Vector<T> coefficients = xGELSx(B,field_neighbors);
			if (B.getFlag() == -1)
			{
				m_log->ERROR(B.getInfo());
			}
			result[i] = coefficients;
		}
		return result;
	}
	//----------------------------------------------------------------------------
	//----------------------------------------------------------------------------
	//  Moving least squares
	//----------------------------------------------------------------------------
	//----------------------------------------------------------------------------
	//  scalarDerivativeMLSPoint - approximate the derivative for a point
	//                             using the MLS method
	//  Arguments:  ugrid        - UGrid<T> pointer
	//              field        - ScalarField<T> pointer
	//              index        - uint64_t
	//              n            - order of the derivative
	//
	//  Returns:    std::vector<T> of the gradient.
	//----------------------------------------------------------------------------
	template<typename T>
	Vector<T>
	Interpolator<T>::scalarDerivativeMLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
													 const std::shared_ptr<ScalarField<T>> field,
													 uint64_t index,
													 uint32_t n)
	{
		//	Generate a monomial up to order n for the vanilla method
		Monomial mono(ugrid->getDim(),n);
		//	Query the nearest neighbors of ugrid
		ugrid->queryNeighbors(_params.k);
		//	Get the nearest neighbors for the point 'index'
		std::vector<uint64_t> neighbors = ugrid->getNeighbors(index);
		//	Construct the taylor matrix B
		Matrix<T> B = constructTaylorMatrix(ugrid,neighbors,index,mono);
		//	Construct the vector of corresponding field values for neighboPointrs
		Vector<T> field_neighbors(_params.k,0.0);
		for (uint32_t i = 0; i < _params.k; i++)
		{
			field_neighbors(i) = (*field)(neighbors[i]);
		}
		//	Complete the vanilla least squares to get the coefficients
		//	compute the transpose of B
		Matrix<T> B_T = B.transpose();
		//	compute the product (B^T*B)
		Matrix<T> B_TB = B_T * B;
		//	compute the inverse of the product
		Matrix<T> B_TB_inv = DGETRI(B_TB);
		//	finally take the product
		Vector<T> coefficients = B_TB_inv * B_T * field_neighbors;
		if (B.getFlag() == -1)
		{
			m_log->ERROR(B.getInfo());
			m_log->ERROR(B_TB.getInfo());
		}
		return coefficients;
	}
	//----------------------------------------------------------------------------
	//  scalarDerivativeMLSPoint - approximate the derivative for a point
	//                             using the MLS method
	//  Arguments:  ugrid        - UGrid<T> pointer
	//              field        - ScalarField<T> pointer
	//              point        - coordinates of the point of interest
	//              n            - order of the derivative
	//
	//  Returns:    std::vector<T> of the gradient.
	//----------------------------------------------------------------------------
	template<typename T>
	Vector<T>
	Interpolator<T>::scalarDerivativeMLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
													 const std::shared_ptr<ScalarField<T>> field,
													 std::vector<T> point,
													 uint32_t n)
	{
		//	Generate a monomial up to order n for the vanilla method
		Monomial mono(ugrid->getDim(),n);
		//	Get the nearest neighbors for the point 'index'
		std::vector<uint64_t> neighbors = ugrid->queryNeighbors(point,_params.k);
		//	Construct the taylor matrix B
		Matrix<T> B = constructTaylorMatrix(ugrid,neighbors,point,n);
		//	Construct the vector of corresponding field values for neighboPointrs
		Vector<T> field_neighbors(_params.k,0.0);
		for (uint32_t i = 0; i < _params.k; i++)
		{
			field_neighbors(i) = (*field)(neighbors[i]);
		}
		//	Complete the vanilla least squares to get the coefficients
		//	compute the transpose of B
		Matrix<T> B_T = B.transpose();
		//	compute the product (B^T*B)
		Matrix<T> B_TB = B_T * B;
		//	compute the inverse of the product
		Matrix<T> B_TB_inv = DGETRI(B_TB);
		//	finally take the product
		Vector<T> coefficients = B_TB_inv * B_T * field_neighbors;
		if (B.getFlag() == -1)
		{
			m_log->ERROR(B.getInfo());
			m_log->ERROR(B_TB.getInfo());
		}
		return coefficients;
	}
	//----------------------------------------------------------------------------
	//  scalarDerivativeMLS      - approximate the derivative for an entire
	//                             field using the MLS method
	//  Arguments:  ugrid        - UGrid<T> pointer
	//              field        - ScalarField<T> pointer
	//              n            - order of the derivative
	//
	//  Returns:    std::vector<T> of the gradient.
	//----------------------------------------------------------------------------
	template<typename T>
	std::vector<Vector<T>>
	Interpolator<T>::scalarDerivativeMLS(const std::shared_ptr<UGrid<T>> ugrid,
											const std::shared_ptr<ScalarField<T>> field,
											uint32_t n)
	{
		std::vector<Vector<T>> result(field->getN());
		//	Generate a monomial up to order n for the vanilla method
		Monomial mono(ugrid->getDim(),n);
		//	Query the nearest neighbors of ugrid
		ugrid->queryNeighbors(_params.k);
		for (uint64_t i = 0; i < field->getN(); i++)
		{
			//	Get the nearest neighbors for the point 'index'
			std::vector<uint64_t> neighbors = ugrid->getNeighbors(i);
			//	Construct the taylor matrix B
			Matrix<T> B = constructTaylorMatrix(ugrid,neighbors,i,mono);
			//	Construct the vector of corresponding field values for neighbors
			Vector<T> field_neighbors(_params.k,0.0);
			for (uint32_t i = 0; i < _params.k; i++)
			{
				field_neighbors(i) = (*field)(neighbors[i]);
			}
			//	Complete the vanilla least squares to get the coefficients
			//	Complete the vanilla least squares to get the coefficients
			//	compute the transpose of B
			Matrix<T> B_T = B.transpose();
			//	compute the product (B^T*B)
			Matrix<T> B_TB = B_T * B;
			//	compute the inverse of the product
			Matrix<T> B_TB_inv = DGETRI(B_TB);
			//	finally take the product
			Vector<T> coefficients = B_TB_inv * B_T * field_neighbors;
			if (B.getFlag() == -1)
			{
				m_log->ERROR(B.getInfo());
				m_log->ERROR(B_TB.getInfo());
			}
			result[i] = coefficients;
		}
		return result;
	}
	//----------------------------------------------------------------------------
	//----------------------------------------------------------------------------
	//  Passing field as const reference
	//----------------------------------------------------------------------------
	//----------------------------------------------------------------------------
	//  scalarDerivativeMLSPoint - approximate the derivative for a point
	//                             using the MLS method
	//  Arguments:  ugrid        - UGrid<T> pointer
	//              field        - const ScalarField<T>& reference
	//              index        - uint64_t
	//              n            - order of the derivative
	//
	//  Returns:    std::vector<T> of the gradient.
	//----------------------------------------------------------------------------
	template<typename T>
	Vector<T>
	Interpolator<T>::scalarDerivativeMLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
													 const ScalarField<T>& field,
													 uint64_t index,
													 uint32_t n)
	{
		//	Generate a monomial up to order n for the vanilla method
		Monomial mono(ugrid->getDim(),n);
		//	Query the nearest neighbors of ugrid
		ugrid->queryNeighbors(_params.k);
		//	Get the nearest neighbors for the point 'index'
		std::vector<uint64_t> neighbors = ugrid->getNeighbors(index);
		//	Construct the taylor matrix B
		Matrix<T> B = constructTaylorMatrix(ugrid,neighbors,index,mono);
		//	Construct the vector of corresponding field values for neighboPointrs
		Vector<T> field_neighbors(_params.k,0.0);
		for (uint32_t i = 0; i < _params.k; i++)
		{
			field_neighbors(i) = field(neighbors[i]);
		}
		//	Complete the vanilla least squares to get the coefficients
		//	compute the transpose of B
		Matrix<T> B_T = B.transpose();
		//	compute the product (B^T*B)
		Matrix<T> B_TB = B_T * B;
		//	compute the inverse of the product
		Matrix<T> B_TB_inv = DGETRI(B_TB);
		//	finally take the product
		Vector<T> coefficients = B_TB_inv * B_T * field_neighbors;
		if (B.getFlag() == -1)
		{
			m_log->ERROR(B.getInfo());
			m_log->ERROR(B_TB.getInfo());
		}
		return coefficients;
	}
	//----------------------------------------------------------------------------
	//  scalarDerivativeMLSPoint - approximate the derivative for a point
	//                             using the MLS method
	//  Arguments:  ugrid        - UGrid<T> pointer
	//              field        - const ScalarField<T>& reference
	//              point        - coordinates of the point of interest
	//              n            - order of the derivative
	//
	//  Returns:    std::vector<T> of the gradient.
	//----------------------------------------------------------------------------
	template<typename T>
	Vector<T>
	Interpolator<T>::scalarDerivativeMLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
													 const ScalarField<T>& field,
													 std::vector<T> point,
													 uint32_t n)
	{
		//	Generate a monomial up to order n for the vanilla method
		Monomial mono(ugrid->getDim(),n);
		//	Get the nearest neighbors for the point 'index'
		std::vector<uint64_t> neighbors = ugrid->queryNeighbors(point,_params.k);
		//	Construct the taylor matrix B
		Matrix<T> B = constructTaylorMatrix(ugrid,neighbors,point,n);
		//	Construct the vector of corresponding field values for neighboPointrs
		Vector<T> field_neighbors(_params.k,0.0);
		for (uint32_t i = 0; i < _params.k; i++)
		{
			field_neighbors(i) = field(neighbors[i]);
		}
		//	Complete the vanilla least squares to get the coefficients
		//	compute the transpose of B
		Matrix<T> B_T = B.transpose();
		//	compute the product (B^T*B)
		Matrix<T> B_TB = B_T * B;
		//	compute the inverse of the product
		Matrix<T> B_TB_inv = DGETRI(B_TB);
		//	finally take the product
		Vector<T> coefficients = B_TB_inv * B_T * field_neighbors;
		if (B.getFlag() == -1)
		{
			m_log->ERROR(B.getInfo());
			m_log->ERROR(B_TB.getInfo());
		}
		return coefficients;
	}
	//----------------------------------------------------------------------------
	//  scalarDerivativeLS      - approximate the derivative for an entire
	//                             field using the vanilla LS method
	//  Arguments:  ugrid        - UGrid<T> pointer
	//              field        - const ScalarField<T>& reference
	//              n            - order of the derivative
	//
	//  Returns:    std::vector<T> of the gradient.
	//----------------------------------------------------------------------------
	template<typename T>
	std::vector<Vector<T>>
	Interpolator<T>::scalarDerivativeMLS(const std::shared_ptr<UGrid<T>> ugrid,
											const ScalarField<T>& field,
											uint32_t n)
	{
		std::vector<Vector<T>> result(field.getN());
		//	Generate a monomial up to order n for the vanilla method
		Monomial mono(ugrid->getDim(),n);
		//	Query the nearest neighbors of ugrid
		ugrid->queryNeighbors(_params.k);
		for (uint64_t i = 0; i < field.getN(); i++)
		{
			//	Get the nearest neighbors for the point 'index'
			std::vector<uint64_t> neighbors = ugrid->getNeighbors(i);
			//	Construct the taylor matrix B
			Matrix<T> B = constructTaylorMatrix(ugrid,neighbors,i,mono);
			//	Construct the vector of corresponding field values for neighbors
			Vector<T> field_neighbors(_params.k,0.0);
			for (uint32_t i = 0; i < _params.k; i++)
			{
				field_neighbors(i) = field(neighbors[i]);
			}
			//	Complete the vanilla least squares to get the coefficients
			//	Complete the vanilla least squares to get the coefficients
			//	compute the transpose of B
			Matrix<T> B_T = B.transpose();
			//	compute the product (B^T*B)
			Matrix<T> B_TB = B_T * B;
			//	compute the inverse of the product
			Matrix<T> B_TB_inv = DGETRI(B_TB);
			//	finally take the product
			Vector<T> coefficients = B_TB_inv * B_T * field_neighbors;
			if (B.getFlag() == -1)
			{
				m_log->ERROR(B.getInfo());
				m_log->ERROR(B_TB.getInfo());
			}
			result[i] = coefficients;
		}
		return result;
	}
	//----------------------------------------------------------------------------
	//----------------------------------------------------------------------------
	//  Weighted Moving least squares
	//----------------------------------------------------------------------------
	//----------------------------------------------------------------------------
	//  scalarDerivativeWMLSPoint - approximate the derivative for a point
	//                             using the WMLS method
	//  Arguments:  ugrid        - UGrid<T> pointer
	//              field        - ScalarField<T> pointer
	//              index        - uint64_t
	//              n            - order of the derivative
	//
	//  Returns:    std::vector<T> of the gradient.
	//----------------------------------------------------------------------------
	template<typename T>
	Vector<T>
	Interpolator<T>::scalarDerivativeWMLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
													 const std::shared_ptr<ScalarField<T>> field,
													 uint64_t index,
													 uint32_t n)
	{
		//	Generate a monomial up to order n for the vanilla method
		Monomial mono(ugrid->getDim(),n);
		//	Query the nearest neighbors of ugrid
		ugrid->queryNeighbors(_params.k);
		//	Get the nearest neighbors for the point 'index'
		std::vector<uint64_t> neighbors = ugrid->getNeighbors(index);
		//	Construct the taylor matrix B
		Matrix<T> B = constructTaylorMatrix(ugrid,neighbors,index,mono);
		//	Construct the weight matrix W
		Matrix<T> W = constructWeightMatrix(ugrid,neighbors,index);
		//	Construct the vector of corresponding field values for neighboPointrs
		Vector<T> field_neighbors(_params.k,0.0);
		for (uint32_t i = 0; i < _params.k; i++)
		{
			field_neighbors(i) = (*field)(neighbors[i]);
		}
		//	Complete the vanilla least squares to get the coefficients
		//	compute the transpose of B
		Matrix<T> B_T = B.transpose();
		//	compute the product (B^T*B)
		Matrix<T> B_TWB = B_T * W * B;
		//	compute the inverse of the product
		Matrix<T> B_TWB_inv = DGETRI(B_TWB);
		//	finally take the product
		Vector<T> coefficients = B_TWB_inv * B_T * W * field_neighbors;
		if (B.getFlag() == -1)
		{
			m_log->ERROR(B.getInfo());
			m_log->ERROR(B_TWB.getInfo());
		}
		return coefficients;
	}
	//----------------------------------------------------------------------------
	//  scalarDerivativeWMLSPoint - approximate the derivative for a point
	//                             using the WMLS method
	//  Arguments:  ugrid        - UGrid<T> pointer
	//              field        - ScalarField<T> pointer
	//              point        - std::vector<T> of the point
	//              n            - order of the derivative
	//
	//  Returns:    std::vector<T> of the gradient.
	//----------------------------------------------------------------------------
	template<typename T>
	Vector<T>
	Interpolator<T>::scalarDerivativeWMLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
													 const std::shared_ptr<ScalarField<T>> field,
													 std::vector<T> point,
													 uint32_t n)
	{
		//	Generate a monomial up to order n for the vanilla method
		Monomial mono(ugrid->getDim(),n);
		//	Get the nearest neighbors for the point 'index'
		std::vector<uint64_t> neighbors = ugrid->queryNeighbors(point,_params.k);
		//	Construct the taylor matrix B
		Matrix<T> B = constructTaylorMatrix(ugrid,neighbors,point,n);
		//	Construct the weight matrix W
		Matrix<T> W = constructWeightMatrix(ugrid,neighbors,point);
		//	Construct the vector of corresponding field values for neighboPointrs
		Vector<T> field_neighbors(_params.k,0.0);
		for (uint32_t i = 0; i < _params.k; i++)
		{
			field_neighbors(i) = (*field)(neighbors[i]);
		}
		//	Complete the vanilla least squares to get the coefficients
		//	compute the transpose of B
		Matrix<T> B_T = B.transpose();
		//	compute the product (B^T*B)
		Matrix<T> B_TWB = B_T * W * B;
		//	compute the inverse of the product
		Matrix<T> B_TWB_inv = DGETRI(B_TWB);
		//	finally take the product
		Vector<T> coefficients = B_TWB_inv * B_T * W * field_neighbors;
		if (B.getFlag() == -1)
		{
			m_log->ERROR(B.getInfo());
			m_log->ERROR(B_TWB.getInfo());
		}
		return coefficients;
	}
	//----------------------------------------------------------------------------
	//  scalarDerivativeWMLS      - approximate the derivative for an entire
	//                             field using the weighted MLS method
	//  Arguments:  ugrid        - UGrid<T> pointer
	//              field        - ScalarField<T> pointer
	//              n            - order of the derivative
	//
	//  Returns:    std::vector<T> of the gradient.
	//----------------------------------------------------------------------------
	template<typename T>
	std::vector<Vector<T>>
	Interpolator<T>::scalarDerivativeWMLS(const std::shared_ptr<UGrid<T>> ugrid,
											const std::shared_ptr<ScalarField<T>> field,
											uint32_t n)
	{
		std::vector<Vector<T>> result(field->getN());
		//	Generate a monomial up to order n for the vanilla method
		Monomial mono(ugrid->getDim(),n);
		//	Query the nearest neighbors of ugrid
		ugrid->queryNeighbors(_params.k);
		for (uint64_t i = 0; i < field->getN(); i++)
		{
			//	Get the nearest neighbors for the point 'index'
			std::vector<uint64_t> neighbors = ugrid->getNeighbors(i);
			//	Construct the taylor matrix B
			Matrix<T> B = constructTaylorMatrix(ugrid,neighbors,i,mono);
			//	Construct the weight matrix W
			Matrix<T> W = constructWeightMatrix(ugrid,neighbors,i);
			//	Construct the vector of corresponding field values for neighbors
			Vector<T> field_neighbors(_params.k,0.0);
			for (uint32_t i = 0; i < _params.k; i++)
			{
				field_neighbors(i) = (*field)(neighbors[i]);
			}
			//	compute the transpose of B
			Matrix<T> B_T = B.transpose();
			//	compute the product (B^T*B)
			Matrix<T> B_TWB = B_T * W * B;
			//	compute the inverse of the product
			Matrix<T> B_TWB_inv = DGETRI(B_TWB);
			//	finally take the product
			Vector<T> coefficients = B_TWB_inv * B_T * W * field_neighbors;
			if (B.getFlag() == -1)
			{
				m_log->ERROR(B.getInfo());
				m_log->ERROR(B_TWB.getInfo());
			}
			result[i] = coefficients;
		}
		return result;
	}
	//----------------------------------------------------------------------------
	//----------------------------------------------------------------------------
	//  Passing field as const reference
	//----------------------------------------------------------------------------
	//----------------------------------------------------------------------------
	//  scalarDerivativeWMLSPoint - approximate the derivative for a point
	//                             using the WMLS method
	//  Arguments:  ugrid        - UGrid<T> pointer
	//              field        - const ScalarField<T>& reference
	//              index        - uint64_t
	//              n            - order of the derivative
	//
	//  Returns:    std::vector<T> of the gradient.
	//----------------------------------------------------------------------------
	template<typename T>
	Vector<T>
	Interpolator<T>::scalarDerivativeWMLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
													 const ScalarField<T>& field,
													 uint64_t index,
													 uint32_t n)
	{
		//	Generate a monomial up to order n for the vanilla method
		Monomial mono(ugrid->getDim(),n);
		//	Query the nearest neighbors of ugrid
		ugrid->queryNeighbors(_params.k);
		//	Get the nearest neighbors for the point 'index'
		std::vector<uint64_t> neighbors = ugrid->getNeighbors(index);
		//	Construct the taylor matrix B
		Matrix<T> B = constructTaylorMatrix(ugrid,neighbors,index,mono);
		//	Construct the weight matrix W
		Matrix<T> W = constructWeightMatrix(ugrid,neighbors,index);
		//	Construct the vector of corresponding field values for neighboPointrs
		Vector<T> field_neighbors(_params.k,0.0);
		for (uint32_t i = 0; i < _params.k; i++)
		{
			field_neighbors(i) = field(neighbors[i]);
		}
		//	compute the transpose of B
		Matrix<T> B_T = B.transpose();
		//	compute the product (B^T*B)
		Matrix<T> B_TWB = B_T * W * B;
		//	compute the inverse of the product
		Matrix<T> B_TWB_inv = DGETRI(B_TWB);
		//	finally take the product
		Vector<T> coefficients = B_TWB_inv * B_T * W * field_neighbors;
		if (B.getFlag() == -1)
		{
			m_log->ERROR(B.getInfo());
			m_log->ERROR(B_TWB.getInfo());
		}
		return coefficients;
	}
	//----------------------------------------------------------------------------
	//  scalarDerivativeWMLSPoint - approximate the derivative for a point
	//                             using the WMLS method
	//  Arguments:  ugrid        - UGrid<T> pointer
	//              field        - const ScalarField<T>& reference
	//              point        - std::vector<T> of the point
	//              n            - order of the derivative
	//
	//  Returns:    std::vector<T> of the gradient.
	//----------------------------------------------------------------------------
	template<typename T>
	Vector<T>
	Interpolator<T>::scalarDerivativeWMLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
													 const ScalarField<T>& field,
													 std::vector<T> point,
													 uint32_t n)
	{
		//	Generate a monomial up to order n for the vanilla method
		Monomial mono(ugrid->getDim(),n);
		//	Get the nearest neighbors for the point 'index'
		std::vector<uint64_t> neighbors = ugrid->queryNeighbors(point,_params.k);
		//	Construct the taylor matrix B
		Matrix<T> B = constructTaylorMatrix(ugrid,neighbors,point,n);
		//	Construct the weight matrix W
		Matrix<T> W = constructWeightMatrix(ugrid,neighbors,point);
		//	Construct the vector of corresponding field values for neighboPointrs
		Vector<T> field_neighbors(_params.k,0.0);
		for (uint32_t i = 0; i < _params.k; i++)
		{
			field_neighbors(i) = field(neighbors[i]);
		}
		//	Complete the vanilla least squares to get the coefficients
		//	compute the transpose of B
		Matrix<T> B_T = B.transpose();
		//	compute the product (B^T*B)
		Matrix<T> B_TWB = B_T * W * B;
		//	compute the inverse of the product
		Matrix<T> B_TWB_inv = DGETRI(B_TWB);
		//	finally take the product
		Vector<T> coefficients = B_TWB_inv * B_T * W * field_neighbors;
		if (B.getFlag() == -1)
		{
			m_log->ERROR(B.getInfo());
			m_log->ERROR(B_TWB.getInfo());
		}
		return coefficients;
	}
	//----------------------------------------------------------------------------
	//  scalarDerivativeWMLS      - approximate the derivative for an entire
	//                             field using the weighted MLS method
	//  Arguments:  ugrid        - UGrid<T> pointer
	//              field        - const ScalarField<T>& reference
	//              n            - order of the derivative
	//
	//  Returns:    std::vector<T> of the gradient.
	//----------------------------------------------------------------------------
	template<typename T>
	std::vector<Vector<T>>
	Interpolator<T>::scalarDerivativeWMLS(const std::shared_ptr<UGrid<T>> ugrid,
											const ScalarField<T>& field,
											uint32_t n)
	{
		std::vector<Vector<T>> result(field.getN());
		//	Generate a monomial up to order n for the vanilla method
		Monomial mono(ugrid->getDim(),n);
		//	Query the nearest neighbors of ugrid
		ugrid->queryNeighbors(_params.k);
		for (uint64_t i = 0; i < field.getN(); i++)
		{
			//	Get the nearest neighbors for the point 'index'
			std::vector<uint64_t> neighbors = ugrid->getNeighbors(i);
			//	Construct the taylor matrix B
			Matrix<T> B = constructTaylorMatrix(ugrid,neighbors,i,mono);
			//	Construct the weight matrix W
			Matrix<T> W = constructWeightMatrix(ugrid,neighbors,i);
			//	Construct the vector of corres
			Vector<T> field_neighbors(_params.k,0.0);
			for (uint32_t i = 0; i < _params.k; i++)
			{
				field_neighbors(i) = field(neighbors[i]);
			}
			//	Complete the vanilla least squares to get the coefficients
			//	compute the transpose of B
			Matrix<T> B_T = B.transpose();
			//	compute the product (B^T*B)
			Matrix<T> B_TWB = B_T * W * B;
			//	compute the inverse of the product
			Matrix<T> B_TWB_inv = DGETRI(B_TWB);
			//	finally take the product
			Vector<T> coefficients = B_TWB_inv * B_T * W * field_neighbors;
			if (B.getFlag() == -1)
			{
				m_log->ERROR(B.getInfo());
				m_log->ERROR(B_TWB.getInfo());
			}
			result[i] = coefficients;
		}
		return result;
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  Helper functions
	//----------------------------------------------------------------------------
	//----------------------------------------------------------------------------
	template<typename T>
	std::vector<Vector<T>>
	Interpolator<T>::xScalarDerivative(const std::shared_ptr<UGrid<T>> ugrid,
																	 	 const std::shared_ptr<ScalarField<T>> field,
																	 	 uint32_t n)
	{
		if (_type == InterpolatorType::MLS)
		{
			return scalarDerivativeMLS(ugrid,field,n);
		}
		else if (_type == InterpolatorType::WMLS)
		{
			return scalarDerivativeWMLS(ugrid,field,n);
		}
		else
		{
			return scalarDerivativeLS(ugrid,field,n);
		}
	}
	//----------------------------------------------------------------------------
	template<typename T>
	std::vector<Vector<T>>
	Interpolator<T>::xScalarDerivative(const std::shared_ptr<UGrid<T>> ugrid,
																	 	 const ScalarField<T>& field,
																	 	 uint32_t n)
	{
		if (_type == InterpolatorType::MLS)
		{
			return scalarDerivativeMLS(ugrid,field,n);
		}
		else if (_type == InterpolatorType::WMLS)
		{
			return scalarDerivativeWMLS(ugrid,field,n);
		}
		else
		{
			return scalarDerivativeLS(ugrid,field,n);
		}
	}
	//----------------------------------------------------------------------------
	template<typename T>
	Vector<T>
	Interpolator<T>::xScalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
																	 const std::shared_ptr<ScalarField<T>> field,
																	 uint64_t index, uint32_t n)
	{
		if (_type == InterpolatorType::MLS)
		{
			return scalarDerivativeMLSPoint(ugrid,field,index,n);
		}
		else if (_type == InterpolatorType::WMLS)
		{
			return scalarDerivativeWMLSPoint(ugrid,field,index,n);
		}
		else
		{
			return scalarDerivativeLSPoint(ugrid,field,index,n);
		}
	}
	//----------------------------------------------------------------------------
	template<typename T>
	Vector<T>
	Interpolator<T>::xScalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
																	 const ScalarField<T>& field,
											             uint64_t index, uint32_t n)
	{
		if (_type == InterpolatorType::MLS)
		{
			return scalarDerivativeMLSPoint(ugrid,field,index,n);
		}
		else if (_type == InterpolatorType::WMLS)
		{
			return scalarDerivativeWMLSPoint(ugrid,field,index,n);
		}
		else
		{
			return scalarDerivativeLSPoint(ugrid,field,index,n);
		}
	}
	//----------------------------------------------------------------------------
	template<typename T>
	Vector<T>
	Interpolator<T>::xScalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
																	 const std::shared_ptr<ScalarField<T>> field,
																	 std::vector<T> point, uint32_t n)
	{
		if (_type == InterpolatorType::MLS)
		{
			return scalarDerivativeMLSPoint(ugrid,field,point,n);
		}
		else if (_type == InterpolatorType::WMLS)
		{
			return scalarDerivativeWMLSPoint(ugrid,field,point,n);
		}
		else
		{
			return scalarDerivativeLSPoint(ugrid,field,point,n);
		}
	}
	//----------------------------------------------------------------------------
	template<typename T>
	Vector<T>
	Interpolator<T>::xScalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
																	 const ScalarField<T>& field,
																	 std::vector<T> point, uint32_t n)
	{
		if (_type == InterpolatorType::MLS)
		{
			return scalarDerivativeMLSPoint(ugrid,field,point,n);
		}
		else if (_type == InterpolatorType::WMLS)
		{
			return scalarDerivativeWMLSPoint(ugrid,field,point,n);
		}
		else
		{
			return scalarDerivativeLSPoint(ugrid,field,point,n);
		}
	}
	//----------------------------------------------------------------------------
	template<typename T>
	Vector<T> Interpolator<T>::xGELSx(Matrix<T> B, Vector<T> u)
	{
		if (m_lsdriver == LSDriver::xGELS) {
			return DGELS(B,u);
		}
		else if (m_lsdriver == LSDriver::xGELSY) {
			return DGELSY(B,u);
		}
		else if (m_lsdriver == LSDriver::xGELSD) {
			return DGELSD(B,u);
		}
		else {
			return DGELSS(B,u);
		}
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Derivatives of vector fields
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//
	//----------------------------------------------------------------------------
	template<typename T>
	std::vector<std::vector<T>>
	Interpolator<T>::vectorDerivative(const std::shared_ptr<UGrid<T>> ugrid,
									                  const VectorField<T>& field,
									                  uint32_t dir, uint32_t n)
	{
		std::vector<std::vector<T>> result(field.getN());
		Monomial mono(ugrid->getDim(),n);
		ugrid->queryNeighbors(_params.k);
		for (uint64_t p = 0; p < field.getN(); p++)
		{
			std::vector<uint64_t> neighbors = ugrid->getNeighbors(p);
			Matrix<T> B = constructTaylorMatrix(ugrid,neighbors,p,mono);
			std::vector<T> field_neighbors_x(_params.k);
			std::vector<T> field_neighbors_y(_params.k);
			std::vector<T> field_neighbors_z(_params.k);
			for (uint32_t i = 0; i < _params.k; i++)
			{
				field_neighbors_x[i] = field(neighbors[i],0);
				field_neighbors_y[i] = field(neighbors[i],1);
				field_neighbors_z[i] = field(neighbors[i],2);
			}
			Vector<T> field_vals_x(field_neighbors_x);
			Vector<T> answer_x = xGELSx(B,field_vals_x);
			if (B.getFlag() == -1)
			{
				m_log->ERROR(B.getInfo());
			}
			Vector<T> field_vals_y(field_neighbors_y);
			Vector<T> answer_y = xGELSx(B,field_vals_y);
			if (B.getFlag() == -1)
			{
				m_log->ERROR(B.getInfo());
			}
			Vector<T> field_vals_z(field_neighbors_z);
			Vector<T> answer_z = xGELSx(B,field_vals_z);
			if (B.getFlag() == -1)
			{
				m_log->ERROR(B.getInfo());
			}
			//	Get the index of the monomial expansion corresponding to
			//	the nth-derivative in the 'dir'-direction
			std::vector<uint32_t> deriv(field.getDim(),0);
			deriv[dir] = n;
			uint32_t index = mono.getTaylorIndex(deriv);
			result[p][0] = answer_x(index);
			result[p][1] = answer_y(index);
			result[p][2] = answer_z(index);
		}
		return result;
	}
	//----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Construct the Taylor expansion matrix
  //----------------------------------------------------------------------------
  template<typename T>
  Matrix<T>
	Interpolator<T>::constructTaylorMatrix(				 const std::shared_ptr<UGrid<T>> ugrid,
    const std::vector<uint64_t> neighbors, uint64_t index, uint64_t order)
  {
    Monomial mono(ugrid->getDim(),order);
    std::vector<std::vector<double>> B;
    for (auto i = 0; i < neighbors.size(); i++)
    {
      auto id = neighbors[i];
      std::vector<double>
			temp = mono.taylorMonomialExpansion(ugrid->getPoint(index),
                                          ugrid->getPoint(id));
      B.push_back(temp);
    }
    return Matrix<T>("B",B);
  }
  //----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
  //  Construct the Taylor expansion matrix about a point
  //------------------------------------				 ----------------------------------------
  template<typename T>
  Matrix<T>
	Interpolator<T>::constructTaylorMatrix(const std::shared_ptr<UGrid<T>> ugrid,
    const std::vector<uint64_t> neighbors, std::vector<T> point, uint64_t order)
  {
    Monomial mono(ugrid->getDim(),order);
    std::vector<std::vector<double>> B;
    for (uint64_t i = 0; i < neighbors.size(); i++)
    {
      uint64_t id = neighbors[i];
      std::vector<double>
			temp = mono.taylorMonomialExpansion(point,
                                          ugrid->getPoint(id));
      B.push_back(temp);
    }
    return Matrix<T>("B",B);
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Construct the Taylor expansion matrix
  //  with mono
  //----------------------------------------------------------------------------
  template<typename T>
  Matrix<T>
	Interpolator<T>::constructTaylorMatrix(const std::shared_ptr<UGrid<T>> ugrid,
    const std::vector<uint64_t> neighbors, uint64_t index, Monomial& mono)
  {
    std::vector<std::vector<double>> B;
    for (uint64_t i = 0; i < neighbors.size(); i++)
    {
      uint64_t id = neighbors[i];
      std::vector<double>
			temp = mono.taylorMonomialExpansion(ugrid->getPoint(index),
                                          ugrid->getPoint(id));
      B.push_back(temp);
    }
    return Matrix<T>("B",B);
  }
  //----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Construct weight matrix
	//----------------------------------------------------------------------------
	template<typename T>
	Matrix<T>
	Interpolator<T>::constructWeightMatrix(const std::shared_ptr<UGrid<T>> ugrid,
    const std::vector<uint64_t> neighbors, uint64_t index)
	{
	  //	For a simple gaussian weight, find the distances from index
		//	to each point in neighbors.
		std::vector<double> distances = ugrid->getDistances()[index];
		Matrix<T> W(neighbors.size());
		for (uint32_t i = 0; i < neighbors.size(); i++)
		{
			W(i,i) = exp(-.5*distances[i]);
		}
		return W;
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Construct weight matrix
	//----------------------------------------------------------------------------
	template<typename T>
	Matrix<T>
	Interpolator<T>::constructWeightMatrix(const std::shared_ptr<UGrid<T>> ugrid,
    const std::vector<uint64_t> neighbors, std::vector<T> point)
	{
	  //	For a simple gaussian weight, find the distances from point
		//	to each point in neighbors.
		std::vector<double> distances = ugrid->queryDistances(point,_params.k);
		Matrix<T> W(neighbors.size());
		for (uint32_t i = 0; i < neighbors.size(); i++)
		{
			W(i,i) = exp(-.5*distances[i]);
		}
		return W;
	}
	//----------------------------------------------------------------------------

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
    sum += "\n       name:  '" + mm_name + "'";
    sum += "\n       type:  '" + InterpolatorTypeNameMap[_type] + "'";
    sum += "\n   LSDriver:  '" + LSDriverNameMap[m_lsdriver] + "'";
		sum += "\n          k:   " + std::to_string(_params.k);
    sum += "\n          n:   " + std::to_string(_params.n);
		sum += "\n---------------------------------------------------";
    //sum += "\n   RBF at: " + getMem(*getRadialBasisInterpolator()) + ",";
		//sum += "\n   ref at: " + getMem(m_rbf);
		sum += "\nLogger at: " + getMem(*getLogger()) + ",";
		sum += "\n   ref at: " + getMem(m_log);
		sum += "\n++++++++++++++++++++++++++++++++++++++++++++++++++++";
		return sum;
  }
  //----------------------------------------------------------------------------
}
