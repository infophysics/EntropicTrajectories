//------------------------------------------------------------------------------
//  approximator.cpp
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
#include "approximator.h"


namespace ET
{
  //----------------------------------------------------------------------------
  //  Enum maps
  //----------------------------------------------------------------------------
  //  map for a string to enum of approximator type
  std::map<std::string, ApproxType> ApproxTypeMap =
  {
    { "LS", ApproxType::LS },
		{ "MLS", ApproxType::MLS },
		{ "WMLS", ApproxType::WMLS},
    { "RBF", ApproxType::RBF }
  };
  //  map for a string to enum of LSDriver type
  std::map<std::string, LSDriver> LSDriverMap =
  {
    { "xGELS",  LSDriver::xGELS },
    { "xGELSY", LSDriver::xGELSY },
    { "xGELSD", LSDriver::xGELSD },
    { "xGELSS", LSDriver::xGELSS }
  };
  std::map<ApproxType, std::string> ApproxTypeNameMap =
  {
    { ApproxType::LS, "Vanilla least squares" },
		{ ApproxType::MLS, "Moving least squares" },
		{ ApproxType::WMLS, "Weighted moving least squares"},
    { ApproxType::RBF, "Radial basis functions" }
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
  Approximator<T>::Approximator() : _name("default")
  {
    _type = 0;
    _lsdriver = 0;
    //##########################################################################
		_log = std::make_shared<Log>();
		_log->init("ET:Approximator:default", ".logs/approx_default.txt");
		_log->TRACE("Approximator 'default' created at location "
		            + getMem(*this));
		//##########################################################################
  }
	//----------------------------------------------------------------------------
  //  Default destructor
  //----------------------------------------------------------------------------
  template<typename T>
  Approximator<T>::~Approximator()
  {
		//##########################################################################
		_log->TRACE("Approximator '" + _name
								+ "' destroyed at location " + getMem(*this));
		//##########################################################################
	}
	//----------------------------------------------------------------------------
  //  Various constructors taking in arguments for
	//
  //----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  Constructor with approximator type
	//----------------------------------------------------------------------------
  template<typename T>
  Approximator<T>::Approximator(int type) : _type(type), _name("default")
  {
		//##########################################################################
		_log = std::make_shared<Log>();
		_log->init("ET:Approximator:default", ".logs/approx_default.txt");
		_log->TRACE("Approximator 'default' created at location "
								+ getMem(*this));
		//##########################################################################
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
  //  Constructor taking in approximator type as a string
  //----------------------------------------------------------------------------
  template<typename T>
  Approximator<T>::Approximator(std::string type) : _name("default")
  {
    _type = ApproxTypeMap[type];
		//##########################################################################
		_log = std::make_shared<Log>();
		_log->init("ET:Approximator:default", ".logs/approx_default.txt");
		_log->TRACE("Approximator 'default' created at location "
		            + getMem(*this));
		//##########################################################################
  }
  //----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Various constructors taking in shared logger
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Constructor with shared logger
	//----------------------------------------------------------------------------
	template<typename T>
  Approximator<T>::Approximator(std::shared_ptr<Log> log) : _name("default")
  {
    _type = 0;
    _lsdriver = 0;
    //##########################################################################
		_log = log;
		_log->TRACE("Approximator 'default' created at location "
		            + getMem(*this));
		_log->INFO("Logger passed to Approximator 'default'");
		//##########################################################################
  }
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  Constructor with approximator type
	//----------------------------------------------------------------------------
  template<typename T>
  Approximator<T>::Approximator(int type, std::shared_ptr<Log> log)
	: _type(type), _name("default")
  {
		//##########################################################################
		_log = log;
		_log->TRACE("Approximator 'default' created at location "
								+ getMem(*this));
		_log->INFO("Logger passed to Approximator 'default'");
		//##########################################################################
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
  //  Constructor taking in approximator type as a string
  //----------------------------------------------------------------------------
  template<typename T>
  Approximator<T>::Approximator(std::string type, std::shared_ptr<Log> log)
	: _name("default")
  {
    _type = ApproxTypeMap[type];
		//##########################################################################
		_log = log;
		_log->TRACE("Approximator 'default' created at location "
		            + getMem(*this));
		_log->INFO("Logger passed to Approximator 'default'");
		//##########################################################################
  }
  //----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
  //  Getters and Setters
  //----------------------------------------------------------------------------
  template<typename T>
  int Approximator<T>::getApproxType() const
  {
    return _type;
  }
  template<typename T>
  ApproxParams Approximator<T>::getApproxParams() const
  {
    return _params;
  }
  template<typename T>
  int Approximator<T>::getLSDriver() const
  {
    return _lsdriver;
  }
  template<typename T>
  int Approximator<T>::getFlag() const
  {
    return _flag;
  }
  template<typename T>
  std::string Approximator<T>::getInfo() const
  {
    return _info;
  }
	template<typename T>
	std::shared_ptr<Log> Approximator<T>::getLogger()
	{
		return _log;
	}
  template<typename T>
  void Approximator<T>::setApproxType(std::string type)
  {
    auto res = ApproxTypeMap.find(type);
		if (res == ApproxTypeMap.end())
		{
			//########################################################################
			_log->ERROR("Approximator " + _name + ": Attempted to set ApproxType "
		              + "to " + type + " which is not a valid type");
			return;
			//########################################################################
		}
		else
		{
			_type = ApproxTypeMap[type];
			//########################################################################
			_log->INFO("Approximator " + _name + ": ApproxType set to " + type);
			//########################################################################
		}
  }
  template<typename T>
  void Approximator<T>::setApproxParams(ApproxParams params)
  {
    _params = params;
		//##########################################################################
		_log->INFO("Approximator " + _name + ": ApproxParams set");
		//##########################################################################
  }
  template<typename T>
  void Approximator<T>::setLSDriver(std::string type)
  {
		auto res = LSDriverMap.find(type);
		if (res == LSDriverMap.end())
		{
			//########################################################################
			_log->ERROR("Approximator " + _name + ": Attempted to set LSDriver "
		              + "to " + type + " which is not a valid type");
			return;
			//########################################################################
		}
		else
		{
			_lsdriver = LSDriverMap[type];
			//########################################################################
			_log->INFO("Approximator " + _name + ": LSDriver set to " + type);
			//########################################################################
		}
	}
  template<typename T>
  void Approximator<T>::set_k(uint64_t k)
  {
    _params.k = k;
		//##########################################################################
		_log->INFO("Approximator " + _name + ": k set to " + std::to_string(k));
		//##########################################################################
  }
  template<typename T>
  void Approximator<T>::set_n(uint64_t n)
  {
    _params.n = n;
		//##########################################################################
		_log->INFO("Approximator " + _name + ": n set to " + std::to_string(n));
		//##########################################################################
  }
  template<typename T>
  void Approximator<T>::setFlag(int flag)
  {
    _flag = flag;
  }
  template<typename T>
  void Approximator<T>::setInfo(std::string info)
  {
    _info = info;
  }
  //----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  Gradient functions
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  Scalar fields
	//----------------------------------------------------------------------------
	//----------------------------------------------------------------------------
	//  scalarGradientPoint - approximate the gradient at a point
	//  Arguments:  ugrid   - UGrid<T> pointer
	//              field   - ScalarField<T> pointer
	//              index   - index of the point
	//
	//  Returns:    std::vector<T> of the gradient.
	//----------------------------------------------------------------------------
  template<typename T>
  std::vector<T>
	Approximator<T>::scalarGradientPoint(const std::shared_ptr<UGrid<T>> ugrid,
    const std::shared_ptr<ScalarField<T>> field, uint64_t index)
  {
    if (_type == 0)
      return scalarGradientLSPoint(ugrid, field, index);
  }
  //----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  scalarGradientLSPoint - approximate the gradient at a point
	//                           using the vanilla LS method
	//  Arguments:  ugrid      - UGrid<T> pointer
	//              field      - ScalarField<T> pointer
	//              index      - index of the point
	//
	//  Returns:    std::vector<T> of the gradient.
	//----------------------------------------------------------------------------
  template<typename T>
  std::vector<T>
	Approximator<T>::scalarGradientLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
    const std::shared_ptr<ScalarField<T>> field, uint64_t index)
  {
    std::vector<T> result(field->getDim());
		Monomial mono(ugrid->getDim(),_params.n);
    //  First, find the nearest neighbors associated to the point specified by
    //  index.
    ugrid->queryNeighbors(_params.k);
    std::vector<uint64_t> neighbors = ugrid->getNeighbors(index);
    //  Construct the matrix associated with the ugrid spacing
    Matrix<T> B = constructTaylorMatrix(ugrid,neighbors,index,mono);
    //  Construct the vector of field values associated to each point
    std::vector<T> field_neighbors(_params.k);
    for (uint32_t i = 0; i < _params.k; i++)
    {
      field_neighbors[i] = (*field)(neighbors[i]);
    }
    Vector<T> field_vals(field_neighbors);
    Vector<T> answer = xGELSx(B,field_vals);
		if (B.getFlag() == -1)
		{
			_log->ERROR(B.getInfo());
		}
    return answer.getVec();
  }
  //----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  scalarGradient      - approximate the gradient for an entire field
	//  Arguments:  ugrid   - UGrid<T> pointer
	//              field   - ScalarField<T> pointer
	//
	//  Returns:    std::vector<std::vector<T>> of the gradient.
	//----------------------------------------------------------------------------
  template<typename T>
  std::vector<std::vector<T>>
	Approximator<T>::scalarGradient(const std::shared_ptr<UGrid<T>> ugrid,
                                  const std::shared_ptr<ScalarField<T>> field)
  {
    if (_type == 0)
      return scalarGradientLS(ugrid, field);
  }
  //----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  scalarGradientLS      - approximate the gradient for an entire field
	//                           using the vanilla LS method
	//  Arguments:  ugrid      - UGrid<T> pointer
	//              field      - ScalarField<T> pointer
	//
	//  Returns:    std::vector<std::vector<T>> of the gradient.
	//----------------------------------------------------------------------------
  template<typename T>
  std::vector<std::vector<T>>
	Approximator<T>::scalarGradientLS(const std::shared_ptr<UGrid<T>> ugrid,
    													const std::shared_ptr<ScalarField<T>> field)
  {
    std::vector<std::vector<T>> result(field->getN());
    Monomial mono(ugrid->getDim(),_params.n);
		ugrid->queryNeighbors(_params.k);
    for (uint64_t i = 0; i < field->getN(); i++)
    {
      std::vector<uint64_t> neighbors = ugrid->getNeighbors(i);
      Matrix<T> B = constructTaylorMatrix(ugrid,neighbors,i,mono);
      std::vector<T> field_neighbors(_params.k);
      for (uint32_t i = 0; i < _params.k; i++)
      {
        field_neighbors[i] = (*field)(neighbors[i]);
      }
      Vector<T> field_vals(field_neighbors);
      Vector<T> answer = xGELSx(B,field_vals);
			if (B.getFlag() == -1)
			{
				_log->ERROR(B.getInfo());
			}
      //  Trim result to the first field->getDim() elements
			std::vector<T> v = answer.getVec();
      std::vector<T> u = {v.begin()+1,v.begin()+field->getDim()+1};
      result[i] = u;
    }
		return result;
  }
  //----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
  //  Passing field as a const reference
  //----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  scalarGradientPoint - approximate the gradient at a point
	//  Arguments:  ugrid   - UGrid<T> pointer
	//              field   - const ScalarField<T>& reference
	//              index   - index of the point
	//
	//  Returns:    std::vector<T> of the gradient.
	//----------------------------------------------------------------------------
  template<typename T>
  std::vector<T>
	Approximator<T>::scalarGradientPoint(const std::shared_ptr<UGrid<T>> ugrid,
    const ScalarField<T>& field, uint64_t index)
  {
    if (_type == 0)
      return scalarGradientLSPoint(ugrid, field, index);
  }
  //----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  scalarGradientLSPoint - approximate the gradient at a point
	//                           using the vanilla LS method
	//  Arguments:  ugrid      - UGrid<T> pointer
	//              field      - const ScalarField<T>& reference
	//              index      - index of the point
	//
	//  Returns:    std::vector<T> of the gradient.
	//----------------------------------------------------------------------------
  template<typename T>
  std::vector<T>
	Approximator<T>::scalarGradientLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
    const ScalarField<T>& field, uint64_t index)
  {
    std::vector<T> result(field.getDim());
		Monomial mono(ugrid->getDim(),_params.n);
    //  First, find the nearest neighbors associated to the point specified by
    //  index.
    ugrid->queryNeighbors(_params.k);
    std::vector<uint64_t> neighbors = ugrid->getNeighbors(index);
    //  Construct the matrix associated with the ugrid spacing
    Matrix<T> B = constructTaylorMatrix(ugrid,neighbors,index,mono);
    //  Construct the vector of field values associated to each point
    std::vector<T> field_neighbors(_params.k);
    for (uint32_t i = 0; i < _params.k; i++)
    {
      field_neighbors[i] = field(neighbors[i]);
    }
    Vector<T> field_vals(field_neighbors);
    Vector<T> answer = xGELSx(B,field_vals);
		if (B.getFlag() == -1)
		{
			_log->ERROR(B.getInfo());
		}
    return answer.getVec();
  }
  //----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  scalarGradient      - approximate the gradient for an entire field
	//  Arguments:  ugrid   - UGrid<T> pointer
	//              field   - const ScalarField<T>& reference
	//              index   - index of the point
	//
	//  Returns:    std::vector<std::vector<T>> of the gradient.
	//----------------------------------------------------------------------------
  template<typename T>
  std::vector<std::vector<T>>
	Approximator<T>::scalarGradient(const std::shared_ptr<UGrid<T>> ugrid,
                                  const ScalarField<T>& field)
  {
    if (_type == 0)
      return scalarGradientLS(ugrid, field);
  }
  //----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  scalarGradientLS      - approximate the gradient for an entire field
	//                           using the vanilla LS method
	//  Arguments:  ugrid      - UGrid<T> pointer
	//              field      - const ScalarField<T>& reference
	//              index      - index of the point
	//
	//  Returns:    std::vector<std::vector<T>> of the gradient.
	//----------------------------------------------------------------------------
  template<typename T>
  std::vector<std::vector<T>>
	Approximator<T>::scalarGradientLS(const std::shared_ptr<UGrid<T>> ugrid,
    													const ScalarField<T>& field)
  {
    std::vector<std::vector<T>> result(field.getN());
    Monomial mono(ugrid->getDim(),_params.n);
		ugrid->queryNeighbors(_params.k);
    for (uint64_t i = 0; i < field.getN(); i++)
    {
      std::vector<uint64_t> neighbors = ugrid->getNeighbors(i);
      Matrix<T> B = constructTaylorMatrix(ugrid,neighbors,i,mono);
      std::vector<T> field_neighbors(_params.k);
      for (uint32_t i = 0; i < _params.k; i++)
      {
        field_neighbors[i] = field(neighbors[i]);
      }
      Vector<T> field_vals(field_neighbors);
      Vector<T> answer = xGELSx(B,field_vals);
			if (B.getFlag() == -1)
			{
				_log->ERROR(B.getInfo());
			}
      //  Trim result to the first field->getDim() elements
			std::vector<T> v = answer.getVec();
      std::vector<T> u = {v.begin()+1,v.begin()+field.getDim()+1};
      result[i] = u;
    }
		return result;
  }
  //----------------------------------------------------------------------------

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
	Approximator<T>::scalarDerivative(const std::shared_ptr<UGrid<T>> ugrid,
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
	Approximator<T>::scalarDerivative(const std::shared_ptr<UGrid<T>> ugrid,
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
	Approximator<T>::scalarDerivative(const std::shared_ptr<UGrid<T>> ugrid,
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
	Approximator<T>::scalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
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
	T	Approximator<T>::scalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
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
	//              index      - index of the point
	//              deriv      - vector denoting the direction and order
	//
	//  Returns:    T          - the gradient in direction dir and order n.
	//----------------------------------------------------------------------------
	template<typename T>
	T	Approximator<T>::scalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
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
	Approximator<T>::scalarDerivative(const std::shared_ptr<UGrid<T>> ugrid,
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
	Approximator<T>::scalarDerivative(const std::shared_ptr<UGrid<T>> ugrid,
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
	Approximator<T>::scalarDerivative(const std::shared_ptr<UGrid<T>> ugrid,
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
	Approximator<T>::scalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
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
	T Approximator<T>::scalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
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
	//              field      - const ScalarField<T>& reference
	//              index      - index of the point
	//              deriv      - vector denoting the direction and order
	//
	//  Returns:    T          - the gradient in direction dir and order n.
	//----------------------------------------------------------------------------
	template<typename T>
	T Approximator<T>::scalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
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
	Approximator<T>::scalarDerivativeLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
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
			_log->ERROR(B.getInfo());
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
	Approximator<T>::scalarDerivativeLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
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
			_log->ERROR(B.getInfo());
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
	Approximator<T>::scalarDerivativeLS(const std::shared_ptr<UGrid<T>> ugrid,
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
				_log->ERROR(B.getInfo());
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
	Approximator<T>::scalarDerivativeLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
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
			_log->ERROR(B.getInfo());
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
	Approximator<T>::scalarDerivativeLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
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
			_log->ERROR(B.getInfo());
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
	Approximator<T>::scalarDerivativeLS(const std::shared_ptr<UGrid<T>> ugrid,
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
				_log->ERROR(B.getInfo());
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
	Approximator<T>::scalarDerivativeMLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
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
			_log->ERROR(B.getInfo());
			_log->ERROR(B_TB.getInfo());
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
	Approximator<T>::scalarDerivativeMLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
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
			_log->ERROR(B.getInfo());
			_log->ERROR(B_TB.getInfo());
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
	Approximator<T>::scalarDerivativeMLS(const std::shared_ptr<UGrid<T>> ugrid,
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
				_log->ERROR(B.getInfo());
				_log->ERROR(B_TB.getInfo());
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
	Approximator<T>::scalarDerivativeMLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
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
			_log->ERROR(B.getInfo());
			_log->ERROR(B_TB.getInfo());
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
	Approximator<T>::scalarDerivativeMLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
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
			_log->ERROR(B.getInfo());
			_log->ERROR(B_TB.getInfo());
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
	Approximator<T>::scalarDerivativeMLS(const std::shared_ptr<UGrid<T>> ugrid,
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
				_log->ERROR(B.getInfo());
				_log->ERROR(B_TB.getInfo());
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
	Approximator<T>::scalarDerivativeWMLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
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
			_log->ERROR(B.getInfo());
			_log->ERROR(B_TWB.getInfo());
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
	Approximator<T>::scalarDerivativeWMLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
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
			_log->ERROR(B.getInfo());
			_log->ERROR(B_TWB.getInfo());
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
	Approximator<T>::scalarDerivativeWMLS(const std::shared_ptr<UGrid<T>> ugrid,
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
				_log->ERROR(B.getInfo());
				_log->ERROR(B_TWB.getInfo());
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
	Approximator<T>::scalarDerivativeWMLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
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
			_log->ERROR(B.getInfo());
			_log->ERROR(B_TWB.getInfo());
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
	Approximator<T>::scalarDerivativeWMLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
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
			_log->ERROR(B.getInfo());
			_log->ERROR(B_TWB.getInfo());
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
	Approximator<T>::scalarDerivativeWMLS(const std::shared_ptr<UGrid<T>> ugrid,
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
				_log->ERROR(B.getInfo());
				_log->ERROR(B_TWB.getInfo());
			}
			result[i] = coefficients;
		}
		return result;
	}
	//----------------------------------------------------------------------------
	//----------------------------------------------------------------------------
	//  Radial basis interpolation
	//----------------------------------------------------------------------------
	//----------------------------------------------------------------------------
	//  scalarDerivativeRBFPoint - approximate the derivative for a point
	//                             using the RBF method
	//  Arguments:  ugrid        - UGrid<T> pointer
	//              field        - ScalarField<T> pointer
	//              index        - uint64_t
	//              n            - order of the derivative
	//
	//  Returns:    std::vector<T> of the gradient.
	//----------------------------------------------------------------------------
	template<typename T>
	Vector<T>
	Approximator<T>::scalarDerivativeRBFPoint(const std::shared_ptr<UGrid<T>> ugrid,
													 const std::shared_ptr<ScalarField<T>> field,
													 uint64_t index,
													 uint32_t n)
	{
		//	Query the nearest neighbors of ugrid
		ugrid->queryNeighbors(_params.k);
		//	Get the nearest neighbors for the point 'index'
		std::vector<uint64_t> neighbors = ugrid->getNeighbors(index);
		//	Construct the RBF matrix A
		Matrix<T> A = constructRBFMatrix(ugrid,neighbors,index);
		//	Construct the vector of corresponding field values for neighboPointrs
		Vector<T> field_neighbors(_params.k,0.0);
		for (uint32_t i = 0; i < _params.k; i++)
		{
			field_neighbors(i) = (*field)(neighbors[i]);
		}
		//	Complete the weights using least squares
		//	finally take the product
		Vector<T> weights = DGELS(A,field_neighbors);
		if (A.getFlag() == -1)
		{
			_log->ERROR(A.getInfo());
		}
		return weights;
	}
	//----------------------------------------------------------------------------
	//  Helper functions
	//----------------------------------------------------------------------------
	//----------------------------------------------------------------------------
	template<typename T>
	std::vector<Vector<T>>
	Approximator<T>::xScalarDerivative(const std::shared_ptr<UGrid<T>> ugrid,
																	 	 const std::shared_ptr<ScalarField<T>> field,
																	 	 uint32_t n)
	{
		if (_type == 1)
		{
			return scalarDerivativeMLS(ugrid,field,n);
		}
		else if (_type == 2)
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
	Approximator<T>::xScalarDerivative(const std::shared_ptr<UGrid<T>> ugrid,
																	 	 const ScalarField<T>& field,
																	 	 uint32_t n)
	{
		if (_type == 1)
		{
			return scalarDerivativeMLS(ugrid,field,n);
		}
		else if (_type == 2)
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
	Approximator<T>::xScalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
																	 const std::shared_ptr<ScalarField<T>> field,
																	 uint64_t index, uint32_t n)
	{
		if (_type == 1)
		{
			return scalarDerivativeMLSPoint(ugrid,field,index,n);
		}
		else if (_type == 2)
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
	Approximator<T>::xScalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
																	 const ScalarField<T>& field,
																	 uint64_t index, uint32_t n)
	{
		if (_type == 1)
		{
			return scalarDerivativeMLSPoint(ugrid,field,index,n);
		}
		else if (_type == 2)
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
	Approximator<T>::xScalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
																	 const std::shared_ptr<ScalarField<T>> field,
																	 std::vector<T> point, uint32_t n)
	{
		if (_type == 1)
		{
			return scalarDerivativeMLSPoint(ugrid,field,point,n);
		}
		else if (_type == 2)
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
	Approximator<T>::xScalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
																	 const ScalarField<T>& field,
																	 std::vector<T> point, uint32_t n)
	{
		if (_type == 1)
		{
			return scalarDerivativeMLSPoint(ugrid,field,point,n);
		}
		else if (_type == 2)
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
	Vector<T> Approximator<T>::xGELSx(Matrix<T> B, Vector<T> u)
	{
		if (_lsdriver == 0)
		{
			return DGELS(B,u);
		}
		else if (_lsdriver == 1)
		{
			return DGELSY(B,u);
		}
		else if (_lsdriver == 2)
		{
			return DGELSD(B,u);
		}
		else
		{
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
	Approximator<T>::vectorDerivative(const std::shared_ptr<UGrid<T>> ugrid,
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
				_log->ERROR(B.getInfo());
			}
			Vector<T> field_vals_y(field_neighbors_y);
			Vector<T> answer_y = xGELSx(B,field_vals_y);
			if (B.getFlag() == -1)
			{
				_log->ERROR(B.getInfo());
			}
			Vector<T> field_vals_z(field_neighbors_z);
			Vector<T> answer_z = xGELSx(B,field_vals_z);
			if (B.getFlag() == -1)
			{
				_log->ERROR(B.getInfo());
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
	Approximator<T>::constructTaylorMatrix(const std::shared_ptr<UGrid<T>> ugrid,
    const std::vector<uint64_t> neighbors, uint64_t index, uint64_t order)
  {
    Monomial mono(ugrid->getDim(),order);
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
  //  Construct the Taylor expansion matrix about a point
  //----------------------------------------------------------------------------
  template<typename T>
  Matrix<T>
	Approximator<T>::constructTaylorMatrix(const std::shared_ptr<UGrid<T>> ugrid,
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
	Approximator<T>::constructTaylorMatrix(const std::shared_ptr<UGrid<T>> ugrid,
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
	Approximator<T>::constructWeightMatrix(const std::shared_ptr<UGrid<T>> ugrid,
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
	Approximator<T>::constructWeightMatrix(const std::shared_ptr<UGrid<T>> ugrid,
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
  std::string Approximator<T>::summary()
  {
    std::string s = "\nApproximator type: " + ApproxTypeNameMap[_type];
    if (_type == LS)
    {
      s += "\nLeast squares driver type: " + LSDriverNameMap[_lsdriver];
    }
    s += "\nApproximator parameters - k = " + std::to_string(_params.k);
    s += "\n                          n = " + std::to_string(_params.n);
    return s;
  }
  //----------------------------------------------------------------------------
}
