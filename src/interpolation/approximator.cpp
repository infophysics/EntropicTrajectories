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
    { "MLS", ApproxType::MLS },
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
    { ApproxType::MLS, "Moving least squares" },
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
  //  Methods for Scalar fields
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<T>
	Approximator<T>::scalarGradientPoint(const std::shared_ptr<UGrid<T>> ugrid,
    const std::shared_ptr<ScalarField<T>> field, uint64_t index)
  {
    if (_type == 0)
      return scalarGradientMLSPoint(ugrid, field, index);
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<T>
	Approximator<T>::scalarGradientMLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
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
    Vector<T> answer = DGELS(B,field_vals);
    return answer.getVec();
  }
  //----------------------------------------------------------------------------

  template<typename T>
  std::vector<std::vector<T>>
	Approximator<T>::scalarGradient(const std::shared_ptr<UGrid<T>> ugrid,
                                  const std::shared_ptr<ScalarField<T>> field)
  {
    if (_type == 0)
      return scalarGradientMLS(ugrid, field);
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<std::vector<T>>
	Approximator<T>::scalarGradientMLS(const std::shared_ptr<UGrid<T>> ugrid,
    													const std::shared_ptr<ScalarField<T>> field)
  {
    std::vector<std::vector<T>> result(field->getN());
    Monomial mono(ugrid->getDim(),_params.n);
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
      Vector<T> answer = DGELS(B,field_vals);
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
  template<typename T>
  std::vector<T>
	Approximator<T>::scalarGradientPoint(const std::shared_ptr<UGrid<T>> ugrid,
    const ScalarField<T>& field, uint64_t index)
  {
    if (_type == 0)
      return scalarGradientMLSPoint(ugrid, field, index);
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<T>
	Approximator<T>::scalarGradientMLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
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
    Vector<T> answer = DGELS(B,field_vals);
    return answer.getVec();
  }
  //----------------------------------------------------------------------------

  template<typename T>
  std::vector<std::vector<T>>
	Approximator<T>::scalarGradient(const std::shared_ptr<UGrid<T>> ugrid,
                                  const ScalarField<T>& field)
  {
    if (_type == 0)
      return scalarGradientMLS(ugrid, field);
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<std::vector<T>>
	Approximator<T>::scalarGradientMLS(const std::shared_ptr<UGrid<T>> ugrid,
    													const ScalarField<T>& field)
  {
    std::vector<std::vector<T>> result(field.getN());
    Monomial mono(ugrid->getDim(),_params.n);
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
      Vector<T> answer = DGELS(B,field_vals);
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
	template<typename T>
	std::vector<T>
	Approximator<T>::scalarDerivative(const std::shared_ptr<UGrid<T>> ugrid,
													        	const std::shared_ptr<ScalarField<T>> field,
														        uint32_t dir, uint32_t n)
	{
		std::vector<T> result(field->getN());
		Monomial mono(ugrid->getDim(),n);
		//	Get the index of the monomial expansion corresponding to
		//	the nth-derivative in the 'dir'-direction
		std::vector<uint32_t> deriv(field->getDim(),0);
		deriv[dir] = n;
		uint32_t index = mono.getTaylorIndex(deriv);
		//	Loop over every point
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
			Vector<T> answer = DGELS(B,field_vals);
			//  Trim result to the first field->getDim() elements
			result[i] = answer(index);
		}
		return result;
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  Passing field as a const reference
	//----------------------------------------------------------------------------
	template<typename T>
	std::vector<T>
	Approximator<T>::scalarDerivative(const std::shared_ptr<UGrid<T>> ugrid,
													        	const ScalarField<T>& field,
														        uint32_t dir, uint32_t n)
	{
		std::vector<T> result(field.getN());
		Monomial mono(ugrid->getDim(),n);
		//	Get the index of the monomial expansion corresponding to
		//	the nth-derivative in the 'dir'-direction
		std::vector<uint32_t> deriv(field.getDim(),0);
		deriv[dir] = n;
		uint32_t index = mono.getTaylorIndex(deriv);
		//	Loop over every point
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
			Vector<T> answer = DGELS(B,field_vals);
			//  Trim result to the first field->getDim() elements
			result[i] = answer(index);
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
  //  Various functions
  //----------------------------------------------------------------------------
  template<typename T>
  std::string Approximator<T>::summary()
  {
    std::string s = "\nApproximator type: " + ApproxTypeNameMap[_type];
    if (_type == MLS)
    {
      s += "\nLeast squares driver type: " + LSDriverNameMap[_lsdriver];
    }
    s += "\nApproximator parameters - k = " + std::to_string(_params.k);
    s += "\n                          n = " + std::to_string(_params.n);
    return s;
  }
  //----------------------------------------------------------------------------
}
