//------------------------------------------------------------------------------
//  scalarfield.cpp
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
#include "scalarfield.h"


namespace ET
{
  //----------------------------------------------------------------------------
  //  ScalarField constructors
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Default constructor
  //    sets name = "default", and _dim, _N = 0
  //----------------------------------------------------------------------------
  template<typename T>
  ScalarField<T>::ScalarField() : _dim(0), _N(0), _name("default")
  {
    //##########################################################################
		_log = std::make_shared<Log>();
		_log->init("ET:ScalarField:default", ".logs/scalarfield_default.txt");
		_log->TRACE("Scalar Field 'default' created at location "
		            + getMem(*this));
		//##########################################################################
  }
	//----------------------------------------------------------------------------
  //  Default destructor
  //----------------------------------------------------------------------------
  template<typename T>
  ScalarField<T>::~ScalarField()
  {
		//##########################################################################
		_log->TRACE("Scalar Field '" + _name
								+ "' destroyed at location " + getMem(*this));
		//##########################################################################
  }
	//----------------------------------------------------------------------------
  //  Various constructors taking in arguments for
	//		_dim, _name, _N, _ugrid and _log.
  //----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Constructor with UGrid
	//----------------------------------------------------------------------------
  template<typename T>
  ScalarField<T>::ScalarField(std::shared_ptr<UGrid<T>> ugrid)
	: _ugrid(ugrid), _name("default")
  {
    _N = _ugrid->getN();
    _approx = std::make_shared<Approximator<T>>();
		//##########################################################################
		_log = std::make_shared<Log>();
		_log->init("ET:ScalarField:default", ".logs/scalarfield_default.txt");
		_log->TRACE("Scalar Field 'default' created at location "
		            + getMem(*this));
		//##########################################################################
  }
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  Constructor with name and UGrid
	//----------------------------------------------------------------------------
  template<typename T>
  ScalarField<T>::ScalarField(std::string name, std::shared_ptr<UGrid<T>> ugrid)
  : _name(name), _ugrid(ugrid)
  {
    _N = _ugrid->getN();
    _approx = std::make_shared<Approximator<T>>();
		//##########################################################################
		_log = std::make_shared<Log>();
		_log->init("ET:ScalarField:" + _name, ".logs/scalarfield_default.txt");
		_log->TRACE("Scalar Field '" + _name + "' created at location "
		            + getMem(*this));
		//##########################################################################
  }
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  Constructor with UGrid and field
	//----------------------------------------------------------------------------
  template<typename T>
  ScalarField<T>::ScalarField(std::shared_ptr<UGrid<T>> ugrid,
		                          std::vector<T> field)
  : _ugrid(ugrid), _field(field), _name("default")
  {
    _N = ugrid->getN();
    _approx = std::make_shared<Approximator<T>>();
		//##########################################################################
		_log = std::make_shared<Log>();
		_log->init("ET:ScalarField:default", ".logs/scalarfield_default.txt");
		_log->TRACE("Scalar Field 'default' created at location "
		            + getMem(*this));
		//##########################################################################
  }
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  Constructor with name, UGrid and field
	//----------------------------------------------------------------------------
  template<typename T>
  ScalarField<T>::ScalarField(std::string name, std::shared_ptr<UGrid<T>> ugrid,
                              std::vector<T> field)
  : _name(name), _ugrid(ugrid), _field(field)
  {
    _N = ugrid->getN();
    _approx = std::make_shared<Approximator<T>>();
		//##########################################################################
		_log = std::make_shared<Log>();
		_log->init("ET:ScalarField:" + _name, ".logs/scalarfield_default.txt");
		_log->TRACE("Scalar Field '" + _name + "' created at location "
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
	ScalarField<T>::ScalarField(std::shared_ptr<Log> log)
	: _name("default"), _dim(0), _N(0)
	{
		//##########################################################################
		_log = log;
		_log->TRACE("Scalar Field 'default' created at location "
								+ getMem(*this));
		_log->INFO("Logger passed to Scalar Field 'default'");
		//##########################################################################
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Constructor with UGrid and shared logger
	//----------------------------------------------------------------------------
	template<typename T>
	ScalarField<T>::ScalarField(std::shared_ptr<UGrid<T>> ugrid,
		                          std::shared_ptr<Log> log)
	: _ugrid(ugrid), _name("default")
	{
		_N = _ugrid->getN();
		_approx = std::make_shared<Approximator<T>>();
		//##########################################################################
		_log = log;
		_log->TRACE("Scalar Field 'default' created at location "
								+ getMem(*this));
		_log->INFO("Logger passed to Scalar Field 'default'");
		//##########################################################################
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  Constructor with name and UGrid
	//----------------------------------------------------------------------------
	template<typename T>
	ScalarField<T>::ScalarField(std::string name, std::shared_ptr<UGrid<T>> ugrid,
	                            std::shared_ptr<Log> log)
	: _name(name), _ugrid(ugrid)
	{
		_N = _ugrid->getN();
		_approx = std::make_shared<Approximator<T>>();
		//##########################################################################
		_log = log;
		_log->TRACE("Scalar Field '" + _name + "' created at location "
								+ getMem(*this));
		_log->INFO("Logger passed to Scalar Field '" + _name + "'");
		//##########################################################################
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  Constructor with UGrid and field
	//----------------------------------------------------------------------------
	template<typename T>
	ScalarField<T>::ScalarField(std::shared_ptr<UGrid<T>> ugrid,
															std::vector<T> field, std::shared_ptr<Log> log)
	: _ugrid(ugrid), _field(field), _name("default")
	{
		_N = ugrid->getN();
		_approx = std::make_shared<Approximator<T>>();
		//##########################################################################
		_log = log;
		_log->TRACE("Scalar Field 'default' created at location "
								+ getMem(*this));
		_log->INFO("Logger passed to Scalar Field 'default'");
		//##########################################################################
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  Constructor with name, UGrid and field
	//----------------------------------------------------------------------------
	template<typename T>
	ScalarField<T>::ScalarField(std::string name, std::shared_ptr<UGrid<T>> ugrid,
															std::vector<T> field, std::shared_ptr<Log> log)
	: _name(name), _ugrid(ugrid), _field(field)
	{
		_N = ugrid->getN();
		_approx = std::make_shared<Approximator<T>>();
		//##########################################################################
		_log = log;
		_log->TRACE("Scalar Field '" + _name + "' created at location "
								+ getMem(*this));
		_log->INFO("Logger passed to Scalar Field '" + _name + "'");
		//##########################################################################
	}
	//----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Getters and Setters
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<T> ScalarField<T>::getField() const
  {
    return _field;
  }
  template<typename T>
  std::vector<T>* ScalarField<T>::accessField()
  {
    return &_field;
  }
  template<typename T>
  T* ScalarField<T>::data()
  {
    return _field.data();
  }
  template<typename T>
  std::string ScalarField<T>::getName() const
  {
    return _name;
  }
  template<typename T>
  uint64_t ScalarField<T>::getN() const
  {
    return _N;
  }
  template<typename T>
  uint32_t ScalarField<T>::getDim() const
  {
    return _dim;
  }
  template<typename T>
  std::shared_ptr<Approximator<T>> ScalarField<T>::getApproximator() const
  {
    return _approx;
  }
  template<typename T>
  int ScalarField<T>::getFlag() const
  {
    return _flag;
  }
  template<typename T>
  std::string ScalarField<T>::getInfo() const
  {
    return _info;
  }
	template<typename T>
	std::shared_ptr<Log> ScalarField<T>::getLogger()
	{
		return _log;
	}
  template<typename T>
  void ScalarField<T>::setUGrid(std::shared_ptr<UGrid<T>> ugrid)
  {
    _ugrid = ugrid;
  }
  template<typename T>
  void ScalarField<T>::setField(std::vector<T> field)
  {
    _field = field;
  }
  template<typename T>
  void ScalarField<T>::setName(std::string name)
  {
    _name = name;
  }
  template<typename T>
  void ScalarField<T>::setApproxType(std::string type)
  {
    _approx->setApproxType(type);
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Operator overloads
  //----------------------------------------------------------------------------
  template<typename T>
  T& ScalarField<T>::operator()(const uint32_t& i)
  {
		if (i >= _N )
		{
			//########################################################################
			_log->ERROR("ScalarField " + _name
									+ ": Attempted to access _field array of size "
									+ std::to_string(_N) + " with index "
									+ std::to_string(i));
			//########################################################################
			if(_field.size() > 0)
			{
				//######################################################################
				_log->INFO("ScalarField "+ _name +": Returning the element at index 0");
				//######################################################################
				return _field[0];
			}
			else
			{
				//######################################################################
				_log->INFO("ScalarField " + _name + ": Terminating program");
				//######################################################################
				exit(0);
			}
		}
    return _field[i];
  }
  template<typename T>
  const T& ScalarField<T>::operator()(const uint32_t& i) const
  {
		if (i >= _N )
		{
			//########################################################################
			_log->ERROR("ScalarField " + _name
									+ ": Attempted to access _field array of size "
									+ std::to_string(_N) + " with index "
									+ std::to_string(i));
			//########################################################################
			if(_field.size() > 0)
			{
				//######################################################################
				_log->INFO("ScalarField "+ _name +": Returning the element at index 0");
				//######################################################################
				return _field[0];
			}
			else
			{
				//######################################################################
				_log->INFO("ScalarField " + _name + ": Terminating program");
				//######################################################################
				exit(0);
			}
		}
    return _field[i];
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Methods for calculating derivatives
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<std::vector<T>> ScalarField<T>::gradient()
  {
    return _approx->scalarGradient(_ugrid,std::make_shared<ScalarField<T>>(*this));
  }
  //----------------------------------------------------------------------------
}
