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
		_ugrid = std::make_shared<UGrid<T>>();
		_approx = std::make_shared<Approximator<T>>();
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
		_dim = _ugrid->getDim();
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
		_dim = _ugrid->getDim();
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
		_dim = _ugrid->getDim();
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
		_dim = _ugrid->getDim();
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
		_dim = _ugrid->getDim();
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
		_dim = _ugrid->getDim();
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
		_dim = _ugrid->getDim();
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
		_dim = _ugrid->getDim();
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
	std::shared_ptr<UGrid<T>> ScalarField<T>::getUGrid() const
	{
		return _ugrid;
	}
  template<typename T>
	std::shared_ptr<Log> ScalarField<T>::getLogger()
	{
		return _log;
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
	ScalarField<T> ScalarField<T>::operator+(const ScalarField<T>& scalar) const
	{
		if (_N != scalar.getN())
		{
			//########################################################################
			_log->ERROR("Attempted to add scalar fields " + _name  + " and "
			            + scalar.getName() + " with sizes " + std::to_string(_N)
									+ " and " + std::to_string(scalar.getN()));
			//########################################################################
			return *this;
		}
		if (_dim != scalar.getDim())
		{
			//########################################################################
			_log->ERROR("Attempted to add scalar fields " + _name  + " and "
			            + scalar.getName() + " with dimensions "
									+ std::to_string(_dim) + " and "
									+ std::to_string(scalar.getDim()));
			//########################################################################
			return *this;
		}
		if (_ugrid != scalar.getUGrid())
		{
			//########################################################################
			_log->WARN("UGrids for scalar fields " + _name  + " and "
			            + scalar.getName() + " do not match");
			//########################################################################
		}
		std::vector<T> _copy = _field;
		std::transform(_copy.begin(), _copy.end(), scalar.getField().begin(),
	                 _copy.begin(), std::plus<T>());
		std::string name;
		if (_name != " " && scalar.getName() != " ")
		{
			name += "(" + _name + " + " + scalar.getName() + ")";
		}
		return ScalarField<T>(name,_ugrid,_copy,_log);
	}
	//----------------------------------------------------------------------------
	template<typename T>
	ScalarField<T> ScalarField<T>::operator-(const ScalarField<T>& scalar) const
	{
		if (_N != scalar.getN())
		{
			//########################################################################
			_log->ERROR("Attempted to subtract scalar fields " + _name  + " and "
			            + scalar.getName() + " with sizes " + std::to_string(_N)
									+ " and " + std::to_string(scalar.getN()));
			//########################################################################
			return *this;
		}
		if (_dim != scalar.getDim())
		{
			//########################################################################
			_log->ERROR("Attempted to subtract scalar fields " + _name  + " and "
			            + scalar.getName() + " with dimensions "
									+ std::to_string(_dim) + " and "
									+ std::to_string(scalar.getDim()));
			//########################################################################
			return *this;
		}
		if (_ugrid != scalar.getUGrid())
		{
			//########################################################################
			_log->WARN("UGrids for scalar fields " + _name  + " and "
			            + scalar.getName() + " do not match");
			//########################################################################
		}
		std::vector<T> _copy = _field;
		std::transform(_copy.begin(), _copy.end(), scalar.getField().begin(),
	                 _copy.begin(), std::minus<T>());
		std::string name;
		if (_name != " " && scalar.getName() != " ")
		{
			name += "(" + _name + " - " + scalar.getName() + ")";
		}
		return ScalarField<T>(name,_ugrid,_copy,_log);
	}
	//----------------------------------------------------------------------------
	template<typename T>
	ScalarField<T> ScalarField<T>::operator*(const ScalarField<T>& scalar) const
	{
		if (_N != scalar.getN())
		{
			//########################################################################
			_log->ERROR("Attempted to multiply scalar fields " + _name  + " and "
			            + scalar.getName() + " with sizes " + std::to_string(_N)
									+ " and " + std::to_string(scalar.getN()));
			//########################################################################
			return *this;
		}
		if (_dim != scalar.getDim())
		{
			//########################################################################
			_log->ERROR("Attempted to multiply scalar fields " + _name  + " and "
			            + scalar.getName() + " with dimensions "
									+ std::to_string(_dim) + " and "
									+ std::to_string(scalar.getDim()));
			//########################################################################
			return *this;
		}
		if (_ugrid != scalar.getUGrid())
		{
			//########################################################################
			_log->WARN("UGrids for scalar fields " + _name  + " and "
			            + scalar.getName() + " do not match");
			//########################################################################
		}
		std::vector<T> _copy = _field;
		std::transform(_copy.begin(), _copy.end(), scalar.getField().begin(),
	                 _copy.begin(), std::multiplies<T>());
		std::string name;
		if (_name != " " && scalar.getName() != " ")
		{
			name += "(" + _name + " * " + scalar.getName() + ")";
		}
		return ScalarField<T>(name,_ugrid,_copy,_log);
	}
	//----------------------------------------------------------------------------
	template<typename T>
	ScalarField<T> ScalarField<T>::operator/(const ScalarField<T>& scalar) const
	{
		if (_N != scalar.getN())
		{
			//########################################################################
			_log->ERROR("Attempted to divide scalar fields " + _name  + " and "
			            + scalar.getName() + " with sizes " + std::to_string(_N)
									+ " and " + std::to_string(scalar.getN()));
			//########################################################################
			return *this;
		}
		if (_dim != scalar.getDim())
		{
			//########################################################################
			_log->ERROR("Attempted to divide scalar fields " + _name  + " and "
			            + scalar.getName() + " with dimensions "
									+ std::to_string(_dim) + " and "
									+ std::to_string(scalar.getDim()));
			//########################################################################
			return *this;
		}
		if (_ugrid != scalar.getUGrid())
		{
			//########################################################################
			_log->WARN("UGrids for scalar fields " + _name  + " and "
			            + scalar.getName() + " do not match");
			//########################################################################
		}
		std::vector<T> _copy = _field;
		std::transform(_copy.begin(), _copy.end(), scalar.getField().begin(),
	                 _copy.begin(), std::divides<T>());
		std::string name;
		if (_name != " " && scalar.getName() != " ")
		{
			name += "(" + _name + " / " + scalar.getName() + ")";
		}
		return ScalarField<T>(name,_ugrid,_copy,_log);
	}
	//----------------------------------------------------------------------------
	template<typename T>
	ScalarField<T>& ScalarField<T>::operator+=(const ScalarField<T>& scalar)
	{
		if (_N != scalar.getN())
		{
			//########################################################################
			_log->ERROR("Attempted to add scalar fields " + _name  + " and "
			            + scalar.getName() + " with sizes " + std::to_string(_N)
									+ " and " + std::to_string(scalar.getN()));
			//########################################################################
			return *this;
		}
		if (_dim != scalar.getDim())
		{
			//########################################################################
			_log->ERROR("Attempted to add scalar fields " + _name  + " and "
			            + scalar.getName() + " with dimensions "
									+ std::to_string(_dim) + " and "
									+ std::to_string(scalar.getDim()));
			//########################################################################
			return *this;
		}
		if (_ugrid != scalar.getUGrid())
		{
			//########################################################################
			_log->WARN("UGrids for scalar fields " + _name  + " and "
			            + scalar.getName() + " do not match");
			//########################################################################
		}
		std::vector<T> _copy = _field;
		std::transform(_copy.begin(), _copy.end(), scalar.getField().begin(),
	                 _copy.begin(), std::plus<T>());
		std::string name;
		if (_name != " " && scalar.getName() != " ")
		{
			name += "(" + _name + " + " + scalar.getName() + ")";
		}
		_name = name;
		_field = _copy;
		return *this;
	}
	//----------------------------------------------------------------------------
	template<typename T>
	ScalarField<T>& ScalarField<T>::operator-=(const ScalarField<T>& scalar)
	{
		if (_N != scalar.getN())
		{
			//########################################################################
			_log->ERROR("Attempted to subtract scalar fields " + _name  + " and "
									+ scalar.getName() + " with sizes " + std::to_string(_N)
									+ " and " + std::to_string(scalar.getN()));
			//########################################################################
			return *this;
		}
		if (_dim != scalar.getDim())
		{
			//########################################################################
			_log->ERROR("Attempted to subtract scalar fields " + _name  + " and "
									+ scalar.getName() + " with dimensions "
									+ std::to_string(_dim) + " and "
									+ std::to_string(scalar.getDim()));
			//########################################################################
			return *this;
		}
		if (_ugrid != scalar.getUGrid())
		{
			//########################################################################
			_log->WARN("UGrids for scalar fields " + _name  + " and "
			            + scalar.getName() + " do not match");
			//########################################################################
		}
		std::vector<T> _copy = _field;
		std::transform(_copy.begin(), _copy.end(), scalar.getField().begin(),
									 _copy.begin(), std::minus<T>());
		std::string name;
		if (_name != " " && scalar.getName() != " ")
		{
			name += "(" + _name + " - " + scalar.getName() + ")";
		}
		_name = name;
		_field = _copy;
		return *this;
	}
	//----------------------------------------------------------------------------
	template<typename T>
	ScalarField<T>& ScalarField<T>::operator*=(const ScalarField<T>& scalar)
	{
		if (_N != scalar.getN())
		{
			//########################################################################
			_log->ERROR("Attempted to multiply scalar fields " + _name  + " and "
									+ scalar.getName() + " with sizes " + std::to_string(_N)
									+ " and " + std::to_string(scalar.getN()));
			//########################################################################
			return *this;
		}
		if (_dim != scalar.getDim())
		{
			//########################################################################
			_log->ERROR("Attempted to multiply scalar fields " + _name  + " and "
									+ scalar.getName() + " with dimensions "
									+ std::to_string(_dim) + " and "
									+ std::to_string(scalar.getDim()));
			//########################################################################
			return *this;
		}
		if (_ugrid != scalar.getUGrid())
		{
			//########################################################################
			_log->WARN("UGrids for scalar fields " + _name  + " and "
			            + scalar.getName() + " do not match");
			//########################################################################
		}
		std::vector<T> _copy = _field;
		std::transform(_copy.begin(), _copy.end(), scalar.getField().begin(),
									 _copy.begin(), std::multiplies<T>());
		std::string name;
		if (_name != " " && scalar.getName() != " ")
		{
			name += "(" + _name + " * " + scalar.getName() + ")";
		}
		_name = name;
		_field = _copy;
		return *this;
	}
	//----------------------------------------------------------------------------
	template<typename T>
	ScalarField<T>& ScalarField<T>::operator/=(const ScalarField<T>& scalar)
	{
		if (_N != scalar.getN())
		{
			//########################################################################
			_log->ERROR("Attempted to divide scalar fields " + _name  + " and "
									+ scalar.getName() + " with sizes " + std::to_string(_N)
									+ " and " + std::to_string(scalar.getN()));
			//########################################################################
			return *this;
		}
		if (_dim != scalar.getDim())
		{
			//########################################################################
			_log->ERROR("Attempted to divide scalar fields " + _name  + " and "
									+ scalar.getName() + " with dimensions "
									+ std::to_string(_dim) + " and "
									+ std::to_string(scalar.getDim()));
			//########################################################################
			return *this;
		}
		if (_ugrid != scalar.getUGrid())
		{
			//########################################################################
			_log->WARN("UGrids for scalar fields " + _name  + " and "
			            + scalar.getName() + " do not match");
			//########################################################################
		}
		std::vector<T> _copy = _field;
		std::transform(_copy.begin(), _copy.end(), scalar.getField().begin(),
									 _copy.begin(), std::divides<T>());
		std::string name;
		if (_name != " " && scalar.getName() != " ")
		{
			name += "(" + _name + " / " + scalar.getName() + ")";
		}
		_name = name;
		_field = _copy;
		return *this;
	}
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
    return _approx->scalarGradient(_ugrid,(*this));
  }
  //----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	nth-derivative for every point in all directions
	//----------------------------------------------------------------------------
	template<typename T>
	std::vector<std::vector<T>>
	ScalarField<T>::derivative(uint32_t n)
	{
		return _approx->scalarDerivative(_ugrid, (*this), n);
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	nth-derivative in the ith-direction for every point
	//----------------------------------------------------------------------------
	template<typename T>
	std::vector<T> ScalarField<T>::derivative(uint32_t dir, uint32_t n)
	{
		return _approx->scalarDerivative(_ugrid, (*this), dir, n);
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	nth-derivative in the ith-direction for every point
	//----------------------------------------------------------------------------
	template<typename T>
	std::vector<T> ScalarField<T>::derivative(std::vector<uint32_t> deriv)
	{
		return _approx->scalarDerivative(_ugrid, (*this), deriv);
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	nth-derivative for a single point
	//----------------------------------------------------------------------------
	template<typename T>
	std::vector<T> ScalarField<T>::derivativePoint(uint64_t index, uint32_t n)
	{
		return _approx->scalarDerivative(_ugrid, (*this), index, n);
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	nth-derivative in the ith-direction for a single point
	//----------------------------------------------------------------------------
	template<typename T>
	T ScalarField<T>::derivativePoint(uint64_t index, uint32_t dir, uint32_t n)
	{
		return _approx->scalarDerivative(_ugrid, (*this), index, dir, n);
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	nth-derivative in the ith-direction for a single point
	//----------------------------------------------------------------------------
	template<typename T>
	T ScalarField<T>::derivativePoint(uint64_t index, std::vector<uint32_t> deriv)
	{
		return _approx->scalarDerivative(_ugrid, (*this), index, deriv);
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Laplacian for every point
	//----------------------------------------------------------------------------
	template<typename T>
	std::vector<T> ScalarField<T>::laplacian()
	{
		std::vector<T> result(_N);
		std::vector<std::vector<T>> second_derivative = derivative(2);
		for (uint64_t i = 0; i < _N; i++)
		{
			result[i] = std::accumulate(second_derivative[i].begin(),
		                              second_derivative[i].end(),0);
		}
		return result;
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Laplacian for a single point
	//----------------------------------------------------------------------------
	template<typename T>
	T ScalarField<T>::laplacian(uint64_t index)
	{
		std::vector<T> second_derivative = derivative(index, 2);
		return std::accumulate(second_derivative.begin(),second_derivative.end(),0);
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Various functions
	//----------------------------------------------------------------------------
	template<typename T>
	std::string ScalarField<T>::summary()
	{
		std::string sum = "---------------------------------------------------";
		sum += "\n<ET::ScalarField<"+ type_name<decltype(_field[0])>();
		sum += "> object at " + getMem(this) + ">";
		sum += "\n---------------------------------------------------";
    sum += "\n   name: '" + _name + "'";
		sum += "\n    dim: " + std::to_string(_dim);
		sum += "\n      N: " + std::to_string(_N);
		sum += "\n---------------------------------------------------";
		std::string grid = "UGrid '" + _ugrid->getName() + "' at: ";
		int grid_colon_loc = grid.length()-8;
		sum += "\n" + grid + getMem(*getLogger()) + ",";
		std::string grid_ref;
		grid_ref.resize(grid_colon_loc,' ');
		sum += "\n" + grid_ref;
		sum += "ref at: " + getMem(_ugrid);
		sum += "\nApproximator at: " + getMem(*getApproximator()) + ",";
		sum += "\n         ref at: " + getMem(_approx);
		sum += "\nLogger at: " + getMem(*getLogger()) + ",";
		sum += "\n   ref at: " + getMem(_log);
		sum += "\n++++++++++++++++++++++++++++++++++++++++++++++++++++";
		return sum;
	}
	//----------------------------------------------------------------------------
}
