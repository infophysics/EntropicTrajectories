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
  template<typename T>
  ScalarField<T>::ScalarField()
  : Field();
  {
    m_dim = 1;
  }
  template<typename T>
  ScalarField<T>::~ScalarField()
  : ~Field()
  {
    m_dim = 1;
  }
  template<typename T>
  ScalarField<T>::ScalarField(std::shared_ptr<Log> t_log)
  : Field(t_log)
  {
    m_dim = 1;
  }
  template<typename T>
  ScalarField<T>::ScalarField(std::shared_ptr<UGrid<T>> t_ugrid)
  : Field(t_ugrid)
  {
    m_dim = 1;
  }
  template<typename T>
  ScalarField<T>::ScalarField(std::shared_ptr<Interpolator<T>> t_interpolator)
  : Field(t_interpolator)
  {
    m_dim = 1;
  }
  template<typename T>
  ScalarField<T>::ScalarField(std::shared_ptr<UGrid<T>> t_ugrid,
                              std::shared_ptr<Log> t_log)
  : Field(t_ugrid, t_log)
  {
    m_dim = 1;
  }
  template<typename T>
  ScalarField<T>::ScalarField(std::vector<T> t_field)
  : Field()
  {
    m_field = t_field;
    m_dim = 1;
    m_N = m_field.size();
  }
  template<typename T>
  ScalarField<T>::ScalarField(std::vector<T> t_field,
                              std::shared_ptr<Log> t_log)
  : Field(t_log)
  {
    m_field = t_field;
    m_dim = 1;
    m_N = m_field.size();
  }
  template<typename T>
  ScalarField<T>::ScalarField(std::vector<T> t_field,
                              std::shared_ptr<UGrid<T>> t_ugrid)
  : Field(t_ugrid)
  {
    m_field = t_field;
    m_dim = 1;
    m_N = m_field.size();
  }
  template<typename T>
  ScalarField<T>::ScalarField(std::vector<T> t_field,
                              std::shared_ptr<Interpolator<T>> t_interpolator)
  : Field(t_interpolator)
  {
    m_field = t_field;
    m_dim = 1;
    m_N = m_field.size();
  }
  template<typename T>
  ScalarField<T>::ScalarField(std::vector<T> t_field,
                              std::shared_ptr<UGrid<T>> t_ugrid,
                              std::shared_ptr<Log> t_log)
  : Field(t_ugrid, t_log)
  {
    m_field = t_field;
    m_dim = 1;
    m_N = m_field.size();
  }

  template<typename T>
  std::vector<T> ScalarField<T>::getField() const
  {
    return m_field;
  }
  template<typename T>
  std::vector<T>* ScalarField<T>::accessField()
  {
    return &m_field;
  }
  template<typename T>
  T* ScalarField<T>::data()
  {
    return m_field.data();
  }
  template<typename T>
  void ScalarField<T>::setField(std::vector<T> t_field)
  {
    m_field = t_field;
    m_N = m_field.size();
  }

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

  template<typename T>
  Vector<T> ScalarField<T>::constructLocalFieldValues(size_t t_index)
  {
    //  Get the nearest neighbors for the index
    u_grid->queryNeighbors(t_index);
    std::vector<size_t> neighbors = ugrid->getNeighbors(t_index);
    //  Create the empty vector
    Vector<T> f(neighbors.size());
    for (auto i = 0; i < neighbors.size(); i++) {
      f(i) = m_field[neighbors[i]];
    }
    return f;
  }

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
		return _approx->scalarDerivativePoint(_ugrid, (*this), index, n);
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	nth-derivative in the ith-direction for a single point
	//----------------------------------------------------------------------------
	template<typename T>
	T ScalarField<T>::derivativePoint(uint64_t index, uint32_t dir, uint32_t n)
	{
		return _approx->scalarDerivativePoint(_ugrid, (*this), index, dir, n);
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	nth-derivative in the ith-direction for a single point
	//----------------------------------------------------------------------------
	template<typename T>
	T ScalarField<T>::derivativePoint(uint64_t index, std::vector<uint32_t> deriv)
	{
		return _approx->scalarDerivativePoint(_ugrid, (*this), index, deriv);
	}
	//----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
	//	nth-derivative for a single point
	//----------------------------------------------------------------------------
	template<typename T>
	std::vector<T>
  ScalarField<T>::derivativePoint(std::vector<T> point, uint32_t n)
	{
		return _approx->scalarDerivativePoint(_ugrid, (*this), point, n);
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	nth-derivative in the ith-direction for a single point
	//----------------------------------------------------------------------------
	template<typename T>
	T ScalarField<T>::derivativePoint(std::vector<T> point, uint32_t dir,
                                    uint32_t n)
	{
		return _approx->scalarDerivativePoint(_ugrid, (*this), point, dir, n);
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	nth-derivative in the ith-direction for a single point
	//----------------------------------------------------------------------------
	template<typename T>
	T ScalarField<T>::derivativePoint(std::vector<T> point,
                                    std::vector<uint32_t> deriv)
	{
		return _approx->scalarDerivativePoint(_ugrid, (*this), point, deriv);
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
	const std::string ScalarField<T>::summary()
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
		sum += "\nInterpolator at: " + getMem(*getInterpolator()) + ",";
		sum += "\n         ref at: " + getMem(_approx);
		sum += "\nLogger at: " + getMem(*getLogger()) + ",";
		sum += "\n   ref at: " + getMem(_log);
		sum += "\n++++++++++++++++++++++++++++++++++++++++++++++++++++";
		return sum;
	}
	//----------------------------------------------------------------------------

}
