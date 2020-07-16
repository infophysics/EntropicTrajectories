//------------------------------------------------------------------------------
//  t_scalarfield.cpp
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
  : Field<T>()
  {
    this->m_dim = 1;
  }
  template<typename T>
  ScalarField<T>::~ScalarField()
  {
  }
  template<typename T>
  ScalarField<T>::ScalarField(std::shared_ptr<Log> t_log)
  : Field<T>(t_log)
  {
    this->m_dim = 1;
  }
  template<typename T>
  ScalarField<T>::ScalarField(std::shared_ptr<UGrid<T>> t_ugrid)
  : Field<T>(t_ugrid)
  {
    this->m_dim = 1;
  }
  template<typename T>
  ScalarField<T>::ScalarField(std::shared_ptr<Interpolator<T>> t_interpolator)
  : Field<T>(t_interpolator)
  {
    this->m_dim = 1;
  }
  template<typename T>
  ScalarField<T>::ScalarField(std::shared_ptr<UGrid<T>> t_ugrid,
                              std::shared_ptr<Log> t_log)
  : Field<T>(t_ugrid, t_log)
  {
    this->m_dim = 1;
  }
  template<typename T>
  ScalarField<T>::ScalarField(std::vector<T> t_field)
  : Field<T>()
  {
    m_field = t_field;
    this->m_dim = 1;
    this->m_N = m_field.size();
  }
  template<typename T>
  ScalarField<T>::ScalarField(std::vector<T> t_field,
                              std::shared_ptr<Log> t_log)
  : Field<T>(t_log)
  {
    m_field = t_field;
    this->m_dim = 1;
    this->m_N = m_field.size();
  }
  template<typename T>
  ScalarField<T>::ScalarField(std::vector<T> t_field,
                              std::shared_ptr<UGrid<T>> t_ugrid)
  : Field<T>(t_ugrid)
  {
    m_field = t_field;
    this->m_dim = 1;
    this->m_N = m_field.size();
  }
  template<typename T>
  ScalarField<T>::ScalarField(std::vector<T> t_field,
                              std::shared_ptr<Interpolator<T>> t_interpolator)
  : Field<T>(t_interpolator)
  {
    m_field = t_field;
    this->m_dim = 1;
    this->m_N = m_field.size();
  }
  template<typename T>
  ScalarField<T>::ScalarField(std::vector<T> t_field,
                              std::shared_ptr<UGrid<T>> t_ugrid,
                              std::shared_ptr<Log> t_log)
  : Field<T>(t_ugrid, t_log)
  {
    m_field = t_field;
    this->m_dim = 1;
    this->m_N = m_field.size();
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
    this->m_N = m_field.size();
  }

  //----------------------------------------------------------------------------
  //  Operator overloads
  //----------------------------------------------------------------------------
	template<typename T>
	ScalarField<T> ScalarField<T>::operator+(const ScalarField<T>& t_scalar) const
	{
		if (this->m_N != t_scalar.getN()) {
			this->m_log->ERROR("Attempted to add t_scalar fields " + this->m_name  + " and "
			            + t_scalar.getName() + " with sizes " + std::to_string(this->m_N)
									+ " and " + std::to_string(t_scalar.getN()));
			return *this;
		}
		if (this->m_dim != t_scalar.getDim()) {
			this->m_log->ERROR("Attempted to add t_scalar fields " + this->m_name  + " and "
			            + t_scalar.getName() + " with dimensions "
									+ std::to_string(this->m_dim) + " and "
									+ std::to_string(t_scalar.getDim()));
			return *this;
		}
		if (this->m_ugrid != t_scalar.getUGrid()) {
			this->m_log->WARN("UGrids for t_scalar fields " + this->m_name  + " and "
			            + t_scalar.getName() + " do not match");
		}
		std::vector<T> copy = m_field;
		std::transform(copy.begin(), copy.end(), t_scalar.getField().begin(),
	                 copy.begin(), std::plus<T>());
		std::string name;
		if (this->m_name != " " && t_scalar.getName() != " ") {
			name += "(" + this->m_name + " + " + t_scalar.getName() + ")";
		}
		return ScalarField<T>(copy,this->m_ugrid,this->m_log);
	}
	//----------------------------------------------------------------------------
	template<typename T>
	ScalarField<T> ScalarField<T>::operator-(const ScalarField<T>& t_scalar) const
	{
		if (this->m_N != t_scalar.getN()) {
			this->m_log->ERROR("Attempted to subtract t_scalar fields " + this->m_name  + " and "
			            + t_scalar.getName() + " with sizes " + std::to_string(this->m_N)
									+ " and " + std::to_string(t_scalar.getN()));
			return *this;
		}
		if (this->m_dim != t_scalar.getDim()) {
			this->m_log->ERROR("Attempted to subtract t_scalar fields " + this->m_name  + " and "
			            + t_scalar.getName() + " with dimensions "
									+ std::to_string(this->m_dim) + " and "
									+ std::to_string(t_scalar.getDim()));
			return *this;
		}
		if (this->m_ugrid != t_scalar.getUGrid()) {
			this->m_log->WARN("UGrids for t_scalar fields " + this->m_name  + " and "
			            + t_scalar.getName() + " do not match");
		}
		std::vector<T> copy = m_field;
		std::transform(copy.begin(), copy.end(), t_scalar.getField().begin(),
	                 copy.begin(), std::minus<T>());
		std::string name;
		if (this->m_name != " " && t_scalar.getName() != " ") {
			name += "(" + this->m_name + " - " + t_scalar.getName() + ")";
		}
		return ScalarField<T>(copy,this->m_ugrid,this->m_log);
	}
	//----------------------------------------------------------------------------
	template<typename T>
	ScalarField<T> ScalarField<T>::operator*(const ScalarField<T>& t_scalar) const
	{
		if (this->m_N != t_scalar.getN()) {
			this->m_log->ERROR("Attempted to multiply t_scalar fields " + this->m_name  + " and "
			            + t_scalar.getName() + " with sizes " + std::to_string(this->m_N)
									+ " and " + std::to_string(t_scalar.getN()));
			return *this;
		}
		if (this->m_dim != t_scalar.getDim()) {
			this->m_log->ERROR("Attempted to multiply t_scalar fields " + this->m_name  + " and "
			            + t_scalar.getName() + " with dimensions "
									+ std::to_string(this->m_dim) + " and "
									+ std::to_string(t_scalar.getDim()));
			return *this;
		}
		if (this->m_ugrid != t_scalar.getUGrid()) {
			this->m_log->WARN("UGrids for t_scalar fields " + this->m_name  + " and "
			            + t_scalar.getName() + " do not match");
		}
		std::vector<T> copy = m_field;
		std::transform(copy.begin(), copy.end(), t_scalar.getField().begin(),
	                 copy.begin(), std::multiplies<T>());
		std::string name;
		if (this->m_name != " " && t_scalar.getName() != " ") {
			name += "(" + this->m_name + " * " + t_scalar.getName() + ")";
		}
		return ScalarField<T>(copy,this->m_ugrid,this->m_log);
	}
	//----------------------------------------------------------------------------
	template<typename T>
	ScalarField<T> ScalarField<T>::operator/(const ScalarField<T>& t_scalar) const
	{
		if (this->m_N != t_scalar.getN()) {
			this->m_log->ERROR("Attempted to divide t_scalar fields " + this->m_name  + " and "
			            + t_scalar.getName() + " with sizes " + std::to_string(this->m_N)
									+ " and " + std::to_string(t_scalar.getN()));
			return *this;
		}
		if (this->m_dim != t_scalar.getDim()) {
			this->m_log->ERROR("Attempted to divide t_scalar fields " + this->m_name  + " and "
			            + t_scalar.getName() + " with dimensions "
									+ std::to_string(this->m_dim) + " and "
									+ std::to_string(t_scalar.getDim()));
			return *this;
		}
		if (this->m_ugrid != t_scalar.getUGrid()) {
			this->m_log->WARN("UGrids for t_scalar fields " + this->m_name  + " and "
			            + t_scalar.getName() + " do not match");
		}
		std::vector<T> copy = m_field;
		std::transform(copy.begin(), copy.end(), t_scalar.getField().begin(),
	                 copy.begin(), std::divides<T>());
		std::string name;
		if (this->m_name != " " && t_scalar.getName() != " ") {
			name += "(" + this->m_name + " / " + t_scalar.getName() + ")";
		}
		return ScalarField<T>(copy,this->m_ugrid,this->m_log);
	}
	//----------------------------------------------------------------------------
	template<typename T>
	ScalarField<T>& ScalarField<T>::operator+=(const ScalarField<T>& t_scalar)
	{
		if (this->m_N != t_scalar.getN()) {
			this->m_log->ERROR("Attempted to add t_scalar fields " + this->m_name  + " and "
			            + t_scalar.getName() + " with sizes " + std::to_string(this->m_N)
									+ " and " + std::to_string(t_scalar.getN()));
			return *this;
		}
		if (this->m_dim != t_scalar.getDim()) {
			this->m_log->ERROR("Attempted to add t_scalar fields " + this->m_name  + " and "
			            + t_scalar.getName() + " with dimensions "
									+ std::to_string(this->m_dim) + " and "
									+ std::to_string(t_scalar.getDim()));
			return *this;
		}
		if (this->m_ugrid != t_scalar.getUGrid()) {
			this->m_log->WARN("UGrids for t_scalar fields " + this->m_name  + " and "
			            + t_scalar.getName() + " do not match");
		}
		std::vector<T> copy = m_field;
		std::transform(copy.begin(), copy.end(), t_scalar.getField().begin(),
	                 copy.begin(), std::plus<T>());
		std::string name;
		if (this->m_name != " " && t_scalar.getName() != " ") {
			name += "(" + this->m_name + " + " + t_scalar.getName() + ")";
		}
		this->m_name = name;
		m_field = copy;
		return *this;
	}
	//----------------------------------------------------------------------------
	template<typename T>
	ScalarField<T>& ScalarField<T>::operator-=(const ScalarField<T>& t_scalar)
	{
		if (this->m_N != t_scalar.getN()) {
			this->m_log->ERROR("Attempted to subtract t_scalar fields " + this->m_name  + " and "
									+ t_scalar.getName() + " with sizes " + std::to_string(this->m_N)
									+ " and " + std::to_string(t_scalar.getN()));
			return *this;
		}
		if (this->m_dim != t_scalar.getDim()) {
			this->m_log->ERROR("Attempted to subtract t_scalar fields " + this->m_name  + " and "
									+ t_scalar.getName() + " with dimensions "
									+ std::to_string(this->m_dim) + " and "
									+ std::to_string(t_scalar.getDim()));
			return *this;
		}
		if (this->m_ugrid != t_scalar.getUGrid()) {
			this->m_log->WARN("UGrids for t_scalar fields " + this->m_name  + " and "
			            + t_scalar.getName() + " do not match");
		}
		std::vector<T> copy = m_field;
		std::transform(copy.begin(), copy.end(), t_scalar.getField().begin(),
									 copy.begin(), std::minus<T>());
		std::string name;
		if (this->m_name != " " && t_scalar.getName() != " ") {
			name += "(" + this->m_name + " - " + t_scalar.getName() + ")";
		}
		this->m_name = name;
		m_field = copy;
		return *this;
	}
	//----------------------------------------------------------------------------
	template<typename T>
	ScalarField<T>& ScalarField<T>::operator*=(const ScalarField<T>& t_scalar)
	{
		if (this->m_N != t_scalar.getN()) {
			this->m_log->ERROR("Attempted to multiply t_scalar fields " + this->m_name  + " and "
									+ t_scalar.getName() + " with sizes " + std::to_string(this->m_N)
									+ " and " + std::to_string(t_scalar.getN()));
			return *this;
		}
		if (this->m_dim != t_scalar.getDim()) {
			this->m_log->ERROR("Attempted to multiply t_scalar fields " + this->m_name  + " and "
									+ t_scalar.getName() + " with dimensions "
									+ std::to_string(this->m_dim) + " and "
									+ std::to_string(t_scalar.getDim()));
			return *this;
		}
		if (this->m_ugrid != t_scalar.getUGrid()) {
			this->m_log->WARN("UGrids for t_scalar fields " + this->m_name  + " and "
			            + t_scalar.getName() + " do not match");
		}
		std::vector<T> copy = m_field;
		std::transform(copy.begin(), copy.end(), t_scalar.getField().begin(),
									 copy.begin(), std::multiplies<T>());
		std::string name;
		if (this->m_name != " " && t_scalar.getName() != " ") {
			name += "(" + this->m_name + " * " + t_scalar.getName() + ")";
		}
		this->m_name = name;
		m_field = copy;
		return *this;
	}
	//----------------------------------------------------------------------------
	template<typename T>
	ScalarField<T>& ScalarField<T>::operator/=(const ScalarField<T>& t_scalar)
	{
		if (this->m_N != t_scalar.getN()) {
			this->m_log->ERROR("Attempted to divide t_scalar fields " + this->m_name  + " and "
									+ t_scalar.getName() + " with sizes " + std::to_string(this->m_N)
									+ " and " + std::to_string(t_scalar.getN()));
			return *this;
		}
		if (this->m_dim != t_scalar.getDim()) {
			this->m_log->ERROR("Attempted to divide t_scalar fields " + this->m_name  + " and "
									+ t_scalar.getName() + " with dimensions "
									+ std::to_string(this->m_dim) + " and "
									+ std::to_string(t_scalar.getDim()));
			return *this;
		}
		if (this->m_ugrid != t_scalar.getUGrid()) {
			this->m_log->WARN("UGrids for t_scalar fields " + this->m_name  + " and "
			            + t_scalar.getName() + " do not match");
		}
		std::vector<T> copy = m_field;
		std::transform(copy.begin(), copy.end(), t_scalar.getField().begin(),
									 copy.begin(), std::divides<T>());
		std::string name;
		if (this->m_name != " " && t_scalar.getName() != " ") {
			name += "(" + this->m_name + " / " + t_scalar.getName() + ")";
		}
		this->m_name = name;
		m_field = copy;
		return *this;
	}
	//----------------------------------------------------------------------------
	template<typename T>
  T& ScalarField<T>::operator()(const size_t& i)
  {
		if (i >= this->m_N) {
			this->m_log->ERROR("ScalarField " + this->m_name
									+ ": Attempted to access m_field array of size "
									+ std::to_string(this->m_N) + " with index "
									+ std::to_string(i));
			if(m_field.size() > 0) {
				this->m_log->INFO("ScalarField "+ this->m_name +": Returning the element at index 0");
				return m_field[0];
			}
			else {
				this->m_log->INFO("ScalarField " + this->m_name + ": Terminating program");
				exit(0);
			}
		}
    return m_field[i];
  }
  template<typename T>
  const T& ScalarField<T>::operator()(const size_t& i) const
  {
		if (i >= this->m_N) {
			this->m_log->ERROR("ScalarField " + this->m_name
									+ ": Attempted to access m_field array of size "
									+ std::to_string(this->m_N) + " with index "
									+ std::to_string(i));
			if(m_field.size() > 0) {
				this->m_log->INFO("ScalarField "+ this->m_name +": Returning the element at index 0");
				return m_field[0];
			}
			else {
				this->m_log->INFO("ScalarField " + this->m_name + ": Terminating program");
				exit(0);
			}
		}
    return m_field[i];
  }
  //----------------------------------------------------------------------------

  template<typename T>
  Vector<T> ScalarField<T>::constructLocalFieldValues(size_t t_index)
  {
    //  Get the nearest neighbors for the index
    this->m_ugrid->queryNeighbors(t_index);
    std::vector<size_t> neighbors = this->m_ugrid->getNeighbors(t_index);
    //  Create the empty vector
    Vector<T> f(neighbors.size());
    for (auto i = 0; i < neighbors.size(); i++) {
      f(i) = m_field[neighbors[i]];
    }
    return f;
  }
  template<typename T>
  Vector<T>
  ScalarField<T>::constructLocalFieldValues(const std::vector<T>& t_point,
                                            size_t t_k)
  {
    //  Get the nearest neighbors for the index
    std::vector<size_t> neighbors = this->m_ugrid->queryNeighbors(t_point, t_k);
    //  Create the empty vector
    Vector<T> f(neighbors.size());
    for (auto i = 0; i < neighbors.size(); i++) {
      f(i) = m_field[neighbors[i]];
    }
    return f;
  }

  // template<typename T>
  // std::vector<std::vector<T>> ScalarField<T>::gradient()
  // {
  //   return _approx->t_scalarGradient(this->m_ugrid,(*this));
  // }
  // //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //	nth-derivative for every point in all directions
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<std::vector<T>>
	// ScalarField<T>::derivative(size_t n)
	// {
	// 	return _approx->t_scalarDerivative(this->m_ugrid, (*this), n);
	// }
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //	nth-derivative in the ith-direction for every point
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<T> ScalarField<T>::derivative(size_t dir, size_t n)
	// {
	// 	return _approx->t_scalarDerivative(this->m_ugrid, (*this), dir, n);
	// }
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //	nth-derivative in the ith-direction for every point
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<T> ScalarField<T>::derivative(std::vector<size_t> deriv)
	// {
	// 	return _approx->t_scalarDerivative(this->m_ugrid, (*this), deriv);
	// }
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //	nth-derivative for a single point
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<T> ScalarField<T>::derivativePoint(size_t index, size_t n)
	// {
	// 	return _approx->t_scalarDerivativePoint(this->m_ugrid, (*this), index, n);
	// }
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //	nth-derivative in the ith-direction for a single point
	// //----------------------------------------------------------------------------
	// template<typename T>
	// T ScalarField<T>::derivativePoint(size_t index, size_t dir, size_t n)
	// {
	// 	return _approx->t_scalarDerivativePoint(this->m_ugrid, (*this), index, dir, n);
	// }
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //	nth-derivative in the ith-direction for a single point
	// //----------------------------------------------------------------------------
	// template<typename T>
	// T ScalarField<T>::derivativePoint(size_t index, std::vector<size_t> deriv)
	// {
	// 	return _approx->t_scalarDerivativePoint(this->m_ugrid, (*this), index, deriv);
	// }
	// //----------------------------------------------------------------------------
  //
  // //----------------------------------------------------------------------------
	// //	nth-derivative for a single point
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<T>
  // ScalarField<T>::derivativePoint(std::vector<T> point, size_t n)
	// {
	// 	return _approx->t_scalarDerivativePoint(this->m_ugrid, (*this), point, n);
	// }
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //	nth-derivative in the ith-direction for a single point
	// //----------------------------------------------------------------------------
	// template<typename T>
	// T ScalarField<T>::derivativePoint(std::vector<T> point, size_t dir,
  //                                   size_t n)
	// {
	// 	return _approx->t_scalarDerivativePoint(this->m_ugrid, (*this), point, dir, n);
	// }
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //	nth-derivative in the ith-direction for a single point
	// //----------------------------------------------------------------------------
	// template<typename T>
	// T ScalarField<T>::derivativePoint(std::vector<T> point,
  //                                   std::vector<size_t> deriv)
	// {
	// 	return _approx->t_scalarDerivativePoint(this->m_ugrid, (*this), point, deriv);
	// }
	// //----------------------------------------------------------------------------
  //
  //
	// //----------------------------------------------------------------------------
	// //	Laplacian for every point
	// //----------------------------------------------------------------------------
	// template<typename T>
	// std::vector<T> ScalarField<T>::laplacian()
	// {
	// 	std::vector<T> result(this->m_N);
	// 	std::vector<std::vector<T>> second_derivative = derivative(2);
	// 	for (size_t i = 0; i < _N; i++)
	// 	{
	// 		result[i] = std::accumulate(second_derivative[i].begin(),
	// 	                              second_derivative[i].end(),0);
	// 	}
	// 	return result;
	// }
	// //----------------------------------------------------------------------------
  //
	// //----------------------------------------------------------------------------
	// //	Laplacian for a single point
	// //----------------------------------------------------------------------------
	// template<typename T>
	// T ScalarField<T>::laplacian(size_t index)
	// {
	// 	std::vector<T> second_derivative = derivative(index, 2);
	// 	return std::accumulate(second_derivative.begin(),second_derivative.end(),0);
	// }
	// //----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Various functions
	//----------------------------------------------------------------------------
	template<typename T>
	const std::string ScalarField<T>::summary()
	{
		std::string sum = "---------------------------------------------------";
		sum += "\n<ET::ScalarField<"+ type_name<decltype(m_field[0])>();
		sum += "> object at " + getMem(this) + ">";
		sum += "\n---------------------------------------------------";
    sum += "\n   name: '" + this->m_name + "'";
		sum += "\n    dim: " + std::to_string(this->m_dim);
		sum += "\n      N: " + std::to_string(this->m_N);
		sum += "\n---------------------------------------------------";
		std::string grid = "UGrid '" + this->m_ugrid->getName() + "' at: ";
		int grid_colon_loc = grid.length()-8;
		sum += "\n" + grid + getMem(this->m_log) + ",";
		std::string grid_ref;
		grid_ref.resize(grid_colon_loc,' ');
		sum += "\n" + grid_ref;
		sum += "ref at: " + getMem(this->m_ugrid);
		sum += "\nInterpolator at: " + getMem(*this->m_interpolator) + ",";
		sum += "\n         ref at: " + getMem(this->m_interpolator);
		sum += "\nLogger at: " + getMem(*this->m_log) + ",";
		sum += "\n   ref at: " + getMem(this->m_log);
		sum += "\n++++++++++++++++++++++++++++++++++++++++++++++++++++";
		return sum;
	}
	//----------------------------------------------------------------------------

}
