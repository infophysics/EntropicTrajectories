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
    this->m_Interpolator->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_DiffEQ->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_Integrator->setField(std::make_shared<ScalarField<T>>(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  ScalarField<T>::~ScalarField()
  {
  }
  //----------------------------------------------------------------------------
  template<typename T>
  ScalarField<T>::ScalarField(std::string t_name)
  : Field<T>(t_name)
  {
    this->m_dim = 1;
    this->m_Interpolator->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_DiffEQ->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_Integrator->setField(std::make_shared<ScalarField<T>>(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  ScalarField<T>::ScalarField(std::shared_ptr<Log> t_log)
  : Field<T>(t_log)
  {
    this->m_dim = 1;
    this->m_Interpolator->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_DiffEQ->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_Integrator->setField(std::make_shared<ScalarField<T>>(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  ScalarField<T>::ScalarField(std::string t_name, std::shared_ptr<Log> t_log)
  : Field<T>(t_name, t_log)
  {
    this->m_dim = 1;
    this->m_Interpolator->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_DiffEQ->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_Integrator->setField(std::make_shared<ScalarField<T>>(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  ScalarField<T>::ScalarField(std::shared_ptr<Grid<T>> t_grid)
  : Field<T>(t_grid)
  {
    this->m_dim = 1;
    this->m_Interpolator->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_DiffEQ->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_Integrator->setField(std::make_shared<ScalarField<T>>(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  ScalarField<T>::ScalarField(std::string t_name, std::shared_ptr<Grid<T>> t_grid)
  : Field<T>(t_name, t_grid)
  {
    this->m_dim = 1;
    this->m_Interpolator->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_DiffEQ->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_Integrator->setField(std::make_shared<ScalarField<T>>(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  ScalarField<T>::ScalarField(std::shared_ptr<Interpolator<T>> t_interpolator)
  : Field<T>(t_interpolator)
  {
    this->m_dim = 1;
    this->m_Interpolator->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_DiffEQ->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_Integrator->setField(std::make_shared<ScalarField<T>>(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  ScalarField<T>::ScalarField(std::string t_name,
                              std::shared_ptr<Interpolator<T>> t_interpolator)
  : Field<T>(t_name, t_interpolator)
  {
    this->m_dim = 1;
    this->m_Interpolator->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_DiffEQ->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_Integrator->setField(std::make_shared<ScalarField<T>>(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  ScalarField<T>::ScalarField(std::shared_ptr<Grid<T>> t_grid,
                              std::shared_ptr<Log> t_log)
  : Field<T>(t_grid, t_log)
  {
    this->m_dim = 1;
    this->m_Interpolator->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_DiffEQ->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_Integrator->setField(std::make_shared<ScalarField<T>>(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  ScalarField<T>::ScalarField(std::string t_name,
                              std::shared_ptr<Grid<T>> t_grid,
                              std::shared_ptr<Log> t_log)
  : Field<T>(t_name, t_grid, t_log)
  {
    this->m_dim = 1;
    this->m_Interpolator->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_DiffEQ->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_Integrator->setField(std::make_shared<ScalarField<T>>(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  ScalarField<T>::ScalarField(std::vector<T> t_field)
  : Field<T>()
  {
    m_field = t_field;
    this->m_dim = 1;
    this->m_N = m_field.size();
    this->m_Interpolator->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_DiffEQ->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_Integrator->setField(std::make_shared<ScalarField<T>>(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  ScalarField<T>::ScalarField(std::string t_name, std::vector<T> t_field)
  : Field<T>(t_name)
  {
    m_field = t_field;
    this->m_dim = 1;
    this->m_N = m_field.size();
    this->m_Interpolator->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_DiffEQ->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_Integrator->setField(std::make_shared<ScalarField<T>>(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  ScalarField<T>::ScalarField(std::vector<T> t_field,
                              std::shared_ptr<Log> t_log)
  : Field<T>(t_log)
  {
    m_field = t_field;
    this->m_dim = 1;
    this->m_N = m_field.size();
    this->m_Interpolator->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_DiffEQ->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_Integrator->setField(std::make_shared<ScalarField<T>>(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  ScalarField<T>::ScalarField(std::string t_name,
                              std::vector<T> t_field,
                              std::shared_ptr<Log> t_log)
  : Field<T>(t_name, t_log)
  {
    m_field = t_field;
    this->m_dim = 1;
    this->m_N = m_field.size();
    this->m_Interpolator->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_DiffEQ->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_Integrator->setField(std::make_shared<ScalarField<T>>(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  ScalarField<T>::ScalarField(std::vector<T> t_field,
                              std::shared_ptr<Grid<T>> t_grid)
  : Field<T>(t_grid)
  {
    m_field = t_field;
    this->m_dim = 1;
    this->m_N = m_field.size();
    this->m_Interpolator->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_DiffEQ->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_Integrator->setField(std::make_shared<ScalarField<T>>(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  ScalarField<T>::ScalarField(std::string t_name,
                              std::vector<T> t_field,
                              std::shared_ptr<Grid<T>> t_grid)
  : Field<T>(t_name, t_grid)
  {
    m_field = t_field;
    this->m_dim = 1;
    this->m_N = m_field.size();
    this->m_Interpolator->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_DiffEQ->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_Integrator->setField(std::make_shared<ScalarField<T>>(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  ScalarField<T>::ScalarField(std::vector<T> t_field,
                              std::shared_ptr<Interpolator<T>> t_interpolator)
  : Field<T>(t_interpolator)
  {
    m_field = t_field;
    this->m_dim = 1;
    this->m_N = m_field.size();
    this->m_Interpolator->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_DiffEQ->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_Integrator->setField(std::make_shared<ScalarField<T>>(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  ScalarField<T>::ScalarField(std::string t_name,
                              std::vector<T> t_field,
                              std::shared_ptr<Interpolator<T>> t_interpolator)
  : Field<T>(t_name, t_interpolator)
  {
    m_field = t_field;
    this->m_dim = 1;
    this->m_N = m_field.size();
    this->m_Interpolator->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_DiffEQ->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_Integrator->setField(std::make_shared<ScalarField<T>>(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  ScalarField<T>::ScalarField(std::vector<T> t_field,
                              std::shared_ptr<Grid<T>> t_grid,
                              std::shared_ptr<Log> t_log)
  : Field<T>(t_grid, t_log)
  {
    m_field = t_field;
    this->m_dim = 1;
    this->m_N = m_field.size();
    this->m_Interpolator->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_DiffEQ->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_Integrator->setField(std::make_shared<ScalarField<T>>(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  ScalarField<T>::ScalarField(std::string t_name,
                              std::vector<T> t_field,
                              std::shared_ptr<Grid<T>> t_grid,
                              std::shared_ptr<Log> t_log)
  : Field<T>(t_name, t_grid, t_log)
  {
    m_field = t_field;
    this->m_dim = 1;
    this->m_N = m_field.size();
    this->m_Interpolator->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_DiffEQ->setField(std::make_shared<ScalarField<T>>(*this));
    // this->m_Integrator->setField(std::make_shared<ScalarField<T>>(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  enum FieldType ScalarField<T>::getType() const
  {
    return m_type;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<T> ScalarField<T>::getField() const
  {
    return m_field;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<T>* ScalarField<T>::accessField()
  {
    return &m_field;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  T* ScalarField<T>::data()
  {
    return m_field.data();
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void ScalarField<T>::setField(std::vector<T> t_field)
  {
    m_field = t_field;
    this->m_N = m_field.size();
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void
  ScalarField<T>::setInterpolator(std::shared_ptr<Interpolator<T>> t_interpolator)
  {
    this->m_Interpolator = t_interpolator;
    //  m_Grid takes presidence over the grid from t_interpolator
    this->m_Interpolator->setGrid(this->m_Grid);
    this->m_Interpolator->setField(std::make_shared<ScalarField<T>>(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void
  ScalarField<T>::setDiffEQ(std::shared_ptr<DiffEQ<T>> t_diffeq)
  {
    this->m_DiffEQ = t_diffeq;
    //  m_Grid takes presidence over the grid from t_interpolator
    // this->m_DiffEQ->setGrid(this->m_Grid);
    // this->m_DiffEQ->setField(std::make_shared<ScalarField<T>>(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void
  ScalarField<T>::setIntegrator(std::shared_ptr<Integrator<T>> t_integrator)
  {
    this->m_Integrator = t_integrator;
    //  m_Grid takes presidence over the grid from t_interpolator
    // this->m_Integrator->setGrid(this->m_Grid);
    // this->m_Integrator->setField(std::make_shared<ScalarField<T>>(*this));
  }
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
		if (this->m_Grid != t_scalar.getGrid()) {
			this->m_log->WARN("Grids for t_scalar fields " + this->m_name  + " and "
			            + t_scalar.getName() + " do not match");
		}
		std::vector<T> copy = m_field;
		std::transform(copy.begin(), copy.end(), t_scalar.getField().begin(),
	                 copy.begin(), std::plus<T>());
		std::string name;
		if (this->m_name != " " && t_scalar.getName() != " ") {
			name += "(" + this->m_name + " + " + t_scalar.getName() + ")";
		}
		return ScalarField<T>(copy,this->m_Grid,this->m_log);
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
		if (this->m_Grid != t_scalar.getGrid()) {
			this->m_log->WARN("Grids for t_scalar fields " + this->m_name  + " and "
			            + t_scalar.getName() + " do not match");
		}
		std::vector<T> copy = m_field;
		std::transform(copy.begin(), copy.end(), t_scalar.getField().begin(),
	                 copy.begin(), std::minus<T>());
		std::string name;
		if (this->m_name != " " && t_scalar.getName() != " ") {
			name += "(" + this->m_name + " - " + t_scalar.getName() + ")";
		}
		return ScalarField<T>(copy,this->m_Grid,this->m_log);
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
		if (this->m_Grid != t_scalar.getGrid()) {
			this->m_log->WARN("Grids for t_scalar fields " + this->m_name  + " and "
			            + t_scalar.getName() + " do not match");
		}
		std::vector<T> copy = m_field;
		std::transform(copy.begin(), copy.end(), t_scalar.getField().begin(),
	                 copy.begin(), std::multiplies<T>());
		std::string name;
		if (this->m_name != " " && t_scalar.getName() != " ") {
			name += "(" + this->m_name + " * " + t_scalar.getName() + ")";
		}
		return ScalarField<T>(copy,this->m_Grid,this->m_log);
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
		if (this->m_Grid != t_scalar.getGrid()) {
			this->m_log->WARN("Grids for t_scalar fields " + this->m_name  + " and "
			            + t_scalar.getName() + " do not match");
		}
		std::vector<T> copy = m_field;
		std::transform(copy.begin(), copy.end(), t_scalar.getField().begin(),
	                 copy.begin(), std::divides<T>());
		std::string name;
		if (this->m_name != " " && t_scalar.getName() != " ") {
			name += "(" + this->m_name + " / " + t_scalar.getName() + ")";
		}
		return ScalarField<T>(copy,this->m_Grid,this->m_log);
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
		if (this->m_Grid != t_scalar.getGrid()) {
			this->m_log->WARN("Grids for t_scalar fields " + this->m_name  + " and "
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
		if (this->m_Grid != t_scalar.getGrid()) {
			this->m_log->WARN("Grids for t_scalar fields " + this->m_name  + " and "
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
		if (this->m_Grid != t_scalar.getGrid()) {
			this->m_log->WARN("Grids for t_scalar fields " + this->m_name  + " and "
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
		if (this->m_Grid != t_scalar.getGrid()) {
			this->m_log->WARN("Grids for t_scalar fields " + this->m_name  + " and "
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
    this->m_Grid->getKDTree()->queryNeighbors(this->m_Grid->getKDTree()->getCurrentGlobalK());
    std::vector<size_t>
    neighbors = this->m_Grid->getKDTree()->getCurrentNeighborIndices(t_index);
    //  Create the empty vector
    Vector<T> f(neighbors.size());
    for (auto i = 0; i < neighbors.size(); i++) {
      f(i) = m_field[neighbors[i]];
    }
    return f;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Vector<T>
  ScalarField<T>::constructLocalFieldValues(const std::vector<T>& t_point,
                                                  size_t t_k)
  {
    //  Get the nearest neighbors for the index
    std::vector<size_t>
    neighbors = this->m_Grid->getKDTree()->queryNeighbors(t_point, t_k);
    //  Create the empty vector
    Vector<T> f(neighbors.size());
    for (auto i = 0; i < neighbors.size(); i++) {
      f(i) = m_field[neighbors[i]];
    }
    return f;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Vector<T> ScalarField<T>::derivative(const size_t t_index,
                                       const size_t t_degree)
  {
    return this->m_Interpolator->derivative(t_index, t_degree);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  T ScalarField<T>::derivative(const size_t t_index,
                               const size_t t_degree,
                               const size_t t_direction)
  {
    return this->m_Interpolator->derivative(t_index, t_degree, t_direction);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Vector<T> ScalarField<T>::derivative(const std::vector<T>& t_point,
                                       const size_t t_degree)
  {
    return this->m_Interpolator->derivative(t_point, t_degree);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  T ScalarField<T>::derivative(const std::vector<T>& t_point,
                                       const size_t t_degree,
                                       const size_t t_direction)
  {
    return this->m_Interpolator->derivative(t_point, t_degree, t_direction);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<Vector<T>>
  ScalarField<T>::fieldDerivative(const size_t t_degree)
  {
    return this->m_Interpolator->fieldDerivative(t_degree);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<T> ScalarField<T>::fieldDerivative(const size_t t_degree,
                                    const size_t t_direction)
  {
    return this->m_Interpolator->fieldDerivative(t_degree, t_direction);
  }
  //----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Various functions
	//----------------------------------------------------------------------------
	template<typename T>
	const std::string ScalarField<T>::summary()
	{
		std::string sum = "---------------------------------------------------";
		sum += "\n<ET::ScalarField<"+ type_name<decltype(m_field[0])>();
		sum += "> object at " + address_to_string(this) + ">";
		sum += "\n---------------------------------------------------";
    sum += "\n   name: '" + this->m_name + "'";
		sum += "\n    dim: " + std::to_string(this->m_dim);
		sum += "\n      N: " + std::to_string(this->m_N);
		sum += "\n---------------------------------------------------";
		std::string grid = "Grid '" + this->m_Grid->getName() + "' at: ";
		int grid_colon_loc = grid.length()-8;
		sum += "\n" + grid + address_to_string(this->m_log) + ",";
		std::string grid_ref;
		grid_ref.resize(grid_colon_loc,' ');
		sum += "\n" + grid_ref;
		sum += "ref at: " + address_to_string(this->m_Grid);
		sum += "\nInterpolator at: " + address_to_string(*this->m_Interpolator) + ",";
		sum += "\n         ref at: " + address_to_string(this->m_Interpolator);
		sum += "\nLogger at: " + address_to_string(*this->m_log) + ",";
		sum += "\n   ref at: " + address_to_string(this->m_log);
		sum += "\n++++++++++++++++++++++++++++++++++++++++++++++++++++";
		return sum;
	}
	//----------------------------------------------------------------------------

}
