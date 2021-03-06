//------------------------------------------------------------------------------
//  field.cpp
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
#include "field.h"

namespace ET
{
  //----------------------------------------------------------------------------
  std::map<std::string, BoundaryConditionType> BoundaryConditionTypeMap =
  {
    { "Dirichlet", BoundaryConditionType::DIRICHLET },
    { "Neumann",   BoundaryConditionType::NEUMANN },
    { "Robin",     BoundaryConditionType::ROBIN },
    { "Mixed",     BoundaryConditionType::MIXED },
    { "Cauchy",    BoundaryConditionType::CAUCHY },
  };
  //----------------------------------------------------------------------------
  std::map<BoundaryConditionType, std::string> BoundaryConditionTypeNameMap =
  {
    { BoundaryConditionType::DIRICHLET, "Dirichlet" },
    { BoundaryConditionType::NEUMANN,   "Neumann" },
    { BoundaryConditionType::ROBIN,     "Robin" },
    { BoundaryConditionType::MIXED,     "Mixed" },
    { BoundaryConditionType::CAUCHY,    "Cauchy" },
  };
  //----------------------------------------------------------------------------
  template<typename T>
  Field<T>::Field()
  : m_name("default"), m_dim(0), m_N(0)
  {
    m_log = std::make_shared<Log>();
		m_log->init(NAME(), ".logs/field_" + m_name + ".txt");
		m_log->TRACE(NAME() + "Field '" + m_name + "' created at location "
		            + address_to_string(*this));
    m_Grid = std::make_shared<Grid<T>>(m_log);
    m_Interpolator = std::make_shared<Interpolator<T>>(m_Grid,m_log);
    m_DiffEQ = std::make_shared<DiffEQ<T>>(m_Grid,m_log);
    m_Integrator = std::make_shared<Integrator<T>>(m_Grid,m_log);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Field<T>::~Field()
  {
    m_log->TRACE(NAME() + "Field '" + m_name
								+ "' destroyed at location " + address_to_string(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Field<T>::Field(std::string t_name)
  : m_name(t_name), m_dim(0), m_N(0)
  {
    m_log = std::make_shared<Log>();
		m_log->init(NAME(), ".logs/field_" + m_name + ".txt");
		m_log->TRACE(NAME() + "Field '" + m_name + "' created at location "
		            + address_to_string(*this));
    m_Grid = std::make_shared<Grid<T>>(m_log);
    m_Interpolator = std::make_shared<Interpolator<T>>(m_Grid,m_log);
    m_DiffEQ = std::make_shared<DiffEQ<T>>(m_Grid,m_log);
    m_Integrator = std::make_shared<Integrator<T>>(m_Grid,m_log);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Field<T>::Field(std::shared_ptr<Log> t_log)
  : m_name("default"), m_dim(0), m_N(0)
  {
    m_log = t_log;
		m_log->TRACE(NAME() + "Field '" + m_name + "' created at location "
		            + address_to_string(*this));
		m_log->INFO(NAME() + "Logger passed to Field 'default'");
    m_Grid = std::make_shared<Grid<T>>(m_log);
    m_Interpolator = std::make_shared<Interpolator<T>>(m_Grid,m_log);
    m_DiffEQ = std::make_shared<DiffEQ<T>>(m_Grid,m_log);
    m_Integrator = std::make_shared<Integrator<T>>(m_Grid,m_log);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Field<T>::Field(std::string t_name, std::shared_ptr<Log> t_log)
  : m_name(t_name), m_dim(0), m_N(0)
  {
    m_log = t_log;
		m_log->TRACE(NAME() + "Field '" + m_name + "' created at location "
		            + address_to_string(*this));
		m_log->INFO(NAME() + "Logger passed to Field 'default'");
    m_Grid = std::make_shared<Grid<T>>(m_log);
    m_Interpolator = std::make_shared<Interpolator<T>>(m_Grid,m_log);
    m_DiffEQ = std::make_shared<DiffEQ<T>>(m_Grid,m_log);
    m_Integrator = std::make_shared<Integrator<T>>(m_Grid,m_log);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Field<T>::Field(std::shared_ptr<Grid<T>> t_grid)
  : m_name("default"), m_dim(t_grid->getDim()), m_N(t_grid->getN())
  {
    m_log = std::make_shared<Log>();
    m_Grid = t_grid;
		m_log->init(NAME(), ".logs/field_" + m_name + ".txt");
		m_log->TRACE(NAME() + "Field '" + m_name + "' created at location "
		            + address_to_string(*this));
    m_Interpolator = std::make_shared<Interpolator<T>>(m_Grid,m_log);
    m_DiffEQ = std::make_shared<DiffEQ<T>>(m_Grid,m_log);
    m_Integrator = std::make_shared<Integrator<T>>(m_Grid,m_log);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Field<T>::Field(std::string t_name, std::shared_ptr<Grid<T>> t_grid)
  : m_name(t_name), m_dim(t_grid->getDim()), m_N(t_grid->getN())
  {
    m_log = std::make_shared<Log>();
    m_Grid = t_grid;
		m_log->init(NAME(), ".logs/field_" + m_name + ".txt");
		m_log->TRACE(NAME() + "Field '" + m_name + "' created at location "
		            + address_to_string(*this));
    m_Interpolator = std::make_shared<Interpolator<T>>(m_Grid,m_log);
    m_DiffEQ = std::make_shared<DiffEQ<T>>(m_Grid,m_log);
    m_Integrator = std::make_shared<Integrator<T>>(m_Grid,m_log);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Field<T>::Field(std::shared_ptr<Interpolator<T>> t_interpolator)
  : m_name("default"), m_dim(t_interpolator->getGrid()->getDim()),
    m_N(t_interpolator->getGrid()->getN())
  {
    m_log = std::make_shared<Log>();
    m_Interpolator = t_interpolator;
		m_log->init(NAME(), ".logs/field_" + m_name + ".txt");
		m_log->TRACE(NAME() + "Field '" + m_name + "' created at location "
		            + address_to_string(*this));
    m_Grid = m_Interpolator->getGrid();
    m_DiffEQ = std::make_shared<DiffEQ<T>>(m_Grid,m_log);
    m_Integrator = std::make_shared<Integrator<T>>(m_Grid,m_log);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Field<T>::Field(std::string t_name,
                  std::shared_ptr<Interpolator<T>> t_interpolator)
  : m_name(t_name), m_dim(t_interpolator->getGrid()->getDim()),
    m_N(t_interpolator->getGrid()->getN())
  {
    m_log = std::make_shared<Log>();
    m_Interpolator = t_interpolator;
		m_log->init(NAME(), ".logs/field_" + m_name + ".txt");
		m_log->TRACE(NAME() + "Field '" + m_name + "' created at location "
		            + address_to_string(*this));
    m_Grid = m_Interpolator->getGrid();
    m_DiffEQ = std::make_shared<DiffEQ<T>>(m_Grid,m_log);
    m_Integrator = std::make_shared<Integrator<T>>(m_Grid,m_log);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Field<T>::Field(std::shared_ptr<Grid<T>> t_grid, std::shared_ptr<Log> t_log)
  : m_name("default"), m_dim(t_grid->getDim()), m_N(t_grid->getN())
  {
    m_log = t_log;
    m_Grid = t_grid;
		m_log->TRACE(NAME() + "Field '" + m_name + "' created at location "
		            + address_to_string(*this));
		m_log->INFO(NAME() + "Logger passed to Field 'default'");
    m_Interpolator = std::make_shared<Interpolator<T>>(m_Grid,m_log);
    m_DiffEQ = std::make_shared<DiffEQ<T>>(m_Grid,m_log);
    m_Integrator = std::make_shared<Integrator<T>>(m_Grid,m_log);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Field<T>::Field(std::string t_name, std::shared_ptr<Grid<T>> t_grid,
                  std::shared_ptr<Log> t_log)
  : m_name(t_name), m_dim(t_grid->getDim()), m_N(t_grid->getN())
  {
    m_log = t_log;
    m_Grid = t_grid;
    m_log->TRACE(NAME() + "Field '" + m_name + "' created at location "
                + address_to_string(*this));
    m_log->INFO(NAME() + "Logger passed to Field 'default'");
    m_Interpolator = std::make_shared<Interpolator<T>>(m_Grid,m_log);
    m_DiffEQ = std::make_shared<DiffEQ<T>>(m_Grid,m_log);
    m_Integrator = std::make_shared<Integrator<T>>(m_Grid,m_log);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::string Field<T>::getName() const
  {
    return m_name;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  size_t Field<T>::getDim() const
  {
    return m_dim;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  size_t Field<T>::getN() const
  {
    return m_N;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::shared_ptr<Log> Field<T>::getLog() const
  {
    return m_log;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::shared_ptr<Grid<T>> Field<T>::getGrid() const
  {
    return m_Grid;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::shared_ptr<Interpolator<T>> Field<T>::getInterpolator() const
  {
    return m_Interpolator;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::shared_ptr<DiffEQ<T>> Field<T>::getDiffEQ() const
  {
    return m_DiffEQ;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::shared_ptr<Integrator<T>> Field<T>::getIntegrator() const
  {
    return m_Integrator;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<struct BoundaryCondition<T>> Field<T>::getBoundaryConditions() const
  {
    return m_boundary_conditions;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  struct BoundaryCondition<T> Field<T>::getBoundaryCondition(const size_t t_i) const
  {
    if (t_i >= m_boundary_conditions.size()) {
      m_log->ERROR(NAME() + "Attempted to access boundary condition at index "
                   + std::to_string(t_i) + ", while there are only "
                   + std::to_string(m_boundary_conditions.size()));
      m_log->TRACE("Returning empty BoundaryCondition");
      return BoundaryCondition<T>();
    }
    else {
      return m_boundary_conditions[t_i];
    }
  }
  //----------------------------------------------------------------------------
  template<typename T>
  int Field<T>::getFlag() const
  {
    return m_flag;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::string Field<T>::getInfo() const
  {
    return m_info;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  enum FieldType Field<T>::getType() const
  {
    return m_type;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Field<T>::setName(std::string t_name)
  {
    m_log->TRACE(NAME() + "Changed name from '" + m_name
                 + "' to '" + t_name + "'");
    m_name = t_name;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Field<T>::setDim(size_t t_dim)
  {
    m_log->TRACE(NAME() + "Changed dim from '" + std::to_string(m_dim)
                 + "' to '" + std::to_string(t_dim) + "'");
    m_dim = t_dim;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Field<T>::setN(size_t t_N)
  {
    m_log->TRACE(NAME() + "Changed N from '" + std::to_string(m_N)
                 + "' to '" + std::to_string(t_N) + "'");
    m_N = t_N;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Field<T>::setGrid(std::shared_ptr<Grid<T>> t_grid)
  {
    m_log->TRACE(NAME() + "Changed Grid from '" + address_to_string(m_Grid)
                 + "' to '" + address_to_string(t_grid) + "'");
    m_Grid = t_grid;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Field<T>::setLog(std::shared_ptr<Log> t_log)
  {
    m_log->TRACE(NAME() + "Changed log from '" + address_to_string(m_log)
                 + "' to '" + address_to_string(t_log) + "'");
    m_log = t_log;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Field<T>::setInterpolator(std::shared_ptr<Interpolator<T>> t_interpolator)
  {
    m_log->TRACE(NAME() + "Changed Interpolator from '"
                 + address_to_string(m_Interpolator)
                 + "' to '" + address_to_string(t_interpolator) + "'");
    m_Interpolator = t_interpolator;
    //  m_Grid takes presidence over the grid from t_interpolator
    m_Interpolator->setGrid(m_Grid);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Field<T>::setDiffEQ(std::shared_ptr<DiffEQ<T>> t_diffeq)
  {
    m_log->TRACE(NAME() + "Changed DiffEQ from '" + address_to_string(m_DiffEQ)
                 + "' to '" + address_to_string(t_diffeq) + "'");
    m_DiffEQ = t_diffeq;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Field<T>::setIntegrator(std::shared_ptr<Integrator<T>> t_integrator)
  {
    m_log->TRACE(NAME() + "Changed Integrator from '"
                 + address_to_string(m_Integrator)
                 + "' to '" + address_to_string(t_integrator) + "'");
    m_Integrator = t_integrator;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void
  Field<T>::setBoundaryConditions(std::vector<struct BoundaryCondition<T>>
                                  t_boundaryConditions)
  {
    m_boundary_conditions = t_boundaryConditions;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Field<T>::setBoundaryCondition(const size_t t_i,
                                      struct BoundaryCondition<T>
                                      t_boundaryCondition)
  {
    if (t_i >= m_boundary_conditions.size()) {
      m_log->ERROR(NAME() + "Attempted to access boundary condition at index "
                   + std::to_string(t_i) + ", while there are only "
                   + std::to_string(m_boundary_conditions.size()));
      return;
    }
    else {
      m_boundary_conditions[t_i] = t_boundaryCondition;
    }
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Field<T>::addBoundaryCondition(struct BoundaryCondition<T>
                                      t_boundaryCondition)
  {
    m_log->TRACE(NAME() + "Adding boundary condition of type "
                 + BoundaryConditionTypeNameMap[t_boundaryCondition.m_type]
                 + " to list of Boundary Conditions.");
    m_boundary_conditions.push_back(t_boundaryCondition);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Field<T>::setFlag(int t_flag)
  {
    m_log->TRACE(NAME() + "Changed flag from '" + std::to_string(m_flag)
                 + "' to '" + std::to_string(t_flag) + "'");
    m_flag = t_flag;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Field<T>::setInfo(std::string t_info)
  {
    m_log->TRACE(NAME() + "Changed info from '" + m_info
                 + "' to '" + t_info + "'");
    m_info = t_info;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Vector<T> Field<T>::constructLocalFieldValues(size_t t_index)
  {
    m_log->WARN(NAME() + "Running 'constructLocalFieldValues' from Base class");
    m_log->WARN(NAME() + "Returning empty Vector<T>");
    return Vector<T>();
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Vector<T> Field<T>::constructLocalFieldValues(const std::vector<T>& t_point,
                                              size_t t_k)
  {
    m_log->WARN(NAME() + "Running 'constructLocalFieldValues' from Base class");
    m_log->WARN(NAME() + "Returning empty Vector<T>");
    return Vector<T>();
  }
  //----------------------------------------------------------------------------
}
