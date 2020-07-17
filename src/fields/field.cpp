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
  template<typename T>
  Field<T>::Field()
  : m_name("default"), m_dim(0), m_N(0)
  {
    m_log = std::make_shared<Log>();
		m_log->init("ET:Field:default", ".logs/field_default.txt");
		m_log->TRACE("Field 'default' created at location "
		            + getMem(*this));
    m_ugrid = std::make_shared<UGrid<T>>(m_log);
    m_interpolator = std::make_shared<Interpolator<T>>(m_ugrid,m_log);
    m_diffeq = std::make_shared<DiffEQ<T>>(m_ugrid,m_log);
    m_integrator = std::make_shared<Integrator<T>>(m_ugrid,m_log);
  }
  template<typename T>
  Field<T>::~Field()
  {
    m_log->TRACE("Field '" + m_name
								+ "' destroyed at location " + getMem(*this));
  }
  template<typename T>
  Field<T>::Field(std::shared_ptr<Log> t_log)
  {
    m_log = t_log;
		m_log->TRACE("Field 'default' created at location "
		            + getMem(*this));
		m_log->INFO("Logger passed to Field 'default'");
    m_ugrid = std::make_shared<UGrid<T>>(m_log);
    m_interpolator = std::make_shared<Interpolator<T>>(m_ugrid,m_log);
    m_diffeq = std::make_shared<DiffEQ<T>>(m_ugrid,m_log);
    m_integrator = std::make_shared<Integrator<T>>(m_ugrid,m_log);
  }
  template<typename T>
  Field<T>::Field(std::shared_ptr<UGrid<T>> t_ugrid)
  {
    m_log = std::make_shared<Log>();
    m_ugrid = t_ugrid;
		m_log->init("ET:Field:default", ".logs/field_default.txt");
		m_log->TRACE("Field 'default' created at location "
		            + getMem(*this));
    m_interpolator = std::make_shared<Interpolator<T>>(m_ugrid,m_log);
    m_diffeq = std::make_shared<DiffEQ<T>>(m_ugrid,m_log);
    m_integrator = std::make_shared<Integrator<T>>(m_ugrid,m_log);
  }
  template<typename T>
  Field<T>::Field(std::shared_ptr<Interpolator<T>> t_interpolator)
  {
    m_log = std::make_shared<Log>();
    m_interpolator = t_interpolator;
		m_log->init("ET:Field:default", ".logs/field_default.txt");
		m_log->TRACE("Field 'default' created at location "
		            + getMem(*this));
    m_ugrid = m_interpolator->getUGrid();
    m_diffeq = std::make_shared<DiffEQ<T>>(m_ugrid,m_log);
    m_integrator = std::make_shared<Integrator<T>>(m_ugrid,m_log);
  }
  template<typename T>
  Field<T>::Field(std::shared_ptr<UGrid<T>> t_ugrid, std::shared_ptr<Log> t_log)
  {
    m_log = t_log;
    m_ugrid = t_ugrid;
		m_log->TRACE("Field 'default' created at location "
		            + getMem(*this));
		m_log->INFO("Logger passed to Field 'default'");
    m_interpolator = std::make_shared<Interpolator<T>>(m_ugrid,m_log);
    m_diffeq = std::make_shared<DiffEQ<T>>(m_ugrid,m_log);
    m_integrator = std::make_shared<Integrator<T>>(m_ugrid,m_log);
  }

  template<typename T>
  std::string Field<T>::getName() const
  {
    return m_name;
  }
  template<typename T>
  size_t Field<T>::getDim() const
  {
    return m_dim;
  }
  template<typename T>
  size_t Field<T>::getN() const
  {
    return m_N;
  }
  template<typename T>
  std::shared_ptr<Log> Field<T>::getLog() const
  {
    return m_log;
  }
  template<typename T>
  std::shared_ptr<UGrid<T>> Field<T>::getUGrid() const
  {
    return m_ugrid;
  }
  template<typename T>
  std::shared_ptr<Interpolator<T>> Field<T>::getInterpolator() const
  {
    return m_interpolator;
  }
  template<typename T>
  std::shared_ptr<DiffEQ<T>> Field<T>::getDiffEQ() const
  {
    return m_diffeq;
  }
  template<typename T>
  std::shared_ptr<Integrator<T>> Field<T>::getIntegrator() const
  {
    return m_integrator;
  }
  template<typename T>
  int Field<T>::getFlag() const
  {
    return m_flag;
  }
  template<typename T>
  std::string Field<T>::getInfo() const
  {
    return m_info;
  }
  template<typename T>
  void Field<T>::setName(std::string t_name)
  {
    m_name = t_name;
  }
  template<typename T>
  void Field<T>::setUGrid(std::shared_ptr<UGrid<T>> t_ugrid)
  {
    m_ugrid = t_ugrid;
  }
  template<typename T>
  void Field<T>::setLog(std::shared_ptr<Log> t_log)
  {
    m_log = t_log;
  }
  template<typename T>
  void Field<T>::setInterpolator(std::shared_ptr<Interpolator<T>> t_interpolator)
  {
    m_interpolator = t_interpolator;
  }
  template<typename T>
  void Field<T>::setDiffEQ(std::shared_ptr<DiffEQ<T>> t_diffeq)
  {
    m_diffeq = t_diffeq;
  }
  template<typename T>
  void Field<T>::setIntegrator(std::shared_ptr<Integrator<T>> t_integrator)
  {
    m_integrator = t_integrator;
  }
  template<typename T>
  void Field<T>::setFlag(int t_flag)
  {
    m_flag = t_flag;
  }
  template<typename T>
  void Field<T>::setInfo(std::string t_info)
  {
    m_info = t_info;
  }
  template<typename T>
  Vector<T> Field<T>::constructLocalFieldValues(size_t t_index)
  {
    return Vector<T>();
  }
  template<typename T>
  Vector<T> Field<T>::constructLocalFieldValues(const std::vector<T>& t_point,
                                              size_t t_k)
  {
    return Vector<T>();
  }
}
