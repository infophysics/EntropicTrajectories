//------------------------------------------------------------------------------
//  integrator.cpp
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
#include "integrator.h"

namespace ET
{
  //----------------------------------------------------------------------------
  //  Base class Integrator
  //----------------------------------------------------------------------------
  template<typename T>
  Integrator<T>::Integrator()
  : m_name("default"), m_dim(0)
  {
    m_log = std::make_shared<Log>();
		m_log->init(NAME(), ".logs/integrator_" + m_name + ".txt");
		m_log->TRACE(NAME() + "Integrator '" + m_name + "' created at location "
		            + address_to_string(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Integrator<T>::~Integrator()
  {
    m_log->TRACE(NAME() + "Integrator '" + m_name
								+ "' destroyed at location " + address_to_string(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Integrator<T>::Integrator(std::shared_ptr<Log> t_log)
  : m_name("default"), m_dim(0), m_log(t_log)
  {
    m_log->TRACE(NAME() + "Integrator '" + m_name + "' created at location "
                + address_to_string(*this));
    m_log->INFO(NAME() + "Logger passed to Integrator '" + m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Integrator<T>::Integrator(std::shared_ptr<Grid<T>> t_Grid, std::shared_ptr<Log> t_log)
  : m_name("default"), m_dim(0), m_Grid(t_Grid), m_log(t_log)
  {
    m_log->TRACE(NAME() + "Integrator '" + m_name + "' created at location "
                + address_to_string(*this));
    m_log->INFO(NAME() + "Logger passed to Integrator '" + m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::string Integrator<T>::getName() const
  {
    return m_name;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  size_t Integrator<T>::getDim() const
  {
    return m_dim;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::shared_ptr<Grid<T>> Integrator<T>::getGrid() const
  {
    return m_Grid;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::shared_ptr<DiffEQ<T>> Integrator<T>::getDiffEQ() const
  {
    return m_DiffEQ;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::shared_ptr<Log> Integrator<T>::getLog() const
  {
    return m_log;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  enum IntegratorType Integrator<T>::getType() const
  {
    return m_type;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Integrator<T>::setName(const std::string t_name)
  {
    m_log->TRACE(NAME() + "Changed name from '" + m_name
                 + "' to '" + t_name + "'");
    m_name = t_name;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Integrator<T>::setDim(const size_t t_dim)
  {
    m_log->TRACE(NAME() + "Changed dim from '" + std::to_string(m_dim)
                 + "' to '" + std::to_string(t_dim) + "'");
    m_dim = t_dim;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Integrator<T>::setGrid(const std::shared_ptr<Grid<T>> t_Grid)
  {
    m_Grid = t_Grid;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Integrator<T>::setDiffEQ(const std::shared_ptr<DiffEQ<T>> t_DiffEQ)
  {
    m_DiffEQ = t_DiffEQ;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Integrator<T>::setLog(const std::shared_ptr<Log> t_log)
  {
    m_log = t_log;
  }
  //----------------------------------------------------------------------------
  //  Base class ODEIntegrator
  //----------------------------------------------------------------------------
  template<typename T>
  ODEIntegrator<T>::ODEIntegrator()
  : DiffEQ<T>()
  {
    this->m_log = std::make_shared<Log>();
		this->m_log->init(NAME(), ".logs/odeintegrator_" + this->m_name + ".txt");
		this->m_log->TRACE(NAME() + "ODEIntegrator '" + this->m_name
                       + "' created at location " + address_to_string(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  ODEIntegrator<T>::~ODEIntegrator()
  {
    this->m_log->TRACE(NAME() + "ODEIntegrator '" + this->m_name
								+ "' destroyed at location " + address_to_string(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  ODEIntegrator<T>::ODEIntegrator(std::shared_ptr<Log> t_log)
  : DiffEQ<T>(t_log)
  {
    this->m_log->TRACE(NAME() + "ODEIntegrator '" + this->m_name
                       + "' created at location " + address_to_string(*this));
    this->m_log->INFO(NAME() + "Logger passed to ODEIntegrator '" + this->m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  ODEIntegrator<T>::ODEIntegrator(std::shared_ptr<Grid<T>> t_Grid,
                                  std::shared_ptr<Log> t_log)
  : DiffEQ<T>(t_Grid,t_log)
  {
    this->m_log->TRACE(NAME() + "ODEIntegrator '" + this->m_name
                       + "' created at location " + address_to_string(*this));
    this->m_log->INFO(NAME() + "Logger passed to ODEIntegrator '" + this->m_name + "'");
  }
  //----------------------------------------------------------------------------
  //  Base class ODEIntegrator
  //----------------------------------------------------------------------------
  template<typename T>
  PDEIntegrator<T>::PDEIntegrator()
  : DiffEQ<T>()
  {
    this->m_log = std::make_shared<Log>();
		this->m_log->init(NAME(), ".logs/odeintegrator_" + this->m_name + ".txt");
		this->m_log->TRACE(NAME() + "PDEIntegrator '" + this->m_name
                       + "' created at location " + address_to_string(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  PDEIntegrator<T>::~PDEIntegrator()
  {
    this->m_log->TRACE(NAME() + "PDEIntegrator '" + this->m_name
								+ "' destroyed at location " + address_to_string(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  PDEIntegrator<T>::PDEIntegrator(std::shared_ptr<Log> t_log)
  : DiffEQ<T>(t_log)
  {
    this->m_log->TRACE(NAME() + "PDEIntegrator '" + this->m_name
                       + "' created at location " + address_to_string(*this));
    this->m_log->INFO(NAME() + "Logger passed to PDEIntegrator '" + this->m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  PDEIntegrator<T>::PDEIntegrator(std::shared_ptr<Grid<T>> t_Grid,
                                  std::shared_ptr<Log> t_log)
  : DiffEQ<T>(t_Grid,t_log)
  {
    this->m_log->TRACE(NAME() + "PDEIntegrator '" + this->m_name
                       + "' created at location " + address_to_string(*this));
    this->m_log->INFO(NAME() + "Logger passed to PDEIntegrator '" + this->m_name + "'");
  }
  //----------------------------------------------------------------------------
}
