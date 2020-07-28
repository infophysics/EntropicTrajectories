//------------------------------------------------------------------------------
//  diffeq.cpp
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
#include "diffeq.h"

namespace ET
{
  //----------------------------------------------------------------------------
  //  Base class DiffEQ
  //----------------------------------------------------------------------------
  template<typename T>
  DiffEQ<T>::DiffEQ()
  : m_name("default"), m_dim(0)
  {
    m_log = std::make_shared<Log>();
		m_log->init(NAME(), ".logs/diffeq_" + m_name + ".txt");
		m_log->TRACE(NAME() + "DiffEQ '" + m_name + "' created at location "
		            + address_to_string(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  DiffEQ<T>::~DiffEQ()
  {
    m_log->TRACE(NAME() + "DiffEQ '" + m_name
								+ "' destroyed at location " + address_to_string(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  DiffEQ<T>::DiffEQ(std::shared_ptr<Log> t_log)
  : m_name("default"), m_dim(0), m_log(t_log)
  {
    m_log->TRACE(NAME() + "DiffEQ '" + m_name + "' created at location "
                + address_to_string(*this));
    m_log->INFO(NAME() + "Logger passed to DiffEQ '" + m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  DiffEQ<T>::DiffEQ(std::shared_ptr<Grid<T>> t_Grid, std::shared_ptr<Log> t_log)
  : m_name("default"), m_dim(0), m_Grid(t_Grid), m_log(t_log)
  {
    m_log->TRACE(NAME() + "DiffEQ '" + m_name + "' created at location "
                + address_to_string(*this));
    m_log->INFO(NAME() + "Logger passed to DiffEQ '" + m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  DiffEQ<T>::DiffEQ(std::shared_ptr<Grid<T>> t_Grid,
                    std::shared_ptr<Field<T>> t_Field,
                    std::shared_ptr<Log> t_log)
  : m_name("default"), m_Grid(t_Grid), m_log(t_log), m_Field(t_Field)
  {
    m_log->TRACE(NAME() + "DiffEQ '" + m_name + "' created at location "
                + address_to_string(*this));
    m_log->INFO(NAME() + "Log passed to DiffEQ '" + m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  DiffEQ<T>::DiffEQ(std::string t_name,
                    std::shared_ptr<Grid<T>> t_Grid,
                    std::shared_ptr<Field<T>> t_Field,
                    std::shared_ptr<Log> t_log)
  : m_name(t_name), m_Grid(t_Grid), m_log(t_log), m_Field(t_Field)
  {
    m_log->TRACE(NAME() + "DiffEQ '" + m_name + "' created at location "
                + address_to_string(*this));
    m_log->INFO(NAME() + "Log passed to DiffEQ '" + m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::string DiffEQ<T>::getName() const
  {
    return m_name;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  size_t DiffEQ<T>::getDim() const
  {
    return m_dim;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::shared_ptr<Grid<T>> DiffEQ<T>::getGrid() const
  {
    return m_Grid;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::shared_ptr<Field<T>> DiffEQ<T>::getField() const
  {
    return m_Field;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::shared_ptr<Log> DiffEQ<T>::getLog() const
  {
    return m_log;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  enum DiffEQType DiffEQ<T>::getType() const
  {
    return m_type;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void DiffEQ<T>::setName(const std::string t_name)
  {
    m_log->TRACE(NAME() + "Changed name from '" + m_name
                 + "' to '" + t_name + "'");
    m_name = t_name;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void DiffEQ<T>::setDim(const size_t t_dim)
  {
    m_log->TRACE(NAME() + "Changed dim from '" + std::to_string(m_dim)
                 + "' to '" + std::to_string(t_dim) + "'");
    m_dim = t_dim;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void DiffEQ<T>::setGrid(const std::shared_ptr<Grid<T>> t_Grid)
  {
    m_Grid = t_Grid;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void DiffEQ<T>::setField(const std::shared_ptr<Field<T>> t_Field)
  {
    m_Field = t_Field;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void DiffEQ<T>::setLog(const std::shared_ptr<Log> t_log)
  {
    m_log = t_log;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  T DiffEQ<T>::evaluate(const std::vector<T>& t_point, double t_time)
  {
    m_log->WARN(NAME() + "Called DEFAULT evaluate function for DiffEQ");
    return 0;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<T> DiffEQ<T>::evaluate(const std::vector<std::vector<T>>& t_points,
                                     double t_time)
  {
    m_log->WARN(NAME() + "Called DEFAULT evaluate function for DiffEQ");
   return std::vector<T>();
  }
  //----------------------------------------------------------------------------
  //  Base class FirstOrderODE
  //----------------------------------------------------------------------------
  template<typename T>
  FirstOrderODE<T>::FirstOrderODE()
  : DiffEQ<T>()
  {
    this->m_log = std::make_shared<Log>();
		this->m_log->init(NAME(), ".logs/diffeq_" + this->m_name + ".txt");
		this->m_log->TRACE(NAME() + "FirstOrderODE '" + this->m_name
                       + "' created at location " + address_to_string(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  FirstOrderODE<T>::~FirstOrderODE()
  {
    this->m_log->TRACE(NAME() + "FirstOrderODE '" + this->m_name
								+ "' destroyed at location " + address_to_string(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  FirstOrderODE<T>::FirstOrderODE(std::shared_ptr<Log> t_log)
  : DiffEQ<T>(t_log)
  {
    this->m_log->TRACE(NAME() + "FirstOrderODE '" + this->m_name
                       + "' created at location " + address_to_string(*this));
    this->m_log->INFO(NAME() + "Logger passed to FirstOrderODE '" + this->m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  FirstOrderODE<T>::FirstOrderODE(std::shared_ptr<Grid<T>> t_Grid,
                                  std::shared_ptr<Log> t_log)
  : DiffEQ<T>(t_Grid,t_log)
  {
    this->m_log->TRACE(NAME() + "FirstOrderODE '" + this->m_name
                       + "' created at location " + address_to_string(*this));
    this->m_log->INFO(NAME() + "Logger passed to FirstOrderODE '" + this->m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  T FirstOrderODE<T>::evaluate(const std::vector<T>& t_point, double t_time)
  {
    this->m_log->WARN(NAME() + "Called DEFAULT evaluate function for First Order ODE");
    return 0;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<T>
  FirstOrderODE<T>::evaluate(const std::vector<std::vector<T>>& t_points,
                             double t_time)
  {
    this->m_log->WARN(NAME() + "Called DEFAULT evaluate function for First Order ODE");
    return std::vector<T>();
  }
  //----------------------------------------------------------------------------
  //  Base class SecondOrderODE
  //----------------------------------------------------------------------------
  template<typename T>
  SecondOrderODE<T>::SecondOrderODE()
  : DiffEQ<T>()
  {
    this->m_log = std::make_shared<Log>();
		this->m_log->init(NAME(), ".logs/diffeq_" + this->m_name + ".txt");
		this->m_log->TRACE(NAME() + "SecondOrderODE '" + this->m_name
                       + "' created at location " + address_to_string(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  SecondOrderODE<T>::~SecondOrderODE()
  {
    this->m_log->TRACE(NAME() + "SecondOrderODE '" + this->m_name
								+ "' destroyed at location " + address_to_string(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  SecondOrderODE<T>::SecondOrderODE(std::shared_ptr<Log> t_log)
  : DiffEQ<T>(t_log)
  {
    this->m_log->TRACE(NAME() + "SecondOrderODE '" + this->m_name
                       + "' created at location " + address_to_string(*this));
    this->m_log->INFO(NAME() + "Logger passed to SecondOrderODE '" + this->m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  SecondOrderODE<T>::SecondOrderODE(std::shared_ptr<Grid<T>> t_Grid,
                                  std::shared_ptr<Log> t_log)
  : DiffEQ<T>(t_Grid,t_log)
  {
    this->m_log->TRACE(NAME() + "SecondOrderODE '" + this->m_name
                       + "' created at location " + address_to_string(*this));
    this->m_log->INFO(NAME() + "Logger passed to SecondOrderODE '" + this->m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  T SecondOrderODE<T>::evaluate(const std::vector<T>& t_point, double t_time)
  {
    this->m_log->WARN(NAME() + "Called DEFAULT evaluate function for Second Order ODE");
    return 0;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<T>
  SecondOrderODE<T>::evaluate(const std::vector<std::vector<T>>& t_points,
                             double t_time)
  {
    this->m_log->WARN(NAME() + "Called DEFAULT evaluate function for Second Order ODE");
    return std::vector<T>();
  }
  //----------------------------------------------------------------------------
  //  Base class PDE
  //----------------------------------------------------------------------------
  template<typename T>
  PDE<T>::PDE()
  : DiffEQ<T>()
  {
    this->m_log = std::make_shared<Log>();
		this->m_log->init(NAME(), ".logs/diffeq_" + this->m_name + ".txt");
		this->m_log->TRACE(NAME() + "PDE '" + this->m_name
                       + "' created at location " + address_to_string(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  PDE<T>::~PDE()
  {
    this->m_log->TRACE(NAME() + "PDE '" + this->m_name
								+ "' destroyed at location " + address_to_string(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  PDE<T>::PDE(std::shared_ptr<Log> t_log)
  : DiffEQ<T>(t_log)
  {
    this->m_log->TRACE(NAME() + "PDE '" + this->m_name
                       + "' created at location " + address_to_string(*this));
    this->m_log->INFO(NAME() + "Logger passed to PDE '" + this->m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  PDE<T>::PDE(std::shared_ptr<Grid<T>> t_Grid,
                                  std::shared_ptr<Log> t_log)
  : DiffEQ<T>(t_Grid,t_log)
  {
    this->m_log->TRACE(NAME() + "PDE '" + this->m_name
                       + "' created at location " + address_to_string(*this));
    this->m_log->INFO(NAME() + "Logger passed to PDE '" + this->m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  T PDE<T>::evaluate(const std::vector<T>& t_point, double t_time)
  {
    this->m_log->WARN(NAME() + "Called DEFAULT evaluate function for PDE");
    return 0;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<T>
  PDE<T>::evaluate(const std::vector<std::vector<T>>& t_points,
                             double t_time)
  {
    this->m_log->WARN(NAME() + "Called DEFAULT evaluate function for PDE");
    return std::vector<T>();
  }
  //----------------------------------------------------------------------------


}
