//------------------------------------------------------------------------------
//  grid.cpp
//  The Entropic Trajectories Framework
//  -----------------------------------
//  Copyright (C) [2020] by [N. Carrara]
//  [ncarrara@albany.edu]
//
//  Permission to use, copy, modify, and/or distribute this software for any
//  purpose with or without fee t_is hereby granted.
//
//  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
//  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
//  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
//  SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
//  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
//  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR
//  IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
//------------------------------------------------------------------------------
#include "grid.h"

namespace ET
{
  template<typename T>
  Grid<T>::Grid()
  : m_name("default"), m_dim(0), m_N(0)
  {
    m_log = std::make_shared<Log>();
		m_log->init("ET:Grid:" + m_name, ".logs/grid_" + m_name + ".txt");
		m_log->TRACE("Grid '" + m_name + "' created at location "
		            + getMem(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Grid<T>::~Grid()
  {
    m_log->TRACE("Grid '" + m_name
								+ "' destroyed at location " + getMem(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Grid<T>::Grid(std::shared_ptr<Log> t_log)
  : m_name("default"), m_dim(0), m_N(0), m_log(t_log)
  {
    m_log->TRACE("Grid '" + m_name + "' created at location "
                + getMem(*this));
    m_log->INFO("Logger passed to Grid '" + m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Grid<T>::Grid(std::string t_name)
  : m_name(t_name), m_dim(0), m_N(0)
  {
    m_log = std::make_shared<Log>();
		m_log->init("ET:Grid:" + m_name, ".logs/grid_" + m_name + ".txt");
		m_log->TRACE("Grid '" + m_name + "' created at location "
		            + getMem(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Grid<T>::Grid(std::string t_name, std::shared_ptr<Log> t_log)
  : m_name(t_name), m_dim(0), m_N(0), m_log(t_log)
  {
    m_log->TRACE("Grid '" + m_name + "' created at location "
                + getMem(*this));
    m_log->INFO("Logger passed to Grid '" + m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Grid<T>::Grid(size_t t_dim)
  : m_name("default"), m_dim(t_dim), m_N(0)
  {
    m_log = std::make_shared<Log>();
		m_log->init("ET:Grid:" + m_name, ".logs/grid_" + m_name + ".txt");
		m_log->TRACE("Grid '" + m_name + "' created at location "
		            + getMem(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Grid<T>::Grid(size_t t_dim, std::shared_ptr<Log> t_log)
  : m_name("default"), m_dim(t_dim), m_N(0), m_log(t_log)
  {
    m_log->TRACE("Grid '" + m_name + "' created at location "
                + getMem(*this));
    m_log->INFO("Logger passed to Grid '" + m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Grid<T>::Grid(std::string t_name, size_t t_dim)
  : m_name(t_name), m_dim(t_dim), m_N(0)
  {
    m_log = std::make_shared<Log>();
		m_log->init("ET:Grid:" + m_name, ".logs/grid_" + m_name + ".txt");
		m_log->TRACE("Grid '" + m_name + "' created at location "
		            + getMem(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Grid<T>::Grid(std::string t_name, size_t t_dim, std::shared_ptr<Log> t_log)
  : m_name(t_name), m_dim(t_dim), m_N(0), m_log(t_log)
  {
    m_log->TRACE("Grid '" + m_name + "' created at location "
                + getMem(*this));
    m_log->INFO("Logger passed to Grid '" + m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Grid<T>::Grid(size_t t_dim, size_t t_N)
  : m_name("default"), m_dim(t_dim), m_N(t_N)
  {
    m_log = std::make_shared<Log>();
		m_log->init("ET:Grid:" + m_name, ".logs/grid_" + m_name + ".txt");
		m_log->TRACE("Grid '" + m_name + "' created at location "
		            + getMem(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Grid<T>::Grid(size_t t_dim, size_t t_N, std::shared_ptr<Log> t_log)
  : m_name("default"), m_dim(t_dim), m_N(t_N), m_log(t_log)
  {
    m_log->TRACE("Grid '" + m_name + "' created at location "
                + getMem(*this));
    m_log->INFO("Logger passed to Grid '" + m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Grid<T>::Grid(std::string t_name, size_t t_dim, size_t t_N)
  : m_name(t_name), m_dim(t_dim), m_N(t_N)
  {
    m_log = std::make_shared<Log>();
		m_log->init("ET:Grid:" + m_name, ".logs/grid_" + m_name + ".txt");
		m_log->TRACE("Grid '" + m_name + "' created at location "
		            + getMem(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Grid<T>::Grid(std::string t_name, size_t t_dim, size_t t_N,
                std::shared_ptr<Log> t_log)
  : m_name(t_name), m_dim(t_dim), m_N(t_N), m_log(t_log)
  {
    m_log->TRACE("Grid '" + m_name + "' created at location "
                + getMem(*this));
    m_log->INFO("Logger passed to Grid '" + m_name + "'");
  }
  //----------------------------------------------------------------------------
  //  Getters and Setters
  //----------------------------------------------------------------------------
  template<typename T>
  std::string Grid<T>::getName() const
  {
    return m_name;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  size_t Grid<T>::getDim() const
  {
    return m_dim;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  size_t Grid<T>::getN() const
  {
    return m_N;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::shared_ptr<Log> Grid<T>::getLog() const
  {
    return m_log;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Grid<T>::setName(const std::string t_name)
  {
    m_name = t_name;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Grid<T>::setDim(const size_t t_dim)
  {
    m_dim = t_dim;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Grid<T>::setN(const size_t t_N)
  {
    m_N = t_N;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Grid<T>::setLog(std::shared_ptr<Log> t_log)
  {
    m_log = t_log;
  }
  //----------------------------------------------------------------------------


}
