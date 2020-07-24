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
		m_log->init(NAME(), ".logs/grid_" + m_name + ".txt");
		m_log->TRACE(NAME() + "Grid '" + m_name + "' created at location "
		            + address_to_string(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Grid<T>::~Grid()
  {
    m_log->TRACE(NAME() + "Grid '" + m_name
								+ "' destroyed at location " + address_to_string(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Grid<T>::Grid(std::shared_ptr<Log> t_log)
  : m_name("default"), m_dim(0), m_N(0), m_log(t_log)
  {
    m_log->TRACE(NAME() + "Grid '" + m_name + "' created at location "
                + address_to_string(*this));
    m_log->INFO(NAME() + "Logger passed to Grid '" + m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Grid<T>::Grid(std::string t_name)
  : m_name(t_name), m_dim(0), m_N(0)
  {
    m_log = std::make_shared<Log>();
		m_log->init(NAME(), ".logs/grid_" + m_name + ".txt");
		m_log->TRACE(NAME() + "Grid '" + m_name + "' created at location "
		            + address_to_string(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Grid<T>::Grid(std::string t_name, std::shared_ptr<Log> t_log)
  : m_name(t_name), m_dim(0), m_N(0), m_log(t_log)
  {
    m_log->TRACE(NAME() + "Grid '" + m_name + "' created at location "
                + address_to_string(*this));
    m_log->INFO(NAME() + "Logger passed to Grid '" + m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Grid<T>::Grid(size_t t_dim)
  : m_name("default"), m_dim(t_dim), m_N(0)
  {
    m_log = std::make_shared<Log>();
		m_log->init(NAME(), ".logs/grid_" + m_name + ".txt");
		m_log->TRACE(NAME() + "Grid '" + m_name + "' created at location "
		            + address_to_string(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Grid<T>::Grid(size_t t_dim, std::shared_ptr<Log> t_log)
  : m_name("default"), m_dim(t_dim), m_N(0), m_log(t_log)
  {
    m_log->TRACE(NAME() + "Grid '" + m_name + "' created at location "
                + address_to_string(*this));
    m_log->INFO(NAME() + "Logger passed to Grid '" + m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Grid<T>::Grid(std::string t_name, size_t t_dim)
  : m_name(t_name), m_dim(t_dim), m_N(0)
  {
    m_log = std::make_shared<Log>();
		m_log->init(NAME(), ".logs/grid_" + m_name + ".txt");
		m_log->TRACE(NAME() + "Grid '" + m_name + "' created at location "
		            + address_to_string(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Grid<T>::Grid(std::string t_name, size_t t_dim, std::shared_ptr<Log> t_log)
  : m_name(t_name), m_dim(t_dim), m_N(0), m_log(t_log)
  {
    m_log->TRACE(NAME() + "Grid '" + m_name + "' created at location "
                + address_to_string(*this));
    m_log->INFO(NAME() + "Logger passed to Grid '" + m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Grid<T>::Grid(size_t t_dim, size_t t_N)
  : m_name("default"), m_dim(t_dim), m_N(t_N)
  {
    m_log = std::make_shared<Log>();
		m_log->init(NAME(), ".logs/grid_" + m_name + ".txt");
		m_log->TRACE(NAME() + "Grid '" + m_name + "' created at location "
		            + address_to_string(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Grid<T>::Grid(size_t t_dim, size_t t_N, std::shared_ptr<Log> t_log)
  : m_name("default"), m_dim(t_dim), m_N(t_N), m_log(t_log)
  {
    m_log->TRACE(NAME() + "Grid '" + m_name + "' created at location "
                + address_to_string(*this));m_log = std::make_shared<Log>();
		m_log->init("ET:KDTree:default", ".logs/KDTree_default.txt");
		m_log->TRACE(NAME() + "KDTree 'default' created at location "
		            + address_to_string(*this));
    m_log->INFO(NAME() + "Logger passed to Grid '" + m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Grid<T>::Grid(std::string t_name, size_t t_dim, size_t t_N)
  : m_name(t_name), m_dim(t_dim), m_N(t_N)
  {
    m_log = std::make_shared<Log>();
		m_log->init(NAME(), ".logs/grid_" + m_name + ".txt");
		m_log->TRACE(NAME() + "Grid '" + m_name + "' created at location "
		            + address_to_string(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Grid<T>::Grid(std::string t_name, size_t t_dim, size_t t_N,
                std::shared_ptr<Log> t_log)
  : m_name(t_name), m_dim(t_dim), m_N(t_N), m_log(t_log)
  {
    m_log->TRACE(NAME() + "Grid '" + m_name + "' created at location "
                + address_to_string(*this));
    m_log->INFO(NAME() + "Logger passed to Grid '" + m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Grid<T>::Grid(std::vector<std::vector<T>> t_grid)
  : m_name("default"), m_dim(t_grid[0].size()), m_N(t_grid.size()),
    m_grid(t_grid)
  {
    m_log = std::make_shared<Log>();
		m_log->init(NAME(), ".logs/grid_" + m_name + ".txt");
		m_log->TRACE(NAME() + "Grid '" + m_name + "' created at location "
		            + address_to_string(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Grid<T>::Grid(std::vector<std::vector<T>> t_grid, std::shared_ptr<Log> t_log)
  : m_name("default"), m_dim(t_grid[0].size()), m_N(t_grid.size()),
    m_grid(t_grid)
  {
    m_log->TRACE(NAME() + "Grid '" + m_name + "' created at location "
                + address_to_string(*this));
    m_log->INFO(NAME() + "Logger passed to Grid '" + m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Grid<T>::Grid(std::string t_name, std::vector<std::vector<T>> t_grid)
  : m_name(t_name), m_dim(t_grid[0].size()), m_N(t_grid.size()),
    m_grid(t_grid)
  {
    m_log = std::make_shared<Log>();
		m_log->init(NAME(), ".logs/grid_" + m_name + ".txt");
		m_log->TRACE(NAME() + "Grid '" + m_name + "' created at location "
		            + address_to_string(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Grid<T>::Grid(std::string t_name, std::vector<std::vector<T>> t_grid,
                std::shared_ptr<Log> t_log)
  : m_name(t_name), m_dim(t_grid[0].size()), m_N(t_grid.size()),
    m_grid(t_grid)
  {
    m_log->TRACE(NAME() + "Grid '" + m_name + "' created at location "
                + address_to_string(*this));
    m_log->INFO(NAME() + "Logger passed to Grid '" + m_name + "'");
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
  std::vector<std::vector<T>> Grid<T>::getGrid() const
  {
    return m_grid;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<std::string> Grid<T>::getCoords() const
  {
    return m_coords;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::shared_ptr<Log> Grid<T>::getLog() const
  {
    return m_log;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  enum GridType Grid<T>::getType() const
  {
    return m_type;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<T> Grid<T>::getPoint(const size_t t_i) const
  {
    //  Check that indicies are within the bounds
    if (t_i >= m_N) {
      m_log->ERROR(NAME() + "Attempted to access the (" + std::to_string(t_i)
                   + ")th point of Grid with "
                   + "size (" + std::to_string(m_dim) + " x "
                   + std::to_string(m_N) + ").");
      m_log->INFO(NAME() + "Returning the first element.");
      return m_grid[0];
    }
    //  otherwise return the correct element
    return m_grid[t_i];
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Grid<T>::setName(const std::string t_name)
  {
    m_log->TRACE(NAME() + "Changed name from '" + m_name
                 + "' to '" + t_name + "'");
    m_name = t_name;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Grid<T>::setDim(const size_t t_dim)
  {
    m_log->TRACE(NAME() + "Changed dim from '" + std::to_string(m_dim)
                 + "' to '" + std::to_string(t_dim) + "'");
    m_dim = t_dim;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Grid<T>::setN(const size_t t_N)
  {
    m_log->TRACE(NAME() + "Changed N from '" + std::to_string(m_N)
                 + "' to '" + std::to_string(t_N) + "'");
    m_N = t_N;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Grid<T>::setGrid(const std::vector<std::vector<T>> t_grid)
  {
    m_grid = t_grid;
    m_log->TRACE(NAME() + "Set grid to array of size (" + std::to_string(m_dim)
                + " x " + std::to_string(m_N) + ")");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Grid<T>::setCoords(const std::vector<std::string> t_coords)
  {
    m_coords = t_coords;
    m_log->TRACE(NAME() + "Changed coords attribute");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Grid<T>::setLog(std::shared_ptr<Log> t_log)
  {
    m_log->TRACE(NAME() + "Changed log from '" + address_to_string(m_log)
                 + "' to '" + address_to_string(t_log) + "'");
    m_log = t_log;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Grid<T>& Grid<T>::operator=(const Grid<T>& t_Grid)
  {
    if (&t_Grid == this) {
      return *this;
    }
    m_name = t_Grid.getName();
    m_dim = t_Grid.getDim();
    m_N = t_Grid.getN();
    m_coords = t_Grid.getCoords();
    m_log = t_Grid.getLog();
    m_grid.resize(m_dim);
    for (auto i = 0; i < m_N; i++) {
        m_grid[i] = t_Grid(i);
    }
    return *this;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  bool Grid<T>::operator==(const Grid<T>& t_Grid) const
  {
    return m_grid == t_Grid.getGrid();
  }
  //----------------------------------------------------------------------------
  template<typename T>
  bool Grid<T>::operator!=(const Grid<T>& t_Grid) const
  {
    return m_grid != t_Grid.getGrid();
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Grid<T> Grid<T>::operator-() const
  {
    auto grid = m_grid;
    for (auto i = 0; i < m_N; i++) {
      for (auto j = 0; j < m_dim; j++) {
        grid[i][j] *= -1;
      }
    }
    auto name = "-" + m_name;
    return Grid<T>(name,grid,m_log);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Grid<T> Grid<T>::operator+(const Grid<T>& t_Grid) const
  {
    //  Check that grids are compatible
    if (m_N != t_Grid.getN() || m_dim != t_Grid.getDim()) {
      m_log->ERROR(NAME() + "Attempted to add two grids of size ("
                   + std::to_string(m_dim) + " x " + std::to_string(m_N)
                   + ") and (" + std::to_string(t_Grid.getDim()) + " x "
                   + std::to_string(t_Grid.getN()) + ").");
      m_log->INFO(NAME() + "Returning lvalue Grid.");
      return *this;
    }
    //  If so then add the grids element by element
    auto grid = m_grid;
    for (auto i = 0; i < m_N; i++) {
      for (auto j = 0; j < m_dim; j++) {
        grid[i][j] += t_Grid(i,j);
      }
    }
    auto name = m_name + "+" + t_Grid.getName();
    return Grid<T>(name,grid,m_log);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Grid<T>& Grid<T>::operator+=(const Grid<T>& t_Grid)
  {
    //  Check that grids are compatible
    if (m_N != t_Grid.getN() || m_dim != t_Grid.getDim()) {
      m_log->ERROR(NAME() + "Attempted to add two grids of size ("
                   + std::to_string(m_dim) + " x " + std::to_string(m_N)
                   + ") and (" + std::to_string(t_Grid.getDim()) + " x "
                   + std::to_string(t_Grid.getN()) + ").");
      m_log->INFO(NAME() + "Returning lvalue Grid.");
      return *this;
    }
    //  If so then add the grids element by element
    for (auto i = 0; i < m_N; i++) {
      for (auto j = 0; j < m_dim; j++) {
        m_grid[i][j] += t_Grid(i,j);
      }
    }
    m_name += "+" + t_Grid.getName();
    return *this;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Grid<T> Grid<T>::operator-(const Grid<T>& t_Grid) const
  {
    //  Check that grids are compatible
    if (m_N != t_Grid.getN() || m_dim != t_Grid.getDim()) {
      m_log->ERROR(NAME() + "Attempted to subtract two grids of size ("
                   + std::to_string(m_dim) + " x " + std::to_string(m_N)
                   + ") and (" + std::to_string(t_Grid.getDim()) + " x "
                   + std::to_string(t_Grid.getN()) + ").");
      m_log->INFO(NAME() + "Returning lvalue Grid.");
      return *this;
    }
    //  If so then add the grids element by element
    auto grid = m_grid;
    for (auto i = 0; i < m_N; i++) {
      for (auto j = 0; j < m_dim; j++) {
        grid[i][j] -= t_Grid(i,j);
      }
    }
    auto name = m_name + "-" + t_Grid.getName();
    return Grid<T>(name,grid,m_log);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Grid<T>& Grid<T>::operator-=(const Grid<T>& t_Grid)
  {
    //  Check that grids are compatible
    if (m_N != t_Grid.getN() || m_dim != t_Grid.getDim()) {
      m_log->ERROR(NAME() + "Attempted to subtract two grids of size ("
                   + std::to_string(m_dim) + " x " + std::to_string(m_N)
                   + ") and (" + std::to_string(t_Grid.getDim()) + " x "
                   + std::to_string(t_Grid.getN()) + ").");
      m_log->INFO(NAME() + "Returning lvalue Grid.");
      return *this;
    }
    //  If so then add the grids element by element
    for (auto i = 0; i < m_N; i++) {
      for (auto j = 0; j < m_dim; j++) {
        m_grid[i][j] -= t_Grid(i,j);
      }
    }
    m_name += "-" + t_Grid.getName();
    return *this;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  T& Grid<T>::operator()(const size_t& t_i, const size_t& t_j)
  {
    //  Check that indicies are within the bounds
    if (t_i >= m_N || t_j >= m_dim) {
      m_log->ERROR(NAME() + "Attempted to access the (" + std::to_string(t_i)
                   + "," + std::to_string(t_j) + ")th element of Grid with "
                   + "size (" + std::to_string(m_dim) + " x "
                   + std::to_string(m_N) + ").");
      m_log->INFO(NAME() + "Returning the first element.");
      return m_grid[0][0];
    }
    //  otherwise return the correct element
    return m_grid[t_i][t_j];
  }
  //----------------------------------------------------------------------------
  template<typename T>
  const T& Grid<T>::operator()(const size_t& t_i, const size_t& t_j) const
  {
    //  Check that indicies are within the bounds
    if (t_i >= m_N || t_j >= m_dim) {
      m_log->ERROR(NAME() + "Attempted to access the (" + std::to_string(t_i)
                   + "," + std::to_string(t_j) + ")th element of Grid with "
                   + "size (" + std::to_string(m_dim) + " x "
                   + std::to_string(m_N) + ").");
      m_log->INFO(NAME() + "Returning the first element.");
      return m_grid[0][0];
    }
    //  otherwise return the correct element
    return m_grid[t_i][t_j];
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<T>& Grid<T>::operator()(const size_t& t_i)
  {
    //  Check that indicies are within the bounds
    if (t_i >= m_N) {
      m_log->ERROR(NAME() + "Attempted to access the (" + std::to_string(t_i)
                   + ")th point of Grid with "
                   + "size (" + std::to_string(m_dim) + " x "
                   + std::to_string(m_N) + ").");
      m_log->INFO(NAME() + "Returning the first element.");
      return m_grid[0];
    }
    //  otherwise return the correct element
    return m_grid[t_i];
  }
  //----------------------------------------------------------------------------
  template<typename T>
  const std::vector<T>& Grid<T>::operator()(const size_t& t_i) const
  {
    //  Check that indicies are within the bounds
    if (t_i >= m_N) {
      m_log->ERROR(NAME() + "Attempted to access the (" + std::to_string(t_i)
                   + ")th point of Grid with "
                   + "size (" + std::to_string(m_dim) + " x "
                   + std::to_string(m_N) + ").");
      m_log->INFO(NAME() + "Returning the first element.");
      return m_grid[0];
    }
    //  otherwise return the correct element
    return m_grid[t_i];
  }
  //----------------------------------------------------------------------------
  //  Special functions
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<T> Grid<T>::proj(const size_t t_axis)
  {
    std::vector<T> projection(m_N,0);
    //  Check that the axis exists
    if (t_axis >= m_dim) {
      m_log->ERROR(NAME() + "Attempted to project along axis " + std::to_string(t_axis)
                   + " when Grid has dimension " + std::to_string(m_dim));
      m_log->INFO(NAME() + "Returning empty std::vector");
      return projection;
    }
    //  otherwise construct the projection
    for (auto i = 0; i < m_N; i++) {
      projection[i] = m_grid[i][t_axis];
    }
    return projection;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<std::vector<T>> Grid<T>::proj(const std::vector<size_t> t_axes)
  {
    //  Check that each the axis' exists
    std::vector<std::vector<T>> projection(m_N,std::vector<T>(t_axes.size(),0));
    for (auto i = 0; i < t_axes.size(); i++) {
      if (t_axes[i] >= m_dim) {
        m_log->ERROR(NAME() + "Attempted to project along axis "
                     + std::to_string(t_axes[i]) +" when Grid has dimension "
                     + std::to_string(m_dim));
        m_log->INFO(NAME() + "Returning empty std::vector");
        return projection;
      }
    }
    //  otherwise construct the projection
    for (auto i = 0; i < m_N; i++) {
      for (auto j = 0; j < t_axes.size(); j++) {
        projection[i][j] = m_grid[i][t_axes[j]];
      }
    }
    return projection;
  }

}
