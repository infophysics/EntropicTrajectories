//------------------------------------------------------------------------------
//  ugrid.cpp
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
#include "ugrid.h"

namespace ET
{
  template<typename T>
  UGrid<T>::UGrid()
  : Grid<T>()
  {
		this->m_log->init(NAME(), ".logs/grid_" + this->m_name + ".txt");
		this->m_log->TRACE(NAME() + "UGrid '" + this->m_name + "' created at location "
		            + address_to_string(*this));
    m_kdtree = std::make_shared<KDTree<T>>();
  }
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::~UGrid()
  {
    this->m_log->TRACE(NAME() + "UGrid '" + this->m_name
								+ "' destroyed at location " + address_to_string(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(std::shared_ptr<Log> t_log)
  : Grid<T>(t_log)
  {
    this->m_log->TRACE(NAME() + "UGrid '" + this->m_name + "' created at location "
                + address_to_string(*this));
    this->m_log->INFO(NAME() + "Logger passed to UGrid '" + this->m_name + "'");
    m_kdtree = std::make_shared<KDTree<T>>();
  }
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(std::string t_name)
  : Grid<T>(t_name)
  {
    this->m_log = std::make_shared<Log>();
		this->m_log->init(NAME(), ".logs/grid_" + this->m_name + ".txt");
		this->m_log->TRACE(NAME() + "UGrid '" + this->m_name + "' created at location "
		            + address_to_string(*this));
    m_kdtree = std::make_shared<KDTree<T>>();
  }
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(std::string t_name, std::shared_ptr<Log> t_log)
  : Grid<T>(t_name, t_log)
  {
    this->m_log->TRACE(NAME() + "UGrid '" + this->m_name + "' created at location "
                + address_to_string(*this));
    this->m_log->INFO(NAME() + "Logger passed to UGrid '" + this->m_name + "'");
    m_kdtree = std::make_shared<KDTree<T>>();
  }
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(size_t t_dim)
  : Grid<T>(t_dim)
  {
    this->m_log = std::make_shared<Log>();
		this->m_log->init(NAME(), ".logs/grid_" + this->m_name + ".txt");
		this->m_log->TRACE(NAME() + "UGrid '" + this->m_name + "' created at location "
		            + address_to_string(*this));
    m_kdtree = std::make_shared<KDTree<T>>();
  }
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(size_t t_dim, std::shared_ptr<Log> t_log)
  : Grid<T>(t_dim, t_log)
  {
    this->m_log->TRACE(NAME() + "UGrid '" + this->m_name + "' created at location "
                + address_to_string(*this));
    this->m_log->INFO(NAME() + "Logger passed to UGrid '" + this->m_name + "'");
    m_kdtree = std::make_shared<KDTree<T>>();
  }
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(std::string t_name, size_t t_dim)
  : Grid<T>(t_name, t_dim)
  {
    this->m_log = std::make_shared<Log>();
		this->m_log->init(NAME(), ".logs/grid_" + this->m_name + ".txt");
		this->m_log->TRACE(NAME() + "UGrid '" + this->m_name + "' created at location "
		            + address_to_string(*this));
    m_kdtree = std::make_shared<KDTree<T>>();
  }
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(std::string t_name, size_t t_dim, std::shared_ptr<Log> t_log)
  : Grid<T>(t_name, t_dim, t_log)
  {
    this->m_log->TRACE(NAME() + "UGrid '" + this->m_name + "' created at location "
                + address_to_string(*this));
    this->m_log->INFO(NAME() + "Logger passed to UGrid '" + this->m_name + "'");
    m_kdtree = std::make_shared<KDTree<T>>();
  }
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(size_t t_dim, size_t t_N)
  : Grid<T>(t_dim, t_N)
  {
    this->m_log = std::make_shared<Log>();
		this->m_log->init(NAME(), ".logs/grid_" + this->m_name + ".txt");
		this->m_log->TRACE(NAME() + "UGrid '" + this->m_name + "' created at location "
		            + address_to_string(*this));
    m_kdtree = std::make_shared<KDTree<T>>();
  }
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(size_t t_dim, size_t t_N, std::shared_ptr<Log> t_log)
  : Grid<T>(t_dim, t_N, t_log)
  {
    this->m_log->TRACE(NAME() + "UGrid '" + this->m_name + "' created at location "
                + address_to_string(*this));
    this->m_log->INFO(NAME() + "Logger passed to UGrid '" + this->m_name + "'");
    m_kdtree = std::make_shared<KDTree<T>>();
  }
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(std::string t_name, size_t t_dim, size_t t_N)
  : Grid<T>(t_name, t_dim, t_N)
  {
    this->m_log = std::make_shared<Log>();
		this->m_log->init(NAME(), ".logs/grid_" + this->m_name + ".txt");
		this->m_log->TRACE(NAME() + "UGrid '" + this->m_name + "' created at location "
		            + address_to_string(*this));
    m_kdtree = std::make_shared<KDTree<T>>();
  }
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(std::string t_name, size_t t_dim, size_t t_N,
                std::shared_ptr<Log> t_log)
  : Grid<T>(t_name, t_dim, t_N, t_log)
  {
    this->m_log->TRACE(NAME() + "UGrid '" + this->m_name + "' created at location "
                + address_to_string(*this));
    this->m_log->INFO(NAME() + "Logger passed to UGrid '" + this->m_name + "'");
    m_kdtree = std::make_shared<KDTree<T>>();
  }
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(std::vector<std::vector<T>> t_grid, bool move_grid)
  : Grid<T>(t_grid,move_grid)
  {
    this->m_log = std::make_shared<Log>();
		this->m_log->init(NAME(), ".logs/grid_" + this->m_name + ".txt");
		this->m_log->TRACE(NAME() + "UGrid '" + this->m_name + "' created at location "
		            + address_to_string(*this));
    m_kdtree = std::make_shared<KDTree<T>>(std::make_shared<std::vector<std::vector<T>>>(this->m_grid));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(std::vector<std::vector<T>> t_grid, std::shared_ptr<Log> t_log,
                  bool move_grid)
  : Grid<T>(t_grid, t_log, move_grid)
  {
    this->m_log->TRACE(NAME() + "UGrid '" + this->m_name + "' created at location "
                + address_to_string(*this));
    this->m_log->INFO(NAME() + "Logger passed to UGrid '" + this->m_name + "'");
    m_kdtree = std::make_shared<KDTree<T>>(std::make_shared<std::vector<std::vector<T>>>(this->m_grid));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(std::string t_name, std::vector<std::vector<T>> t_grid,
                  bool move_grid)
  : Grid<T>(t_name, t_grid, move_grid)
  {
    this->m_log = std::make_shared<Log>();
		this->m_log->init(NAME(), ".logs/grid_" + this->m_name + ".txt");
		this->m_log->TRACE(NAME() + "UGrid '" + this->m_name + "' created at location "
		            + address_to_string(*this));
    m_kdtree = std::make_shared<KDTree<T>>(std::make_shared<std::vector<std::vector<T>>>(this->m_grid));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(std::string t_name, std::vector<std::vector<T>> t_grid,
                std::shared_ptr<Log> t_log, bool move_grid)
  : Grid<T>(t_name, t_grid, t_log, move_grid)
  {
    this->m_log->TRACE(NAME() + "UGrid '" + this->m_name + "' created at location "
                + address_to_string(*this));
    this->m_log->INFO(NAME() + "Logger passed to UGrid '" + this->m_name + "'");
    m_kdtree = std::make_shared<KDTree<T>>(std::make_shared<std::vector<std::vector<T>>>(this->m_grid));
  }
  //----------------------------------------------------------------------------
  //  Getters and Setters
  //----------------------------------------------------------------------------
  template<typename T>
  enum GridType UGrid<T>::getType() const
  {
    return m_type;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::shared_ptr<KDTree<T>> UGrid<T>::getKDTree() const
  {
    return m_kdtree;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void UGrid<T>::setKDTree(std::shared_ptr<KDTree<T>> t_kdtree)
  {
    m_kdtree = t_kdtree;
  }

	//----------------------------------------------------------------------------
	//	Various functions
	//----------------------------------------------------------------------------
	template<typename T>
	bool UGrid<T>::checkConsistency()
	{
		bool consistency = true;
		if(this->m_grid.size() != this->m_N)
		{
			consistency = false;
			this->m_log->WARN(NAME() + "Attribute this->m_N with value "
			           + std::to_string(this->m_N) + " does not match array this->m_grid"
								 + " with size " + std::to_string(this->m_grid.size()));
		}
		if(this->m_grid[0].size() != this->m_dim)
		{
			consistency = false;
			this->m_log->WARN(NAME() + "Attribute this->m_dim with value "
			           + std::to_string(this->m_dim) + " does not match array this->m_grid"
								 + " with dimension " + std::to_string(this->m_grid[0].size()));
		}
		return consistency;
	}
	//----------------------------------------------------------------------------
  //  KDTree methods
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<std::vector<size_t>> UGrid<T>::getCurrentNeighborIndices()
  {
    return m_kdtree->getCurrentNeighborIndices();
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<std::vector<double>> UGrid<T>::getCurrentNeighborDistances()
  {
    return m_kdtree->getCurrentNeighborDistances();
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<size_t> UGrid<T>::getCurrentNeighborIndices(size_t t_i)
  {
    return m_kdtree->getCurrentNeighborIndices(t_i);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<double> UGrid<T>::getCurrentNeighborDistances(size_t t_i)
  {
    return m_kdtree->getCurrentNeighborDistances(t_i);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  size_t UGrid<T>::getCurrentGlobalK() const
  {
    return m_kdtree->getCurrentGlobalK();
  }
  //----------------------------------------------------------------------------
  template<typename T>
  double UGrid<T>::getCurrentGlobalRadius() const
  {
    return m_kdtree->getCurrentGlobalRadius();
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void UGrid<T>::setCurrentGlobalK(const size_t t_k)
  {
    return m_kdtree->setCurrentGlobalK(t_k);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void UGrid<T>::setCurrentGlobalRadius(const double t_radius)
  {
    return m_kdtree->setCurrentGlobalRadius(t_radius);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void UGrid<T>::setupTree()
  {
    return m_kdtree->setupTree();
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void UGrid<T>::queryNeighbors(size_t t_k)
  {
    return m_kdtree->queryNeighbors(t_k);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void UGrid<T>::queryNeighbors(double t_radius)
  {
    return m_kdtree->queryNeighbors(t_radius);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<size_t>
  UGrid<T>::queryNeighbors(const std::vector<T>& t_point, size_t t_k)
  {
    return m_kdtree->queryNeighbors(t_point,t_k);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<double>
  UGrid<T>::queryDistances(const std::vector<T>& t_point, size_t t_k)
  {
    return m_kdtree->queryDistances(t_point,t_k);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::tuple<std::vector<size_t>,std::vector<double>>
  UGrid<T>::query(const std::vector<T>& t_point, size_t t_k)
  {
    return m_kdtree->query(t_point,t_k);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<std::vector<size_t>>
  UGrid<T>::queryNeighbors(const std::vector<std::vector<T>>& t_points,
                            size_t t_k)
  {
    return m_kdtree->queryNeighbors(t_points,t_k);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<std::vector<double>>
  UGrid<T>::queryDistances(const std::vector<std::vector<T>>& t_points,
                            size_t t_k)
  {
    return m_kdtree->queryDistances(t_points,t_k);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::tuple<std::vector<std::vector<size_t>>,std::vector<std::vector<double>>>
  UGrid<T>::query(const std::vector<std::vector<T>>& t_points,
                            size_t t_k)
  {
    return m_kdtree->query(t_points,t_k);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<size_t>
  UGrid<T>::queryNeighbors(const std::vector<T>& t_point, double t_radius)
  {
    return m_kdtree->queryNeighbors(t_point,t_radius);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<double>
  UGrid<T>::queryDistances(const std::vector<T>& t_point, double t_radius)
  {
    return m_kdtree->queryDistances(t_point,t_radius);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::tuple<std::vector<size_t>,std::vector<double>>
  UGrid<T>::query(const std::vector<T>& t_point, double t_radius)
  {
    return m_kdtree->query(t_point,t_radius);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<std::vector<size_t>>
  UGrid<T>::queryNeighbors(const std::vector<std::vector<T>>& t_points,
                            double t_radius)
  {
    return m_kdtree->queryNeighbors(t_points,t_radius);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<std::vector<double>>
  UGrid<T>::queryDistances(const std::vector<std::vector<T>>& t_points,
                            double t_radius)
  {
    return m_kdtree->queryDistances(t_points,t_radius);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::tuple<std::vector<std::vector<size_t>>,std::vector<std::vector<double>>>
  UGrid<T>::query(const std::vector<std::vector<T>>& t_points,
                   double t_radius)
  {
    return m_kdtree->query(t_points,t_radius);
  }
  //----------------------------------------------------------------------------
}
