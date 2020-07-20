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
    this->m_log = std::make_shared<Log>();
		this->m_log->init("ET:UGrid:" + this->m_name, ".logs/grid_" + this->m_name + ".txt");
		this->m_log->TRACE("UGrid '" + this->m_name + "' created at location "
		            + getMem(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::~UGrid()
  {
    this->m_log->TRACE("UGrid '" + this->m_name
								+ "' destroyed at location " + getMem(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(std::shared_ptr<Log> t_log)
  : Grid<T>(t_log)
  {
    this->m_log->TRACE("UGrid '" + this->m_name + "' created at location "
                + getMem(*this));
    this->m_log->INFO("Logger passed to UGrid '" + this->m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(std::string t_name)
  : Grid<T>(t_name)
  {
    this->m_log = std::make_shared<Log>();
		this->m_log->init("ET:UGrid:" + this->m_name, ".logs/grid_" + this->m_name + ".txt");
		this->m_log->TRACE("UGrid '" + this->m_name + "' created at location "
		            + getMem(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(std::string t_name, std::shared_ptr<Log> t_log)
  : Grid<T>(t_name, t_log)
  {
    this->m_log->TRACE("UGrid '" + this->m_name + "' created at location "
                + getMem(*this));
    this->m_log->INFO("Logger passed to UGrid '" + this->m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(size_t t_dim)
  : Grid<T>(t_dim)
  {
    this->m_log = std::make_shared<Log>();
		this->m_log->init("ET:UGrid:" + this->m_name, ".logs/grid_" + this->m_name + ".txt");
		this->m_log->TRACE("UGrid '" + this->m_name + "' created at location "
		            + getMem(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(size_t t_dim, std::shared_ptr<Log> t_log)
  : Grid<T>(t_dim, t_log)
  {
    this->m_log->TRACE("UGrid '" + this->m_name + "' created at location "
                + getMem(*this));
    this->m_log->INFO("Logger passed to UGrid '" + this->m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(std::string t_name, size_t t_dim)
  : Grid<T>(t_name, t_dim)
  {
    this->m_log = std::make_shared<Log>();
		this->m_log->init("ET:UGrid:" + this->m_name, ".logs/grid_" + this->m_name + ".txt");
		this->m_log->TRACE("UGrid '" + this->m_name + "' created at location "
		            + getMem(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(std::string t_name, size_t t_dim, std::shared_ptr<Log> t_log)
  : Grid<T>(t_name, t_dim, t_log)
  {
    this->m_log->TRACE("UGrid '" + this->m_name + "' created at location "
                + getMem(*this));
    this->m_log->INFO("Logger passed to UGrid '" + this->m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(size_t t_dim, size_t t_N)
  : Grid<T>(t_dim, t_N)
  {
    this->m_log = std::make_shared<Log>();
		this->m_log->init("ET:UGrid:" + this->m_name, ".logs/grid_" + this->m_name + ".txt");
		this->m_log->TRACE("UGrid '" + this->m_name + "' created at location "
		            + getMem(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(size_t t_dim, size_t t_N, std::shared_ptr<Log> t_log)
  : Grid<T>(t_dim, t_N, t_log)
  {
    this->m_log->TRACE("UGrid '" + this->m_name + "' created at location "
                + getMem(*this));
    this->m_log->INFO("Logger passed to UGrid '" + this->m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(std::string t_name, size_t t_dim, size_t t_N)
  : Grid<T>(t_name, t_dim, t_N)
  {
    this->m_log = std::make_shared<Log>();
		this->m_log->init("ET:UGrid:" + this->m_name, ".logs/grid_" + this->m_name + ".txt");
		this->m_log->TRACE("UGrid '" + this->m_name + "' created at location "
		            + getMem(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(std::string t_name, size_t t_dim, size_t t_N,
                std::shared_ptr<Log> t_log)
  : Grid<T>(t_name, t_dim, t_N, t_log)
  {
    this->m_log->TRACE("UGrid '" + this->m_name + "' created at location "
                + getMem(*this));
    this->m_log->INFO("Logger passed to UGrid '" + this->m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(std::vector<std::vector<T>> t_grid)
  : Grid<T>(t_grid)
  {
    this->m_log = std::make_shared<Log>();
		this->m_log->init("ET:UGrid:" + this->m_name, ".logs/grid_" + this->m_name + ".txt");
		this->m_log->TRACE("UGrid '" + this->m_name + "' created at location "
		            + getMem(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(std::vector<std::vector<T>> t_grid, std::shared_ptr<Log> t_log)
  : Grid<T>(t_grid, t_log)
  {
    this->m_log->TRACE("UGrid '" + this->m_name + "' created at location "
                + getMem(*this));
    this->m_log->INFO("Logger passed to UGrid '" + this->m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(std::string t_name, std::vector<std::vector<T>> t_grid)
  : Grid<T>(t_name, t_grid)
  {
    this->m_log = std::make_shared<Log>();
		this->m_log->init("ET:UGrid:" + this->m_name, ".logs/grid_" + this->m_name + ".txt");
		this->m_log->TRACE("UGrid '" + this->m_name + "' created at location "
		            + getMem(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(std::string t_name, std::vector<std::vector<T>> t_grid,
                std::shared_ptr<Log> t_log)
  : Grid<T>(t_name, t_grid, t_log)
  {
    this->m_log->TRACE("UGrid '" + this->m_name + "' created at location "
                + getMem(*this));
    this->m_log->INFO("Logger passed to UGrid '" + this->m_name + "'");
  }
  //----------------------------------------------------------------------------
  //  Getters and Setters
  //----------------------------------------------------------------------------
  template<typename T>
  KDTree<T> UGrid<T>::getKDTree() const
  {
    return m_kdtree;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void UGrid<T>::setKDTree(KDTree<T> t_kdtree)
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
			this->m_log->WARN("UGrid " + this->m_name + ": Attribute this->m_N with value "
			           + std::to_string(this->m_N) + " does not match array this->m_grid"
								 + " with size " + std::to_string(this->m_grid.size()));
		}
		if(this->m_grid[0].size() != this->m_dim)
		{
			consistency = false;
			this->m_log->WARN("UGrid " + this->m_name + ": Attribute this->m_dim with value "
			           + std::to_string(this->m_dim) + " does not match array this->m_grid"
								 + " with dimension " + std::to_string(this->m_grid[0].size()));
		}
		return consistency;
	}
	//----------------------------------------------------------------------------
}
