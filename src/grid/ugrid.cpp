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

	template<typename T>
	UGrid<T>::UGrid(std::vector<T> ugrid)
  : Grid<T>()
	{
		this->m_N = ugrid.size();
		m_ugrid.resize(this->m_N);
		this->m_dim = 1;
		for (uint32_t i = 0; i < this->m_N; i++)
		{
			std::vector<T> temp = {ugrid[i]};
			m_ugrid[i] = temp;
		}
		//##########################################################################
		this->m_log = std::make_shared<Log>();
		this->m_log->init("ET:UGrid:default", ".logs/ugrid_default.txt");
		this->m_log->TRACE("Unstructured Grid 'default' created at location "
		            + getMem(*this));
		//##########################################################################
		//  generate KDTree
  	setupTree();
		m_searchFlag = -1;
	}

	template<typename T>
	UGrid<T>::UGrid(std::vector<std::vector<T>> ugrid)
  : Grid<T>()
	{
		this->m_N = ugrid.size();
		this->m_dim = ugrid[0].size();
		m_ugrid = ugrid;
		//##########################################################################
		this->m_log = std::make_shared<Log>();
		this->m_log->init("ET:UGrid:default", ".logs/ugrid_default.txt");
		this->m_log->TRACE("Unstructured Grid 'default' created at location "
		            + getMem(*this));
	  //##########################################################################
		//  generate KDTree
  	setupTree();
		m_searchFlag = -1;
	}

	template<typename T>
	UGrid<T>::UGrid(std::vector<T> ugrid, std::shared_ptr<Log> log)
	: Grid<T>()
	{
		this->m_N = ugrid.size();
		m_ugrid.resize(this->m_N);
		this->m_dim = 1;
		for (uint32_t i = 0; i < this->m_N; i++)
		{
			std::vector<T> temp = {ugrid[i]};
			m_ugrid[i] = temp;
		}
		//##########################################################################
		this->m_log = log;
		this->m_log->TRACE("Unstructured Grid 'default' created at location "
								+ getMem(*this));
		this->m_log->INFO("Logger passed to Unstructured Grid 'default'");
		//##########################################################################
		//  generate KDTree
  	setupTree();
		m_searchFlag = -1;
	}

	template<typename T>
	UGrid<T>::UGrid(std::vector<std::vector<T>> ugrid, std::shared_ptr<Log> log)
	: Grid<T>()
	{
		this->m_N = ugrid.size();
		this->m_dim = ugrid[0].size();
		m_ugrid = ugrid;
		//##########################################################################
		this->m_log = log;
		this->m_log->TRACE("Unstructured Grid 'default' created at location "
								+ getMem(*this));
		this->m_log->INFO("Logger passed to Unstructured Grid 'default'");
		//##########################################################################
		//  generate KDTree
  	setupTree();
		m_searchFlag = -1;
	}
	//----------------------------------------------------------------------------


  template<typename T>
  std::vector<std::vector<T>> UGrid<T>::getUGrid() const
  {
    return m_ugrid;
  }
	template<typename T>
	std::vector<std::vector<size_t>> UGrid<T>::getNeighbors()
	{
		return m_kdt.getNeighbors();
	}
	template<typename T>
	std::vector<std::vector<double>> UGrid<T>::getDistances()
	{
		return m_kdt.getDistances();
	}
	template<typename T>
	std::vector<std::vector<size_t>> UGrid<T>::getNeighborsRadius()
	{
		return m_kdt.getNeighborsRadius();
	}
	template<typename T>
	std::vector<std::vector<double>> UGrid<T>::getDistancesRadius()
	{
		return m_kdt.getDistancesRadius();
	}
	template<typename T>
	std::vector<size_t> UGrid<T>::getNeighbors(uint64_t index)
	{
		if (index > this->m_N)
		{
			//########################################################################
			this->m_log->ERROR("UGrid " + this->m_name
									+ ": Attempted to access neighbors array of size "
									+ std::to_string(this->m_N) + " with index "
									+ std::to_string(index));
			//########################################################################
			if(m_neighbors.size() > 0)
			{
				//######################################################################
				this->m_log->INFO("UGrid " + this->m_name + ": Returning the element at index 0");
				//######################################################################
				return m_neighbors[0];
			}
			else
			{
				//######################################################################
				this->m_log->INFO("UGrid " + this->m_name + ": Returning empty neighbors array");
				//######################################################################
				return std::vector<size_t>(1,0);
			}
		}
		return m_kdt.getNeighbors(index);
	}

  // //  Setters
  // template<typename T>
  // void UGrid<T>::setDim(uint64_t dim)
  // {
	// 	if (dim != m_ugrid[0].size())
	// 	{
	// 		//########################################################################
	// 		this->m_log->WARN("UGrid " + this->m_name + ": New dimension "
	// 	           	+ std::to_string(dim) + " does not match array m_ugrid"
	// 						  + " of dimension " + std::to_string(this->m_dim));
	// 		//########################################################################
	// 	}
  //   this->m_dim = dim;
	// 	//##########################################################################
	// 	this->m_log->INFO("UGrid " + this->m_name + ": Setting dimension this->m_dim to "
	// 	           + std::to_string(this->m_dim));
	// 	//##########################################################################
	// 	m_searchFlag = -1;
  // }
  // template<typename T>
  // void UGrid<T>::setN(uint64_t N)
  // {
	// 	if (N != m_ugrid.size())
	// 	{
	// 		//########################################################################
	// 		this->m_log->WARN("UGrid " + this->m_name + ": New array size "
	// 							+ std::to_string(N) + " does not match array m_ugrid"
	// 							+ " of size " + std::to_string(this->m_N));
	// 		//########################################################################
	// 	}
  //   this->m_N = N;
	// 	//##########################################################################
	// 	this->m_log->INFO("UGrid " + this->m_name+ ": Setting number of elements this->m_N to "
	// 						 + std::to_string(this->m_N));
	// 	//##########################################################################
	// 	m_searchFlag = -1;
  // }
  template<typename T>
  void UGrid<T>::setUGrid(std::vector<std::vector<T>> ugrid)
  {
    m_ugrid = ugrid;
		this->m_N = ugrid.size();
		this->m_dim = ugrid[0].size();
		//##########################################################################
		this->m_log->INFO("UGrid " + this->m_name + ": Setting m_ugrid to array of size "
							 + std::to_string(this->m_N) + " with dimension "
							 + std::to_string(this->m_dim));
		//##########################################################################
    KDTree<T> kdt(std::make_shared<std::vector<std::vector<T>>>(m_ugrid));
    m_kdt = kdt;
		// //  generate KDTree
  	// KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<T>>, T>
		// kdt(this->m_dim, m_ugrid, 16);
    // kdt.index->buildIndex();
		// m_KDTree = std::make_shared<KDTreeVectorOfVectorsAdaptor<
		//                             std::vector<std::vector<T>>, T>>(kdt);
		// m_searchFlag = -1;
  }
  // template<typename T>
  // void UGrid<T>::setName(std::string name)
  // {
	// 	//##########################################################################
	// 	this->m_log->INFO("UGrid " + this->m_name + ": Renaming '"
	// 	           + this->m_name + "' to '" + name + "'");
  //   //##########################################################################
  //   this->m_name = name;
  // }
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  Access operators for ugrid
	//----------------------------------------------------------------------------
	template<typename T>
  T& UGrid<T>::operator()(const uint64_t i, const uint64_t j)
  {
		if (i >= this->m_N || j >= this->m_dim)
		{
			if(i >= this->m_N)
			{
				//######################################################################
				this->m_log->ERROR("UGrid " + this->m_name
										+ ": Attempted to access m_ugrid array of size "
										+ std::to_string(this->m_N) + " with index "
										+ std::to_string(i));
				//######################################################################
			}
			if(j >= this->m_dim)
			{
				//######################################################################
				this->m_log->ERROR("UGrid " + this->m_name
										+ ": Attempted to access m_ugrid array of dimension "
										+ std::to_string(this->m_dim) + " with index "
										+ std::to_string(j));
				//######################################################################
			}
			if(m_ugrid.size() > 0 && m_ugrid[0].size() > 0)
			{
				//######################################################################
				this->m_log->INFO("UGrid " + this->m_name + ": Returning the element at index (0,0)");
				//######################################################################
				return m_ugrid[0][0];
			}
			else
			{
				//######################################################################
				this->m_log->INFO("UGrid " + this->m_name + ": Terminating program");
				//######################################################################
				exit(0);
			}
		}
		m_searchFlag = -1;
    return m_ugrid[i][j];
  }
  template<typename T>
  const T& UGrid<T>::operator()(const uint64_t i, const uint64_t j) const
  {
		if (i >= this->m_N || j >= this->m_dim)
		{
			if(i >= this->m_N)
			{
				//######################################################################
				this->m_log->ERROR("UGrid " + this->m_name
										+ ": Attempted to access m_ugrid array of size "
										+ std::to_string(this->m_N) + " with index "
										+ std::to_string(i));
				//######################################################################
			}
			if(j >= this->m_dim)
			{
				//######################################################################
				this->m_log->ERROR("UGrid " + this->m_name
										+ ": Attempted to access m_ugrid array of dimension "
										+ std::to_string(this->m_dim) + " with index "
										+ std::to_string(j));
				//######################################################################
			}
			if(m_ugrid.size() > 0 && m_ugrid[0].size() > 0)
			{
				//######################################################################
				this->m_log->INFO("UGrid " + this->m_name + ": Returning the element at index (0,0)");
				//######################################################################
				return m_ugrid[0][0];
			}
			else
			{
				//######################################################################
				this->m_log->INFO("UGrid " + this->m_name + ": Terminating program");
				//######################################################################
				exit(0);
			}
		}
    return m_ugrid[i][j];
  }
	template<typename T>
  std::vector<T>& UGrid<T>::operator()(const uint64_t i)
  {
		if (i >= this->m_N )
		{
			//########################################################################
			this->m_log->ERROR("UGrid " + this->m_name
									+ ": Attempted to access m_ugrid array of size "
									+ std::to_string(this->m_N) + " with index "
									+ std::to_string(i));
			//########################################################################
			if(m_ugrid.size() > 0)
			{
				//######################################################################
				this->m_log->INFO("UGrid " + this->m_name + ": Returning the element at index 0");
				//######################################################################
				return m_ugrid[0];
			}
			else
			{
				//######################################################################
				this->m_log->INFO("UGrid " + this->m_name + ": Terminating program");
				//######################################################################
				exit(0);
			}
		}
		m_searchFlag = -1;
    return m_ugrid[i];
  }
  template<typename T>
  const std::vector<T>& UGrid<T>::operator()(const uint64_t i) const
  {
		if (i >= this->m_N )
		{
			//########################################################################
			this->m_log->ERROR("UGrid " + this->m_name
									+ ": Attempted to access m_ugrid array of size "
									+ std::to_string(this->m_N) + " with index "
									+ std::to_string(i));
			//########################################################################
			if(m_ugrid.size() > 0)
			{
				//######################################################################
				this->m_log->INFO("UGrid " + this->m_name + ": Returning the element at index 0");
				//######################################################################
				return m_ugrid[0];
			}
			else
			{
				//######################################################################
				this->m_log->INFO("UGrid " + this->m_name + ": Terminating program");
				//######################################################################
				exit(0);
			}
		}
    return m_ugrid[i];
  }
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  Points and projections
	//----------------------------------------------------------------------------
  template<typename T>
  const std::vector<T>& UGrid<T>::getPoint(const uint64_t i) const
  {
		if (i >= this->m_N )
		{
			//########################################################################
			this->m_log->ERROR("UGrid " + this->m_name
									+ ": Attempted to access m_ugrid array of size "
									+ std::to_string(this->m_N) + " with index "
									+ std::to_string(i));
			//########################################################################
			if(m_neighbors.size() > 0)
			{
				//######################################################################
				this->m_log->INFO("UGrid " + this->m_name + ": Returning the element at index 0");
				//######################################################################
				return m_ugrid[0];
			}
			else
			{
				//######################################################################
				this->m_log->INFO("UGrid " + this->m_name + ": Terminating program");
				//######################################################################
				exit(0);
			}
		}
    return m_ugrid[i];
  }
  template<typename T>
  std::vector<T> UGrid<T>::projection(uint64_t j)
  {
		std::vector<T> result(this->m_N);
		for (uint64_t i = 0; i < this->m_N; i++)
		{
			result[i] = m_ugrid[i][j];
		}
    return result;
  }
  template<typename T>
  void UGrid<T>::setPoint(uint64_t i, std::vector<T> point)
  {
		if (i >= this->m_N )
		{
			//########################################################################
			this->m_log->ERROR("UGrid " + this->m_name
									+ ": Attempted to access m_ugrid array of size "
									+ std::to_string(this->m_N) + " with index "
									+ std::to_string(i));
			//########################################################################
			return;
		}
    m_ugrid[i] = point;
		m_searchFlag = -1;
  }
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	KDTree methods
	//----------------------------------------------------------------------------
	template<typename T>
	void UGrid<T>::setupTree()
	{
    m_kdt = KDTree<T>(std::make_shared<std::vector<std::vector<T>>>(m_ugrid));
	}
	template<typename T>
  void UGrid<T>::queryNeighbors(uint64_t k)
  {
    return m_kdt.queryNeighbors(k);
  }
	//----------------------------------------------------------------------------
	//	Query the tree for neighbors of some point
	//----------------------------------------------------------------------------
	template<typename T>
  std::vector<size_t>
  UGrid<T>::queryNeighbors(const std::vector<T>& point,
		                       uint64_t k)
  {
    return m_kdt.queryNeighbors(point,k);
  }
	//----------------------------------------------------------------------------
	//	Query the tree for neighbors of some point
	//----------------------------------------------------------------------------
	template<typename T>
  std::vector<double>
  UGrid<T>::queryDistances(const std::vector<T>& point,
		                       uint64_t k)
  {
    return m_kdt.queryDistances(point,k);
  }
	//----------------------------------------------------------------------------
	//	Query the tree for neighbors of some set of points
	//----------------------------------------------------------------------------
	template<typename T>
  std::vector<std::vector<size_t>>
  UGrid<T>::queryNeighbors(const std::vector<std::vector<T>>& points,
		                       uint64_t k)
  {
    return m_kdt.queryNeighbors(points,k);
  }
	template<typename T>
  void UGrid<T>::queryRadius(double radius)
  {
    return m_kdt.queryRadius(radius);
  }
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Check the consistency of internal elements
	//----------------------------------------------------------------------------
	template<typename T>
	bool UGrid<T>::checkConsistency()
	{
		bool consistency = true;
		if(m_ugrid.size() != this->m_N)
		{
			consistency = false;
			//########################################################################
			this->m_log->WARN("UGrid " + this->m_name + ": Attribute this->m_N with value "
			           + std::to_string(this->m_N) + " does not match array m_ugrid"
								 + " with size " + std::to_string(m_ugrid.size()));
			//########################################################################
		}
		if(m_ugrid[0].size() != this->m_dim)
		{
			consistency = false;
			//########################################################################
			this->m_log->WARN("UGrid " + this->m_name + ": Attribute this->m_dim with value "
			           + std::to_string(this->m_dim) + " does not match array m_ugrid"
								 + " with dimension " + std::to_string(m_ugrid[0].size()));
			//########################################################################
		}
		return consistency;
	}
	//----------------------------------------------------------------------------
}
