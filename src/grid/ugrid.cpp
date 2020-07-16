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
	//----------------------------------------------------------------------------
  //  UGrid constructors
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Default constructor
  //    sets name = "default", and m_dim, m_N = 0
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid() : m_dim(0), m_N(0), m_name("default")
  {
		//##########################################################################
		m_log = std::make_shared<Log>();
		m_log->init("ET:UGrid:default", ".logs/ugrid_default.txt");
		m_log->TRACE("Unstructured Grid 'default' created at location "
		            + getMem(*this));
		//##########################################################################
		m_searchFlag = -1;
  }
	//----------------------------------------------------------------------------
  //  Default destructor
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::~UGrid()
  {
		//##########################################################################
		m_log->TRACE("Unstructured Grid '" + m_name
		            + "' destroyed at location " + getMem(*this));
		//##########################################################################
  }
	//----------------------------------------------------------------------------
  //  Various constructors taking in arguments for
	//		m_dim, m_name, m_N, m_ugrid and m_log.
  //----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Constructor with dimension
	//----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(uint64_t dim) : m_dim(dim), m_name("default")
  {
		//##########################################################################
		m_log = std::make_shared<Log>();
		m_log->init("ET:UGrid:default", ".logs/ugrid_default.txt");
		m_log->TRACE("Unstructured Grid 'default' created at location "
		            + getMem(*this));
		//##########################################################################
		m_searchFlag = -1;
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  Constructor with name and dimension
	//----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(std::string name, uint64_t dim) : m_name(name), m_dim(dim)
  {
		//##########################################################################
		m_log = std::make_shared<Log>();
		m_log->init("ET:UGrid:" + m_name, ".logs/ugrid_" + m_name
		           + ".txt");
		m_log->TRACE("Unstructured Grid '" + m_name
		            + "' created at location " + getMem(*this));
		//##########################################################################
		m_searchFlag = -1;
  }
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  Constructor with dimension and number of points
	//----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(uint64_t dim, uint64_t N) : m_dim(dim), m_N(N), m_name("default")
  {
		//##########################################################################
		m_log = std::make_shared<Log>();
		m_log->init("ET:UGrid:default", ".logs/ugrid_default.txt");
		m_log->TRACE("Unstructured Grid 'default' created at location "
		            + getMem(*this));
		//##########################################################################
		m_searchFlag = -1;
  }
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  Constructor with name, dimension and number of points
	//----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(std::string name, uint64_t dim, uint64_t N)
  : m_name(name), m_dim(dim), m_N(N)
  {
		//##########################################################################
		m_log = std::make_shared<Log>();
		m_log->init("ET:UGrid:" + m_name, ".logs/ugrid_" + m_name + ".txt");
		m_log->TRACE("Unstructured Grid '" + m_name
		            + "' created at location " + getMem(*this));
	  //##########################################################################
		m_searchFlag = -1;
  }
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Constructor with one dimensional array
	//----------------------------------------------------------------------------
	template<typename T>
	UGrid<T>::UGrid(std::vector<T> ugrid) : m_name("default")
	{
		m_N = ugrid.size();
		m_ugrid.resize(m_N);
		m_dim = 1;
		for (uint32_t i = 0; i < m_N; i++)
		{
			std::vector<T> temp = {ugrid[i]};
			m_ugrid[i] = temp;
		}
		//##########################################################################
		m_log = std::make_shared<Log>();
		m_log->init("ET:UGrid:default", ".logs/ugrid_default.txt");
		m_log->TRACE("Unstructured Grid 'default' created at location "
		            + getMem(*this));
		//##########################################################################
		//  generate KDTree
  	setupTree();
		m_searchFlag = -1;
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Constructor with vector of vectors array
	//----------------------------------------------------------------------------
	template<typename T>
	UGrid<T>::UGrid(std::vector<std::vector<T>> ugrid) : m_name("default")
	{
		m_N = ugrid.size();
		m_dim = ugrid[0].size();
		m_ugrid = ugrid;
		//##########################################################################
		m_log = std::make_shared<Log>();
		m_log->init("ET:UGrid:default", ".logs/ugrid_default.txt");
		m_log->TRACE("Unstructured Grid 'default' created at location "
		            + getMem(*this));
	  //##########################################################################
		//  generate KDTree
  	setupTree();
		m_searchFlag = -1;
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Constructors with shared loggers
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Constructor with logger
	//----------------------------------------------------------------------------
	template<typename T>
	UGrid<T>::UGrid(std::shared_ptr<Log> log) : m_dim(0), m_N(0), m_name("default")
	{
		//##########################################################################
		m_log = log;
		m_log->TRACE("Unstructured Grid 'default' created at location "
								+ getMem(*this));
		m_log->INFO("Logger passed to Unstructured Grid 'default'");
		//##########################################################################
		m_searchFlag = -1;
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Constructor with dimension and logger
	//----------------------------------------------------------------------------
	template<typename T>
	UGrid<T>::UGrid(uint64_t dim, std::shared_ptr<Log> log)
	: m_dim(dim), m_name("default")
	{
		//##########################################################################
		m_log = log;
		m_log->TRACE("Unstructured Grid 'default' created at location "
								+ getMem(*this));
		m_log->INFO("Logger passed to Unstructured Grid 'default'");
		//##########################################################################
		m_searchFlag = -1;
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  Constructor with name, dimension and shared logger
	//----------------------------------------------------------------------------
	template<typename T>
	UGrid<T>::UGrid(std::string name, uint64_t dim, std::shared_ptr<Log> log)
	: m_dim(dim), m_name(name)
	{
		//##########################################################################
		m_log = log;
		m_log->TRACE("Unstructured Grid '" + m_name
								+ "' created at location " + getMem(*this));
		m_log->INFO("Logger passed to Unstructured Grid '" + m_name + "'");
		//##########################################################################
		m_searchFlag = -1;
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  Constructor with dimension, number of points and shared logger
	//----------------------------------------------------------------------------
	template<typename T>
	UGrid<T>::UGrid(uint64_t dim, uint64_t N, std::shared_ptr<Log> log)
	: m_dim(dim), m_N(N), m_name("default")
	{
		//##########################################################################
		m_log = log;
		m_log->TRACE("Unstructured Grid 'default' created at location "
								+ getMem(*this));
		m_log->INFO("Logger passed to Unstructured Grid 'default'");
		//##########################################################################
		m_searchFlag = -1;
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Constructor with name, dimension, number of points and shared logger
	//----------------------------------------------------------------------------
	template<typename T>
	UGrid<T>::UGrid(std::string name, uint64_t dim,
		              uint64_t N, std::shared_ptr<Log> log)
	: m_name(name), m_dim(dim), m_N(N)
	{
		//##########################################################################
		m_log = log;
		m_log->TRACE("Unstructured Grid '" + m_name
								+ "' created at location " + getMem(*this));
		m_log->INFO("Logger passed to Unstructured Grid '" + m_name + "'");
		//##########################################################################
		m_searchFlag = -1;
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Constructor with one dimensional array and shared logger
	//----------------------------------------------------------------------------
	template<typename T>
	UGrid<T>::UGrid(std::vector<T> ugrid, std::shared_ptr<Log> log)
	: m_name("default")
	{
		m_N = ugrid.size();
		m_ugrid.resize(m_N);
		m_dim = 1;
		for (uint32_t i = 0; i < m_N; i++)
		{
			std::vector<T> temp = {ugrid[i]};
			m_ugrid[i] = temp;
		}
		//##########################################################################
		m_log = log;
		m_log->TRACE("Unstructured Grid 'default' created at location "
								+ getMem(*this));
		m_log->INFO("Logger passed to Unstructured Grid 'default'");
		//##########################################################################
		//  generate KDTree
  	setupTree();
		m_searchFlag = -1;
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Constructor with vector of vectors array
	//----------------------------------------------------------------------------
	template<typename T>
	UGrid<T>::UGrid(std::vector<std::vector<T>> ugrid, std::shared_ptr<Log> log)
	: m_name("default")
	{
		m_N = ugrid.size();
		m_dim = ugrid[0].size();
		m_ugrid = ugrid;
		//##########################################################################
		m_log = log;
		m_log->TRACE("Unstructured Grid 'default' created at location "
								+ getMem(*this));
		m_log->INFO("Logger passed to Unstructured Grid 'default'");
		//##########################################################################
		//  generate KDTree
  	setupTree();
		m_searchFlag = -1;
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  Getters and Setters
	//----------------------------------------------------------------------------
  template<typename T>
  const uint64_t UGrid<T>::getDim()
  {
    return m_dim;
  }
  template<typename T>
  const uint64_t UGrid<T>::getN()
  {
    return m_N;
  }
  template<typename T>
  std::vector<std::vector<T> > UGrid<T>::getUGrid()
  {
    return m_ugrid;
  }
  template<typename T>
  const std::string UGrid<T>::getName()
  {
    return m_name;
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
		if (index > m_N)
		{
			//########################################################################
			m_log->ERROR("UGrid " + m_name
									+ ": Attempted to access neighbors array of size "
									+ std::to_string(m_N) + " with index "
									+ std::to_string(index));
			//########################################################################
			if(m_neighbors.size() > 0)
			{
				//######################################################################
				m_log->INFO("UGrid " + m_name + ": Returning the element at index 0");
				//######################################################################
				return m_neighbors[0];
			}
			else
			{
				//######################################################################
				m_log->INFO("UGrid " + m_name + ": Returning empty neighbors array");
				//######################################################################
				return std::vector<size_t>(1,0);
			}
		}
		return m_kdt.getNeighbors(index);
	}
	template<typename T>
	std::shared_ptr<Log> UGrid<T>::getLogger()
	{
		return m_log;
	}
  //  Setters
  template<typename T>
  void UGrid<T>::setDim(uint64_t dim)
  {
		if (dim != m_ugrid[0].size())
		{
			//########################################################################
			m_log->WARN("UGrid " + m_name + ": New dimension "
		           	+ std::to_string(dim) + " does not match array m_ugrid"
							  + " of dimension " + std::to_string(m_dim));
			//########################################################################
		}
    m_dim = dim;
		//##########################################################################
		m_log->INFO("UGrid " + m_name + ": Setting dimension m_dim to "
		           + std::to_string(m_dim));
		//##########################################################################
		m_searchFlag = -1;
  }
  template<typename T>
  void UGrid<T>::setN(uint64_t N)
  {
		if (N != m_ugrid.size())
		{
			//########################################################################
			m_log->WARN("UGrid " + m_name + ": New array size "
								+ std::to_string(N) + " does not match array m_ugrid"
								+ " of size " + std::to_string(m_N));
			//########################################################################
		}
    m_N = N;
		//##########################################################################
		m_log->INFO("UGrid " + m_name+ ": Setting number of elements m_N to "
							 + std::to_string(m_N));
		//##########################################################################
		m_searchFlag = -1;
  }
  template<typename T>
  void UGrid<T>::setUGrid(std::vector<std::vector<T>> ugrid)
  {
    m_ugrid = ugrid;
		m_N = ugrid.size();
		m_dim = ugrid[0].size();
		//##########################################################################
		m_log->INFO("UGrid " + m_name + ": Setting m_ugrid to array of size "
							 + std::to_string(m_N) + " with dimension "
							 + std::to_string(m_dim));
		//##########################################################################
    KDTree<T> kdt(std::make_shared<std::vector<std::vector<T>>>(m_ugrid));
    m_kdt = kdt;
		// //  generate KDTree
  	// KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<T>>, T>
		// kdt(m_dim, m_ugrid, 16);
    // kdt.index->buildIndex();
		// m_KDTree = std::make_shared<KDTreeVectorOfVectorsAdaptor<
		//                             std::vector<std::vector<T>>, T>>(kdt);
		// m_searchFlag = -1;
  }
  template<typename T>
  void UGrid<T>::setName(std::string name)
  {
		//##########################################################################
		m_log->INFO("UGrid " + m_name + ": Renaming '"
		           + m_name + "' to '" + name + "'");
    //##########################################################################
    m_name = name;
  }
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  Access operators for ugrid
	//----------------------------------------------------------------------------
	template<typename T>
  T& UGrid<T>::operator()(const uint64_t i, const uint64_t j)
  {
		if (i >= m_N || j >= m_dim)
		{
			if(i >= m_N)
			{
				//######################################################################
				m_log->ERROR("UGrid " + m_name
										+ ": Attempted to access m_ugrid array of size "
										+ std::to_string(m_N) + " with index "
										+ std::to_string(i));
				//######################################################################
			}
			if(j >= m_dim)
			{
				//######################################################################
				m_log->ERROR("UGrid " + m_name
										+ ": Attempted to access m_ugrid array of dimension "
										+ std::to_string(m_dim) + " with index "
										+ std::to_string(j));
				//######################################################################
			}
			if(m_ugrid.size() > 0 && m_ugrid[0].size() > 0)
			{
				//######################################################################
				m_log->INFO("UGrid " + m_name + ": Returning the element at index (0,0)");
				//######################################################################
				return m_ugrid[0][0];
			}
			else
			{
				//######################################################################
				m_log->INFO("UGrid " + m_name + ": Terminating program");
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
		if (i >= m_N || j >= m_dim)
		{
			if(i >= m_N)
			{
				//######################################################################
				m_log->ERROR("UGrid " + m_name
										+ ": Attempted to access m_ugrid array of size "
										+ std::to_string(m_N) + " with index "
										+ std::to_string(i));
				//######################################################################
			}
			if(j >= m_dim)
			{
				//######################################################################
				m_log->ERROR("UGrid " + m_name
										+ ": Attempted to access m_ugrid array of dimension "
										+ std::to_string(m_dim) + " with index "
										+ std::to_string(j));
				//######################################################################
			}
			if(m_ugrid.size() > 0 && m_ugrid[0].size() > 0)
			{
				//######################################################################
				m_log->INFO("UGrid " + m_name + ": Returning the element at index (0,0)");
				//######################################################################
				return m_ugrid[0][0];
			}
			else
			{
				//######################################################################
				m_log->INFO("UGrid " + m_name + ": Terminating program");
				//######################################################################
				exit(0);
			}
		}
    return m_ugrid[i][j];
  }
	template<typename T>
  std::vector<T>& UGrid<T>::operator()(const uint64_t i)
  {
		if (i >= m_N )
		{
			//########################################################################
			m_log->ERROR("UGrid " + m_name
									+ ": Attempted to access m_ugrid array of size "
									+ std::to_string(m_N) + " with index "
									+ std::to_string(i));
			//########################################################################
			if(m_ugrid.size() > 0)
			{
				//######################################################################
				m_log->INFO("UGrid " + m_name + ": Returning the element at index 0");
				//######################################################################
				return m_ugrid[0];
			}
			else
			{
				//######################################################################
				m_log->INFO("UGrid " + m_name + ": Terminating program");
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
		if (i >= m_N )
		{
			//########################################################################
			m_log->ERROR("UGrid " + m_name
									+ ": Attempted to access m_ugrid array of size "
									+ std::to_string(m_N) + " with index "
									+ std::to_string(i));
			//########################################################################
			if(m_ugrid.size() > 0)
			{
				//######################################################################
				m_log->INFO("UGrid " + m_name + ": Returning the element at index 0");
				//######################################################################
				return m_ugrid[0];
			}
			else
			{
				//######################################################################
				m_log->INFO("UGrid " + m_name + ": Terminating program");
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
		if (i >= m_N )
		{
			//########################################################################
			m_log->ERROR("UGrid " + m_name
									+ ": Attempted to access m_ugrid array of size "
									+ std::to_string(m_N) + " with index "
									+ std::to_string(i));
			//########################################################################
			if(m_neighbors.size() > 0)
			{
				//######################################################################
				m_log->INFO("UGrid " + m_name + ": Returning the element at index 0");
				//######################################################################
				return m_ugrid[0];
			}
			else
			{
				//######################################################################
				m_log->INFO("UGrid " + m_name + ": Terminating program");
				//######################################################################
				exit(0);
			}
		}
    return m_ugrid[i];
  }
  template<typename T>
  std::vector<T> UGrid<T>::projection(uint64_t j)
  {
		std::vector<T> result(m_N);
		for (uint64_t i = 0; i < m_N; i++)
		{
			result[i] = m_ugrid[i][j];
		}
    return result;
  }
  template<typename T>
  void UGrid<T>::setPoint(uint64_t i, std::vector<T> point)
  {
		if (i >= m_N )
		{
			//########################################################################
			m_log->ERROR("UGrid " + m_name
									+ ": Attempted to access m_ugrid array of size "
									+ std::to_string(m_N) + " with index "
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
		if(m_ugrid.size() != m_N)
		{
			consistency = false;
			//########################################################################
			m_log->WARN("UGrid " + m_name + ": Attribute m_N with value "
			           + std::to_string(m_N) + " does not match array m_ugrid"
								 + " with size " + std::to_string(m_ugrid.size()));
			//########################################################################
		}
		if(m_ugrid[0].size() != m_dim)
		{
			consistency = false;
			//########################################################################
			m_log->WARN("UGrid " + m_name + ": Attribute m_dim with value "
			           + std::to_string(m_dim) + " does not match array m_ugrid"
								 + " with dimension " + std::to_string(m_ugrid[0].size()));
			//########################################################################
		}
		return consistency;
	}
	//----------------------------------------------------------------------------
}
