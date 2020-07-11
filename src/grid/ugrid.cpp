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
  //    sets name = "default", and _dim, _N = 0
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid() : _dim(0), _N(0), _name("default")
  {
		//##########################################################################
		_log = std::make_shared<Log>();
		_log->init("ET:UGrid:default", ".logs/ugrid_default.txt");
		_log->TRACE("Unstructured Grid 'default' created at location "
		            + getMem(*this));
		//##########################################################################
		_searchFlag = -1;
  }
	//----------------------------------------------------------------------------
  //  Default destructor
  //----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::~UGrid()
  {
		//##########################################################################
		_log->TRACE("Unstructured Grid '" + _name
		            + "' destroyed at location " + getMem(*this));
		//##########################################################################
  }
	//----------------------------------------------------------------------------
  //  Various constructors taking in arguments for
	//		_dim, _name, _N, _ugrid and _log.
  //----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Constructor with dimension
	//----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(uint64_t dim) : _dim(dim), _name("default")
  {
		//##########################################################################
		_log = std::make_shared<Log>();
		_log->init("ET:UGrid:default", ".logs/ugrid_default.txt");
		_log->TRACE("Unstructured Grid 'default' created at location "
		            + getMem(*this));
		//##########################################################################
		_searchFlag = -1;
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  Constructor with name and dimension
	//----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(std::string name, uint64_t dim) : _name(name), _dim(dim)
  {
		//##########################################################################
		_log = std::make_shared<Log>();
		_log->init("ET:UGrid:" + _name, ".logs/ugrid_" + _name
		           + ".txt");
		_log->TRACE("Unstructured Grid '" + _name
		            + "' created at location " + getMem(*this));
		//##########################################################################
		_searchFlag = -1;
  }
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  Constructor with dimension and number of points
	//----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(uint64_t dim, uint64_t N) : _dim(dim), _N(N), _name("default")
  {
		//##########################################################################
		_log = std::make_shared<Log>();
		_log->init("ET:UGrid:default", ".logs/ugrid_default.txt");
		_log->TRACE("Unstructured Grid 'default' created at location "
		            + getMem(*this));
		//##########################################################################
		_searchFlag = -1;
  }
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  Constructor with name, dimension and number of points
	//----------------------------------------------------------------------------
  template<typename T>
  UGrid<T>::UGrid(std::string name, uint64_t dim, uint64_t N)
  : _name(name), _dim(dim), _N(N)
  {
		//##########################################################################
		_log = std::make_shared<Log>();
		_log->init("ET:UGrid:" + _name, ".logs/ugrid_" + _name + ".txt");
		_log->TRACE("Unstructured Grid '" + _name
		            + "' created at location " + getMem(*this));
	  //##########################################################################
		_searchFlag = -1;
  }
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Constructor with one dimensional array
	//----------------------------------------------------------------------------
	template<typename T>
	UGrid<T>::UGrid(std::vector<T> ugrid) : _name("default")
	{
		_N = ugrid.size();
		_ugrid.resize(_N);
		_dim = 1;
		for (uint32_t i = 0; i < _N; i++)
		{
			std::vector<T> temp = {ugrid[i]};
			_ugrid[i] = temp;
		}
		//##########################################################################
		_log = std::make_shared<Log>();
		_log->init("ET:UGrid:default", ".logs/ugrid_default.txt");
		_log->TRACE("Unstructured Grid 'default' created at location "
		            + getMem(*this));
		//##########################################################################
		//  generate kdtree
  	setupTree();
		_searchFlag = -1;
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Constructor with vector of vectors array
	//----------------------------------------------------------------------------
	template<typename T>
	UGrid<T>::UGrid(std::vector<std::vector<T>> ugrid) : _name("default")
	{
		_N = ugrid.size();
		_dim = ugrid[0].size();
		_ugrid = ugrid;
		//##########################################################################
		_log = std::make_shared<Log>();
		_log->init("ET:UGrid:default", ".logs/ugrid_default.txt");
		_log->TRACE("Unstructured Grid 'default' created at location "
		            + getMem(*this));
	  //##########################################################################
		//  generate kdtree
  	setupTree();
		_searchFlag = -1;
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Constructors with shared loggers
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Constructor with logger
	//----------------------------------------------------------------------------
	template<typename T>
	UGrid<T>::UGrid(std::shared_ptr<Log> log) : _dim(0), _N(0), _name("default")
	{
		//##########################################################################
		_log = log;
		_log->TRACE("Unstructured Grid 'default' created at location "
								+ getMem(*this));
		_log->INFO("Logger passed to Unstructured Grid 'default'");
		//##########################################################################
		_searchFlag = -1;
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Constructor with dimension and logger
	//----------------------------------------------------------------------------
	template<typename T>
	UGrid<T>::UGrid(uint64_t dim, std::shared_ptr<Log> log)
	: _dim(dim), _name("default")
	{
		//##########################################################################
		_log = log;
		_log->TRACE("Unstructured Grid 'default' created at location "
								+ getMem(*this));
		_log->INFO("Logger passed to Unstructured Grid 'default'");
		//##########################################################################
		_searchFlag = -1;
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  Constructor with name, dimension and shared logger
	//----------------------------------------------------------------------------
	template<typename T>
	UGrid<T>::UGrid(std::string name, uint64_t dim, std::shared_ptr<Log> log)
	: _dim(dim), _name(name)
	{
		//##########################################################################
		_log = log;
		_log->TRACE("Unstructured Grid '" + _name
								+ "' created at location " + getMem(*this));
		_log->INFO("Logger passed to Unstructured Grid '" + _name + "'");
		//##########################################################################
		_searchFlag = -1;
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  Constructor with dimension, number of points and shared logger
	//----------------------------------------------------------------------------
	template<typename T>
	UGrid<T>::UGrid(uint64_t dim, uint64_t N, std::shared_ptr<Log> log)
	: _dim(dim), _N(N), _name("default")
	{
		//##########################################################################
		_log = log;
		_log->TRACE("Unstructured Grid 'default' created at location "
								+ getMem(*this));
		_log->INFO("Logger passed to Unstructured Grid 'default'");
		//##########################################################################
		_searchFlag = -1;
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Constructor with name, dimension, number of points and shared logger
	//----------------------------------------------------------------------------
	template<typename T>
	UGrid<T>::UGrid(std::string name, uint64_t dim,
		              uint64_t N, std::shared_ptr<Log> log)
	: _name(name), _dim(dim), _N(N)
	{
		//##########################################################################
		_log = log;
		_log->TRACE("Unstructured Grid '" + _name
								+ "' created at location " + getMem(*this));
		_log->INFO("Logger passed to Unstructured Grid '" + _name + "'");
		//##########################################################################
		_searchFlag = -1;
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Constructor with one dimensional array and shared logger
	//----------------------------------------------------------------------------
	template<typename T>
	UGrid<T>::UGrid(std::vector<T> ugrid, std::shared_ptr<Log> log)
	: _name("default")
	{
		_N = ugrid.size();
		_ugrid.resize(_N);
		_dim = 1;
		for (uint32_t i = 0; i < _N; i++)
		{
			std::vector<T> temp = {ugrid[i]};
			_ugrid[i] = temp;
		}
		//##########################################################################
		_log = log;
		_log->TRACE("Unstructured Grid 'default' created at location "
								+ getMem(*this));
		_log->INFO("Logger passed to Unstructured Grid 'default'");
		//##########################################################################
		//  generate kdtree
  	setupTree();
		_searchFlag = -1;
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Constructor with vector of vectors array
	//----------------------------------------------------------------------------
	template<typename T>
	UGrid<T>::UGrid(std::vector<std::vector<T>> ugrid, std::shared_ptr<Log> log)
	: _name("default")
	{
		_N = ugrid.size();
		_dim = ugrid[0].size();
		_ugrid = ugrid;
		//##########################################################################
		_log = log;
		_log->TRACE("Unstructured Grid 'default' created at location "
								+ getMem(*this));
		_log->INFO("Logger passed to Unstructured Grid 'default'");
		//##########################################################################
		//  generate kdtree
  	setupTree();
		_searchFlag = -1;
	}
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  Getters and Setters
	//----------------------------------------------------------------------------
  template<typename T>
  const uint64_t UGrid<T>::getDim()
  {
    return _dim;
  }
  template<typename T>
  const uint64_t UGrid<T>::getN()
  {
    return _N;
  }
  template<typename T>
  std::vector<std::vector<T> > UGrid<T>::getUGrid()
  {
    return _ugrid;
  }
  template<typename T>
  const std::string UGrid<T>::getName()
  {
    return _name;
  }
	template<typename T>
	std::vector<std::vector<size_t>> UGrid<T>::getNeighbors()
	{
		return _kdt.getNeighbors();
	}
	template<typename T>
	std::vector<std::vector<double>> UGrid<T>::getDistances()
	{
		return _kdt.getDistances();
	}
	template<typename T>
	std::vector<std::vector<size_t>> UGrid<T>::getNeighborsRadius()
	{
		return _kdt.getNeighborsRadius();
	}
	template<typename T>
	std::vector<std::vector<double>> UGrid<T>::getDistancesRadius()
	{
		return _kdt.getDistancesRadius();
	}
	template<typename T>
	std::vector<size_t> UGrid<T>::getNeighbors(uint64_t index)
	{
		if (index > _N)
		{
			//########################################################################
			_log->ERROR("UGrid " + _name
									+ ": Attempted to access neighbors array of size "
									+ std::to_string(_N) + " with index "
									+ std::to_string(index));
			//########################################################################
			if(_neighbors.size() > 0)
			{
				//######################################################################
				_log->INFO("UGrid " + _name + ": Returning the element at index 0");
				//######################################################################
				return _neighbors[0];
			}
			else
			{
				//######################################################################
				_log->INFO("UGrid " + _name + ": Returning empty neighbors array");
				//######################################################################
				return std::vector<size_t>(1,0);
			}
		}
		return _kdt.getNeighbors(index);
	}
	template<typename T>
	std::shared_ptr<Log> UGrid<T>::getLogger()
	{
		return _log;
	}
  //  Setters
  template<typename T>
  void UGrid<T>::setDim(uint64_t dim)
  {
		if (dim != _ugrid[0].size())
		{
			//########################################################################
			_log->WARN("UGrid " + _name + ": New dimension "
		           	+ std::to_string(dim) + " does not match array _ugrid"
							  + " of dimension " + std::to_string(_dim));
			//########################################################################
		}
    _dim = dim;
		//##########################################################################
		_log->INFO("UGrid " + _name + ": Setting dimension _dim to "
		           + std::to_string(_dim));
		//##########################################################################
		_searchFlag = -1;
  }
  template<typename T>
  void UGrid<T>::setN(uint64_t N)
  {
		if (N != _ugrid.size())
		{
			//########################################################################
			_log->WARN("UGrid " + _name + ": New array size "
								+ std::to_string(N) + " does not match array _ugrid"
								+ " of size " + std::to_string(_N));
			//########################################################################
		}
    _N = N;
		//##########################################################################
		_log->INFO("UGrid " + _name+ ": Setting number of elements _N to "
							 + std::to_string(_N));
		//##########################################################################
		_searchFlag = -1;
  }
  template<typename T>
  void UGrid<T>::setUGrid(std::vector<std::vector<T>> ugrid)
  {
    _ugrid = ugrid;
		_N = ugrid.size();
		_dim = ugrid[0].size();
		//##########################################################################
		_log->INFO("UGrid " + _name + ": Setting _ugrid to array of size "
							 + std::to_string(_N) + " with dimension "
							 + std::to_string(_dim));
		//##########################################################################
    kdTree<T> kdt(std::make_shared<std::vector<std::vector<T>>>(_ugrid));
    _kdt = kdt;
		// //  generate kdtree
  	// KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<T>>, T>
		// kdt(_dim, _ugrid, 16);
    // kdt.index->buildIndex();
		// _kdtree = std::make_shared<KDTreeVectorOfVectorsAdaptor<
		//                             std::vector<std::vector<T>>, T>>(kdt);
		// _searchFlag = -1;
  }
  template<typename T>
  void UGrid<T>::setName(std::string name)
  {
		//##########################################################################
		_log->INFO("UGrid " + _name + ": Renaming '"
		           + _name + "' to '" + name + "'");
    //##########################################################################
    _name = name;
  }
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  Access operators for ugrid
	//----------------------------------------------------------------------------
	template<typename T>
  T& UGrid<T>::operator()(const uint64_t i, const uint64_t j)
  {
		if (i >= _N || j >= _dim)
		{
			if(i >= _N)
			{
				//######################################################################
				_log->ERROR("UGrid " + _name
										+ ": Attempted to access _ugrid array of size "
										+ std::to_string(_N) + " with index "
										+ std::to_string(i));
				//######################################################################
			}
			if(j >= _dim)
			{
				//######################################################################
				_log->ERROR("UGrid " + _name
										+ ": Attempted to access _ugrid array of dimension "
										+ std::to_string(_dim) + " with index "
										+ std::to_string(j));
				//######################################################################
			}
			if(_ugrid.size() > 0 && _ugrid[0].size() > 0)
			{
				//######################################################################
				_log->INFO("UGrid " + _name + ": Returning the element at index (0,0)");
				//######################################################################
				return _ugrid[0][0];
			}
			else
			{
				//######################################################################
				_log->INFO("UGrid " + _name + ": Terminating program");
				//######################################################################
				exit(0);
			}
		}
		_searchFlag = -1;
    return _ugrid[i][j];
  }
  template<typename T>
  const T& UGrid<T>::operator()(const uint64_t i, const uint64_t j) const
  {
		if (i >= _N || j >= _dim)
		{
			if(i >= _N)
			{
				//######################################################################
				_log->ERROR("UGrid " + _name
										+ ": Attempted to access _ugrid array of size "
										+ std::to_string(_N) + " with index "
										+ std::to_string(i));
				//######################################################################
			}
			if(j >= _dim)
			{
				//######################################################################
				_log->ERROR("UGrid " + _name
										+ ": Attempted to access _ugrid array of dimension "
										+ std::to_string(_dim) + " with index "
										+ std::to_string(j));
				//######################################################################
			}
			if(_ugrid.size() > 0 && _ugrid[0].size() > 0)
			{
				//######################################################################
				_log->INFO("UGrid " + _name + ": Returning the element at index (0,0)");
				//######################################################################
				return _ugrid[0][0];
			}
			else
			{
				//######################################################################
				_log->INFO("UGrid " + _name + ": Terminating program");
				//######################################################################
				exit(0);
			}
		}
    return _ugrid[i][j];
  }
	template<typename T>
  std::vector<T>& UGrid<T>::operator()(const uint64_t i)
  {
		if (i >= _N )
		{
			//########################################################################
			_log->ERROR("UGrid " + _name
									+ ": Attempted to access _ugrid array of size "
									+ std::to_string(_N) + " with index "
									+ std::to_string(i));
			//########################################################################
			if(_ugrid.size() > 0)
			{
				//######################################################################
				_log->INFO("UGrid " + _name + ": Returning the element at index 0");
				//######################################################################
				return _ugrid[0];
			}
			else
			{
				//######################################################################
				_log->INFO("UGrid " + _name + ": Terminating program");
				//######################################################################
				exit(0);
			}
		}
		_searchFlag = -1;
    return _ugrid[i];
  }
  template<typename T>
  const std::vector<T>& UGrid<T>::operator()(const uint64_t i) const
  {
		if (i >= _N )
		{
			//########################################################################
			_log->ERROR("UGrid " + _name
									+ ": Attempted to access _ugrid array of size "
									+ std::to_string(_N) + " with index "
									+ std::to_string(i));
			//########################################################################
			if(_ugrid.size() > 0)
			{
				//######################################################################
				_log->INFO("UGrid " + _name + ": Returning the element at index 0");
				//######################################################################
				return _ugrid[0];
			}
			else
			{
				//######################################################################
				_log->INFO("UGrid " + _name + ": Terminating program");
				//######################################################################
				exit(0);
			}
		}
    return _ugrid[i];
  }
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//  Points and projections
	//----------------------------------------------------------------------------
  template<typename T>
  const std::vector<T>& UGrid<T>::getPoint(const uint64_t i) const
  {
		if (i >= _N )
		{
			//########################################################################
			_log->ERROR("UGrid " + _name
									+ ": Attempted to access _ugrid array of size "
									+ std::to_string(_N) + " with index "
									+ std::to_string(i));
			//########################################################################
			if(_neighbors.size() > 0)
			{
				//######################################################################
				_log->INFO("UGrid " + _name + ": Returning the element at index 0");
				//######################################################################
				return _ugrid[0];
			}
			else
			{
				//######################################################################
				_log->INFO("UGrid " + _name + ": Terminating program");
				//######################################################################
				exit(0);
			}
		}
    return _ugrid[i];
  }
  template<typename T>
  std::vector<T> UGrid<T>::projection(uint64_t j)
  {
		std::vector<T> result(_N);
		for (uint64_t i = 0; i < _N; i++)
		{
			result[i] = _ugrid[i][j];
		}
    return result;
  }
  template<typename T>
  void UGrid<T>::setPoint(uint64_t i, std::vector<T> point)
  {
		if (i >= _N )
		{
			//########################################################################
			_log->ERROR("UGrid " + _name
									+ ": Attempted to access _ugrid array of size "
									+ std::to_string(_N) + " with index "
									+ std::to_string(i));
			//########################################################################
			return;
		}
    _ugrid[i] = point;
		_searchFlag = -1;
  }
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	KDTree methods
	//----------------------------------------------------------------------------
	template<typename T>
	void UGrid<T>::setupTree()
	{
    kdTree<T> kdt(std::make_shared<std::vector<std::vector<T>>>(_ugrid));
    _kdt = kdt;
	}
	template<typename T>
  void UGrid<T>::queryNeighbors(uint64_t k)
  {
    return _kdt.queryNeighbors(k);
  }
	//----------------------------------------------------------------------------
	//	Query the tree for neighbors of some point
	//----------------------------------------------------------------------------
	template<typename T>
  std::vector<size_t>
  UGrid<T>::queryNeighbors(const std::vector<T>& point,
		                       uint64_t k)
  {
    return _kdt.queryNeighbors(point,k);
  }
	//----------------------------------------------------------------------------
	//	Query the tree for neighbors of some point
	//----------------------------------------------------------------------------
	template<typename T>
  std::vector<double>
  UGrid<T>::queryDistances(const std::vector<T>& point,
		                       uint64_t k)
  {
    return _kdt.queryDistances(point,k);
  }
	//----------------------------------------------------------------------------
	//	Query the tree for neighbors of some set of points
	//----------------------------------------------------------------------------
	template<typename T>
  std::vector<std::vector<size_t>>
  UGrid<T>::queryNeighbors(const std::vector<std::vector<T>>& points,
		                       uint64_t k)
  {
    return _kdt.queryNeighbors(points,k);
  }
	template<typename T>
  void UGrid<T>::queryRadius(double radius)
  {
    return _kdt.queryRadius(radius);
  }
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	//	Check the consistency of internal elements
	//----------------------------------------------------------------------------
	template<typename T>
	bool UGrid<T>::checkConsistency()
	{
		bool consistency = true;
		if(_ugrid.size() != _N)
		{
			consistency = false;
			//########################################################################
			_log->WARN("UGrid " + _name + ": Attribute _N with value "
			           + std::to_string(_N) + " does not match array _ugrid"
								 + " with size " + std::to_string(_ugrid.size()));
			//########################################################################
		}
		if(_ugrid[0].size() != _dim)
		{
			consistency = false;
			//########################################################################
			_log->WARN("UGrid " + _name + ": Attribute _dim with value "
			           + std::to_string(_dim) + " does not match array _ugrid"
								 + " with dimension " + std::to_string(_ugrid[0].size()));
			//########################################################################
		}
		return consistency;
	}
	//----------------------------------------------------------------------------
}
