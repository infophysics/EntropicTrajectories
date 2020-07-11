//------------------------------------------------------------------------------
//  kdtree.h
//  The Entropic Trajectories Framework
//  -----------------------------------
//  Copyright (C) [2020] by [N. Carrara, F. Costa, P. Pessoa]
//  [ncarrara@albany.edu,felipecosta.physics@gmail.com,
//    pedroh.pessoa100@gmail.com]
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

#include "kdtree.h"

namespace ET
{
  template<typename T>
  kdTree<T>::kdTree() : _N(0), _dim(0), _name("default")
  {
    //##########################################################################
		_log = std::make_shared<Log>();
		_log->init("ET:kdTree:default", ".logs/kdtree_default.txt");
		_log->TRACE("kdTree 'default' created at location "
		            + getMem(*this));
		//##########################################################################
		_searchFlag = -1;
  }
  template<typename T>
  kdTree<T>::~kdTree()
  {
    //##########################################################################
		_log->TRACE("kdTree '" + _name
		            + "' destroyed at location " + getMem(*this));
		//##########################################################################
  }
  template<typename T>
  kdTree<T>::kdTree(std::shared_ptr<std::vector<std::vector<T>>> points)
  : _N(points->size()), _dim((*points)[0].size()), _name("default")
  {
    _points = points;
    //##########################################################################
		_log = std::make_shared<Log>();
		_log->init("ET:kdTree:default", ".logs/kdtree_default.txt");
		_log->TRACE("kdTree 'default' created at location "
		            + getMem(*this));
		//##########################################################################
  }

  //----------------------------------------------------------------------------
  //  Getters and Setters
  //----------------------------------------------------------------------------
  template<typename T>
  uint64_t kdTree<T>::getDim()
  {
    return _dim;
  }
  template<typename T>
  uint64_t kdTree<T>::getN()
  {
    return _N;
  }
  template<typename T>
  std::shared_ptr<std::vector<std::vector<T>>> kdTree<T>::getPoints()
  {
    return _points;
  }
  template<typename T>
  std::string kdTree<T>::getName()
  {
    return _name;
  }
  template<typename T>
  std::vector<std::vector<size_t>> kdTree<T>::getNeighbors()
  {
    return _neighbors;
  }
  template<typename T>
  std::vector<std::vector<double>> kdTree<T>::getDistances()
  {
    return _distances;
  }
  template<typename T>
  std::vector<std::vector<size_t>> kdTree<T>::getNeighborsRadius()
  {
    return _neighbors_radius;
  }
  template<typename T>
  std::vector<std::vector<double>> kdTree<T>::getDistancesRadius()
  {
    return _distances_radius;
  }
  template<typename T>
  std::vector<size_t> kdTree<T>::getNeighbors(uint64_t index)
  {
    if (index > _N)
    {
      //########################################################################
      _log->ERROR("kdTree " + _name
                  + ": Attempted to access neighbors array of size "
                  + std::to_string(_N) + " with index "
                  + std::to_string(index));
      //########################################################################
      if(_neighbors.size() > 0)
      {
        //######################################################################
        _log->INFO("kdTree " + _name + ": Returning the element at index 0");
        //######################################################################
        return _neighbors[0];
      }
      else
      {
        //######################################################################
        _log->INFO("kdTree " + _name + ": Returning empty neighbors array");
        //######################################################################
        return std::vector<size_t>(1,0);
      }
    }
    return _neighbors[index];
  }
  template<typename T>
  std::shared_ptr<Log> kdTree<T>::getLogger()
  {
    return _log;
  }
  //----------------------------------------------------------------------------
  //  KDTree methods
  //----------------------------------------------------------------------------
  template<typename T>
  void kdTree<T>::setupTree()
  {
    //  generate kdtree
    KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<T>>, T>
		kdt(_dim, *_points, 16);
    kdt.index->buildIndex();
		_kdtree = std::make_shared<KDTreeVectorOfVectorsAdaptor<
                     std::vector<std::vector<T>>,T>>(kdt);
    //##########################################################################
    _log->INFO("Initialized kdTree with " + std::to_string(_N) + " points "
               + "in " + std::to_string(_dim) + " dimensions");
    //##########################################################################
  }
  template<typename T>
  void kdTree<T>::queryNeighbors(uint64_t k)
  {
    //	check if anything has changed since last query
		if (_searchFlag == 0 && k == _k)
		{
			return;
		}
		if (k >= _N)
		{
			//########################################################################
			_log->ERROR("kdTree " + _name
									+ ": Attempted to query " + std::to_string(k)
									+ " neighbors for points in array _points of size "
									+ std::to_string(_N));
			//########################################################################
			return;
		}
    _neighbors.resize(_N);
		_distances.resize(_N);
		_k = k;
		//##########################################################################
		_log->INFO("kdTree " + _name + ": Querying each point in array _grid of"
	             + " size " + std::to_string(_N) + " and dimension "
						   + std::to_string(_dim) + " for the nearest "
							 + std::to_string(k) + " neighbors");
		//##########################################################################
    //  generate kdtree
    KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<T>>, T>
		kdt(_dim, *_points, 16);
    kdt.index->buildIndex();

    const size_t num_results = k;
    for (uint64_t i = 0; i < _N; i++)
    {
			_neighbors[i].resize(k);
			_distances[i].resize(k);
      std::vector<size_t> ret_indexes(num_results);
      std::vector<double> out_dists_sqr(num_results);

      nanoflann::KNNResultSet<double> resultSet(num_results);
      resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
			kdt.index->findNeighbors(resultSet, &(*_points)[i][0],
				                       nanoflann::SearchParams(10));
			_neighbors[i] = std::move(ret_indexes);
			_distances[i] = std::move(out_dists_sqr);
    }
		_searchFlag = 0;
  }
  template<typename T>
  std::vector<size_t>
  kdTree<T>::queryNeighbors(const std::vector<T>& point, uint64_t k)
  {
    if (k >= _N)
		{
			//########################################################################
			_log->ERROR("kdTree " + _name
									+ ": Attempted to query " + std::to_string(k)
									+ " neighbors for points in array _points of size "
									+ std::to_string(_N));
      _log->INFO("kdTree " + _name + ": Setting k to 3");
      k = 3;
			//########################################################################
		}
		//##########################################################################
		_log->INFO("kdTree " + _name + ": Querying a point in array _grid of"
	             + " size " + std::to_string(_N) + " and dimension "
						   + std::to_string(_dim) + " for the nearest "
							 + std::to_string(k) + " neighbors of a set of points");
		//##########################################################################

    const size_t num_results = k;
    std::vector<size_t> ret_indexes(num_results);
    std::vector<double> out_dists_sqr(num_results);
    nanoflann::KNNResultSet<double> resultSet(num_results);
    resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
    KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<T>>, T>
		kdt(_dim, *_points, 16);
    kdt.index->buildIndex();
		kdt.index->findNeighbors(resultSet, &point[0],
			                       nanoflann::SearchParams(10));
		return ret_indexes;
  }
  template<typename T>
  std::vector<double>
  kdTree<T>::queryDistances(const std::vector<T>& point, uint64_t k)
  {
    if (k >= _N)
		{
      //########################################################################
			_log->ERROR("kdTree " + _name
									+ ": Attempted to query " + std::to_string(k)
									+ " neighbors for points in array _points of size "
									+ std::to_string(_N));
      _log->INFO("kdTree " + _name + ": Setting k to 3");
      k = 3;
			//########################################################################
		}
		std::vector<size_t> neighbors(point.size());
		std::vector<double> distances(point.size());
		//##########################################################################
		_log->INFO("kdTree " + _name + ": Querying each point in array _grid of"
	             + " size " + std::to_string(_N) + " and dimension "
						   + std::to_string(_dim) + " for the nearest "
							 + std::to_string(k) + " neighbors of a set of points");
		//##########################################################################

    const size_t num_results = k;
		neighbors.resize(k);
		distances.resize(k);
    std::vector<size_t> ret_indexes(num_results);
    std::vector<double> out_dists_sqr(num_results);

    nanoflann::KNNResultSet<double> resultSet(num_results);
    resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
		_kdtree->index->findNeighbors(resultSet, &point[0],
			                       nanoflann::SearchParams(10));
		neighbors = std::move(ret_indexes);
		distances = std::move(out_dists_sqr);

		return distances;
  }
  template<typename T>
  std::vector<std::vector<size_t>>
  kdTree<T>::queryNeighbors(const std::vector<std::vector<T>>& points,
                            uint64_t k)
  {
    if (k >= _N)
		{
      //########################################################################
			_log->ERROR("kdTree " + _name
									+ ": Attempted to query " + std::to_string(k)
									+ " neighbors for points in array _points of size "
									+ std::to_string(_N));
      _log->INFO("kdTree " + _name + ": Setting k to 3");
      k = 3;
			//########################################################################
		}
		std::vector<std::vector<size_t>> neighbors(points.size());
		std::vector<std::vector<double>> distances(points.size());
		//##########################################################################
		_log->INFO("kdTree " + _name + ": Querying each point in array _grid of"
	             + " size " + std::to_string(_N) + " and dimension "
						   + std::to_string(_dim) + " for the nearest "
							 + std::to_string(k) + " neighbors of a set of points");
		//##########################################################################
    //  generate kdtree
    KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<T>>, T>
		kdt(_dim, *_points, 16);
    kdt.index->buildIndex();

    const size_t num_results = k;
    for (uint64_t i = 0; i < points.size(); i++)
    {
			neighbors[i].resize(k);
			distances[i].resize(k);
      std::vector<size_t> ret_indexes(num_results);
      std::vector<double> out_dists_sqr(num_results);

      nanoflann::KNNResultSet<double> resultSet(num_results);
      resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
			kdt.index->findNeighbors(resultSet, &points[i][0],
				                       nanoflann::SearchParams(10));
			neighbors[i] = std::move(ret_indexes);
			distances[i] = std::move(out_dists_sqr);
    }
		_kdtree = std::make_shared<KDTreeVectorOfVectorsAdaptor<
                     std::vector<std::vector<T>>,T>>(kdt);
		return neighbors;
  }
  template<typename T>
  void kdTree<T>::queryRadius(double radius)
  {
    //	check if anything has changed since last query
		if (_searchFlag == 0)
		{
			return;
		}
    _neighbors_radius.resize(_N);
		_distances_radius.resize(_N);
		_radius = radius;
		//	the algorithm looks for points that satisfy the squared
		//	distance rather than the square root.
		//##########################################################################
		_log->INFO("kdTree " + _name + ": Querying each point in array _grid of"
	             + " size " + std::to_string(_N) + " and dimension "
						   + std::to_string(_dim) + " for the nearest points within radius"
							 + std::to_string(radius));
		//##########################################################################
		radius *= radius;
		//  generate kdtree
    KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<T>>, T>
		kdt(_dim, *_points, 16);
    kdt.index->buildIndex();

    for (uint64_t i = 0; i < _N; i++)
    {
      std::vector<std::pair<size_t,double> > ret_matches;

			kdt.index->radiusSearch(&(*_points)[i][0], radius, ret_matches,
				                      nanoflann::SearchParams(10));
			std::vector<size_t> indices(ret_matches.size());
			std::vector<double> distances(ret_matches.size());
			for (uint64_t j = 0; j < ret_matches.size(); j++)
			{
				indices[j] = ret_matches[j].first;
				distances[j] = ret_matches[j].second;
			}
			_neighbors_radius[i] = std::move(indices);
			_distances_radius[i] = std::move(distances);
    }
		_searchFlag = 0;
		_kdtree = std::make_shared<KDTreeVectorOfVectorsAdaptor<
                     std::vector<std::vector<T>>,T>>(kdt);
  }
  //----------------------------------------------------------------------------
}
