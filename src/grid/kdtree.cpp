//------------------------------------------------------------------------------
//  kdtree.h
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

#include "kdtree.h"

namespace ET
{
  template<typename T>
  KDTree<T>::KDTree() : m_N(0), m_dim(0), m_name("default")
  {
    m_kdtree = std::make_shared<KDTreeVectorOfVectorsAdaptor<
                                std::vector<std::vector<T>>,T>>
                                (m_dim, *m_points, 16);
    m_log = std::make_shared<Log>();
		m_log->init("ET:KDTree:" + m_name, ".logs/kdtree_" + m_name + ".txt");
		m_log->TRACE("KDTree '" + m_name + "' created at location "
		            + getMem(*this));
  }
  template<typename T>
  KDTree<T>::~KDTree()
  {
		m_log->TRACE("KDTree '" + m_name
		            + "' destroyed at location " + getMem(*this));
  }
  template<typename T>
  KDTree<T>::KDTree(std::shared_ptr<std::vector<std::vector<T>>> t_points)
  : m_N(t_points->size()), m_dim((*t_points)[0].size()), m_name("default")
  {
    m_points = t_points;
    m_kdtree = std::make_shared<KDTreeVectorOfVectorsAdaptor<
                                std::vector<std::vector<T>>,T>>
                                (m_dim, *m_points, 16);
    m_log = std::make_shared<Log>();
		m_log->init("ET:KDTree:" + m_name, ".logs/kdtree_" + m_name + ".txt");
		m_log->TRACE("KDTree '" + m_name + "' created at location "
		            + getMem(*this));
  }

  //----------------------------------------------------------------------------
  //  Getters and Setters
  //----------------------------------------------------------------------------
  template<typename T>
  std::string KDTree<T>::getName() const
  {
    return m_name;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  size_t KDTree<T>::getDim() const
  {
    return m_dim;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  size_t KDTree<T>::getN() const
  {
    return m_N;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::shared_ptr<std::vector<std::vector<T>>> KDTree<T>::getPoints() const
  {
    return m_points;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<std::vector<size_t>> KDTree<T>::getCurrentNeighborIndices()
  {
    return m_currentNeighborIndices;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<std::vector<double>> KDTree<T>::getCurrentNeighborDistances()
  {
    return m_currentNeighborDistances;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<size_t> KDTree<T>::getCurrentNeighborIndices(size_t t_i)
  {
    //  check that index exists
    if (t_i >= m_N) {
      m_log->ERROR("Attempted to access neighbors for index "
                   + std::to_string(t_i) + " for tree with "
                   + std::to_string(m_N) + " points.");
      m_log->INFO("Returning neighbor iindices for first point.");
      return m_currentNeighborIndices[0];
    }
    //  otherwise return the point
    return m_currentNeighborIndices[t_i];
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<double> KDTree<T>::getCurrentNeighborDistances(size_t t_i)
  {
    //  check that index exists
    if (t_i >= m_N) {
      m_log->ERROR("Attempted to access neighbor distances for index "
                   + std::to_string(t_i) + " for tree with "
                   + std::to_string(m_N) + " points.");
      m_log->INFO("Returning neighbor distances for first point.");
      return m_currentNeighborDistances[0];
    }
    //  otherwise return the point
    return m_currentNeighborDistances[t_i];
  }
  //----------------------------------------------------------------------------
  template<typename T>
  size_t KDTree<T>::getCurrentGlobalK() const
  {
    return m_currentGlobalK;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  double KDTree<T>::getCurrentGlobalRadius() const
  {
    return m_currentGlobalRadius;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::shared_ptr<Log> KDTree<T>::getLog() const
  {
    return m_log;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void KDTree<T>::setName(const std::string t_name)
  {
    m_name = t_name;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void KDTree<T>::setDim(const size_t t_dim)
  {
    m_dim = t_dim;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void KDTree<T>::setN(const size_t t_N)
  {
    m_N = t_N;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void
  KDTree<T>::setPoints(const std::shared_ptr<std::vector<std::vector<T>>> t_points)
  {
    m_points = t_points;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void KDTree<T>::setLog(std::shared_ptr<Log> t_log)
  {
    m_log = t_log;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void KDTree<T>::setCurrentGlobalK(const size_t t_k)
  {
    m_currentGlobalK = t_k;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void KDTree<T>::setCurrentGlobalRadius(const double t_radius)
  {
    m_currentGlobalRadius = t_radius;
  }
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //  KDTree methods
  //----------------------------------------------------------------------------
  template<typename T>
  void KDTree<T>::setupTree()
  {
    m_kdtree = std::make_shared<KDTreeVectorOfVectorsAdaptor<
                                std::vector<std::vector<T>>,T>>
                                (m_dim, *m_points, 16);
    m_kdtree->index->buildIndex();
    m_log->INFO("Initialized KDTree with " + std::to_string(m_N) + " points "
               + "in " + std::to_string(m_dim) + " dimensions");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void KDTree<T>::queryNeighbors(size_t t_k)
  {
    //	check if anything has changed since last query
		if (m_searchFlag == 0 && t_k == m_currentGlobalK) {
			return;
		}
		if (t_k >= m_N) {
			m_log->ERROR("KDTree " + m_name
									+ ": Attempted to query " + std::to_string(m_currentGlobalK)
									+ " neighbors for points in array m_points of size "
									+ std::to_string(m_N));
			return;
		}
    m_currentNeighborIndices.resize(m_N);
		m_currentNeighborDistances.resize(m_N);
		m_currentGlobalK = t_k;
		m_log->INFO("KDTree " + m_name + ": Querying each point in array _grid of"
	             + " size " + std::to_string(m_N) + " and dimension "
						   + std::to_string(m_dim) + " for the nearest "
							 + std::to_string(m_currentGlobalK) + " neighbors");
    //  generate KDTree
    // KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<T>>, T>
		// kdt(m_dim, *m_points, 16);
    // m_kdtree->index->buildIndex();

    const size_t num_results = m_currentGlobalK;
    for (auto i = 0; i < m_N; i++) {
			m_currentNeighborIndices[i].resize(m_currentGlobalK);
			m_currentNeighborDistances[i].resize(m_currentGlobalK);
      std::vector<size_t> ret_indexes(num_results);
      std::vector<double> out_dists_sqr(num_results);

      nanoflann::KNNResultSet<double> resultSet(num_results);
      resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
			m_kdtree->index->findNeighbors(resultSet, &(*m_points)[i][0],
				                       nanoflann::SearchParams(10));
			m_currentNeighborIndices[i] = std::move(ret_indexes);
			m_currentNeighborDistances[i] = std::move(out_dists_sqr);
    }
		m_searchFlag = 0;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<size_t>
  KDTree<T>::queryNeighbors(const std::vector<T>& t_point, size_t t_k)
  {
    if (t_k >= m_N) {
			m_log->ERROR("KDTree " + m_name
									+ ": Attempted to query " + std::to_string(t_k)
									+ " neighbors for points in array m_points of size "
									+ std::to_string(m_N));
      m_log->INFO("KDTree " + m_name + ": Setting k to 3");
      m_currentGlobalK = 3;
		}
    else {
      m_currentGlobalK = t_k;
    }
		m_log->INFO("KDTree " + m_name + ": Querying a point in array _grid of"
	             + " size " + std::to_string(m_N) + " and dimension "
						   + std::to_string(m_dim) + " for the nearest "
							 + std::to_string(m_currentGlobalK) + " neighbors of a set of points");

    const size_t num_results = m_currentGlobalK;
    std::vector<size_t> ret_indexes(num_results);
    std::vector<double> out_dists_sqr(num_results);
    nanoflann::KNNResultSet<double> resultSet(num_results);
    resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);

    // KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<T>>, T>
		// kdt(m_dim, *m_points, 16);
    // m_kdtree->index->buildIndex();
		m_kdtree->index->findNeighbors(resultSet, &t_point[0],
			                       nanoflann::SearchParams(10));
		return ret_indexes;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<double>
  KDTree<T>::queryDistances(const std::vector<T>& t_point, size_t t_k)
  {
    if (t_k >= m_N) {
			m_log->ERROR("KDTree " + m_name
									+ ": Attempted to query " + std::to_string(t_k)
									+ " neighbors for points in array m_points of size "
									+ std::to_string(m_N));
      m_log->INFO("KDTree " + m_name + ": Setting k to 3");
      m_currentGlobalK = 3;
		}
    else {
      m_currentGlobalK = t_k;
    }
		std::vector<size_t> neighbors(t_point.size());
		std::vector<double> distances(t_point.size());
		m_log->INFO("KDTree " + m_name + ": Querying each point in array _grid of"
	             + " size " + std::to_string(m_N) + " and dimension "
						   + std::to_string(m_dim) + " for the nearest "
							 + std::to_string(m_currentGlobalK) + " neighbors of a set of points");

    const size_t num_results = m_currentGlobalK;
		neighbors.resize(m_currentGlobalK);
		distances.resize(m_currentGlobalK);
    std::vector<size_t> ret_indexes(num_results);
    std::vector<double> out_dists_sqr(num_results);

    nanoflann::KNNResultSet<double> resultSet(num_results);
    resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
    // KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<T>>, T>
		// kdt(m_dim, *m_points, 16);
    // m_kdtree->index->buildIndex();
		m_kdtree->index->findNeighbors(resultSet, &t_point[0],
			                       nanoflann::SearchParams(10));
		neighbors = std::move(ret_indexes);
		distances = std::move(out_dists_sqr);

		return distances;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<std::vector<size_t>>
  KDTree<T>::queryNeighbors(const std::vector<std::vector<T>>& t_points,
                            size_t t_k)
  {
    if (t_k >= m_N) {
			m_log->ERROR("KDTree " + m_name
									+ ": Attempted to query " + std::to_string(t_k)
									+ " neighbors for points in array m_points of size "
									+ std::to_string(m_N));
      m_log->INFO("KDTree " + m_name + ": Setting k to 3");
      m_currentGlobalK = 3;
		}
    else {
      m_currentGlobalK = t_k;
    }
		std::vector<std::vector<size_t>> neighbors(t_points.size());
		std::vector<std::vector<double>> distances(t_points.size());
		m_log->INFO("KDTree " + m_name + ": Querying each point in array _grid of"
	             + " size " + std::to_string(m_N) + " and dimension "
						   + std::to_string(m_dim) + " for the nearest "
							 + std::to_string(m_currentGlobalK) + " neighbors of a set of points");
    //  generate KDTree
    // KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<T>>, T>
		// kdt(m_dim, *m_points, 16);
    // m_kdtree->index->buildIndex();

    const size_t num_results = m_currentGlobalK;
    for (auto i = 0; i < t_points.size(); i++) {
			neighbors[i].resize(m_currentGlobalK);
			distances[i].resize(m_currentGlobalK);
      std::vector<size_t> ret_indexes(num_results);
      std::vector<double> out_dists_sqr(num_results);

      nanoflann::KNNResultSet<double> resultSet(num_results);
      resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
			m_kdtree->index->findNeighbors(resultSet, &t_points[i][0],
				                       nanoflann::SearchParams(10));
			neighbors[i] = std::move(ret_indexes);
			distances[i] = std::move(out_dists_sqr);
    }
		return neighbors;
  }
  // //----------------------------------------------------------------------------
  // template<typename T>
  // void KDTree<T>::queryRadius(double t_radius)
  // {
  //   //	check if anything has changed since last query
	// 	if (m_searchFlag == 0) {
	// 		return;
	// 	}
  //   m_currentNeighborIndices.resize(m_N);
	// 	m_currentNeighborDistances.resize(m_N);
	// 	m_currentGlobalRadius = t_radius;
	// 	//	the algorithm looks for points that satisfy the squared
	// 	//	distance rather than the square root.
	// 	m_log->INFO("KDTree " + m_name + ": Querying each point in array _grid of"
	//              + " size " + std::to_string(m_N) + " and dimension "
	// 					   + std::to_string(m_dim) + " for the nearest points within radius"
	// 						 + std::to_string(m_currentGlobalRadius));
  //   //  Nanoflann uses r^2 rather than r, since it requires one fewer
  //   //  operations to calculate.
	// 	t_radius *= t_radius;
	// 	//  generate KDTree
  //   KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<T>>, T>
	// 	kdt(m_dim, *m_points, 16);
  //   m_kdtree->index->buildIndex();
  //
  //   for (auto i = 0; i < m_N; i++) {
  //     std::vector<std::pair<size_t,double>> ret_matches;
  //
	// 		m_kdtree->index->radiusSearch(&(*m_points)[i][0], t_radius, ret_matches,
	// 			                      nanoflann::SearchParams(10));
	// 		std::vector<size_t> indices(ret_matches.size());
	// 		std::vector<double> distances(ret_matches.size());
	// 		for (auto j = 0; j < ret_matches.size(); j++) {
	// 			indices[j] = ret_matches[j].first;
	// 			distances[j] = ret_matches[j].second;
	// 		}
	// 		m_currentNeighborIndices[i] = std::move(indices);
	// 		m_currentNeighborDistances[i] = std::move(distances);
  //   }
	// 	m_searchFlag = 0;
  // }
  // //----------------------------------------------------------------------------
}
