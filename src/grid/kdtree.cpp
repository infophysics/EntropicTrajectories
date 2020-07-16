//------------------------------------------------------------------------------
//  KDTree.h
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

#include "KDTree.h"

namespace ET
{
  template<typename T>
  KDTree<T>::KDTree() : m_N(0), m_dim(0), m_name("default")
  {
		m_log = std::make_shared<Log>();
		m_log->init("ET:KDTree:default", ".logs/KDTree_default.txt");
		m_log->TRACE("KDTree 'default' created at location "
		            + getMem(*this));
		m_searchFlag = -1;
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
		m_log = std::make_shared<Log>();
		m_log->init("ET:KDTree:default", ".logs/KDTree_default.txt");
		m_log->TRACE("KDTree 'default' created at location "
		            + getMem(*this));
  }

  //----------------------------------------------------------------------------
  //  Getters and Setters
  //----------------------------------------------------------------------------
  template<typename T>
  uint64_t KDTree<T>::getDim()
  {
    return m_dim;
  }
  template<typename T>
  uint64_t KDTree<T>::getN()
  {
    return m_N;
  }
  template<typename T>
  std::shared_ptr<std::vector<std::vector<T>>> KDTree<T>::getPoints()
  {
    return m_points;
  }
  template<typename T>
  std::string KDTree<T>::getName()
  {
    return m_name;
  }
  template<typename T>
  std::vector<std::vector<size_t>> KDTree<T>::getNeighbors()
  {
    return m_neighbors;
  }
  template<typename T>
  std::vector<std::vector<double>> KDTree<T>::getDistances()
  {
    return m_distances;
  }
  template<typename T>
  std::vector<std::vector<size_t>> KDTree<T>::getNeighborsRadius()
  {
    return m_neighbors_radius;
  }
  template<typename T>
  std::vector<std::vector<double>> KDTree<T>::getDistancesRadius()
  {
    return m_distances_radius;
  }
  template<typename T>
  std::vector<size_t> KDTree<T>::getNeighbors(uint64_t t_index)
  {
    if (t_index > m_N) {
      m_log->ERROR("KDTree " + m_name
                  + ": Attempted to access neighbors array of size "
                  + std::to_string(m_N) + " with index "
                  + std::to_string(t_index));
      if(m_neighbors.size() > 0) {
        m_log->INFO("KDTree " + m_name + ": Returning the element at index 0");
        return m_neighbors[0];
      }
      else {
        m_log->INFO("KDTree " + m_name + ": Returning empty neighbors array");
        return std::vector<size_t>(1,0);
      }
    }
    return m_neighbors[t_index];
  }
  template<typename T>
  std::shared_ptr<Log> KDTree<T>::getLogger()
  {
    return m_log;
  }
  //----------------------------------------------------------------------------
  //  KDTree methods
  //----------------------------------------------------------------------------
  template<typename T>
  void KDTree<T>::setupTree()
  {
    //  generate KDTree
    KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<T>>, T>
		kdt(m_dim, *m_points, 16);
    kdt.index->buildIndex();
		m_KDTree = std::make_shared<KDTreeVectorOfVectorsAdaptor<
                     std::vector<std::vector<T>>,T>>(kdt);
    m_log->INFO("Initialized KDTree with " + std::to_string(m_N) + " points "
               + "in " + std::to_string(m_dim) + " dimensions");
  }
  template<typename T>
  void KDTree<T>::queryNeighbors(uint64_t t_k)
  {
    //	check if anything has changed since last query
		if (m_searchFlag == 0 && t_k == m_k) {
			return;
		}
		if (t_k >= m_N) {
			m_log->ERROR("KDTree " + m_name
									+ ": Attempted to query " + std::to_string(m_k)
									+ " neighbors for points in array m_points of size "
									+ std::to_string(m_N));
			return;
		}
    m_neighbors.resize(m_N);
		m_distances.resize(m_N);
		m_k = t_k;
		m_log->INFO("KDTree " + m_name + ": Querying each point in array _grid of"
	             + " size " + std::to_string(m_N) + " and dimension "
						   + std::to_string(m_dim) + " for the nearest "
							 + std::to_string(m_k) + " neighbors");
    //  generate KDTree
    KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<T>>, T>
		kdt(m_dim, *m_points, 16);
    kdt.index->buildIndex();

    const size_t num_results = m_k;
    for (auto i = 0; i < m_N; i++) {
			m_neighbors[i].resize(m_k);
			m_distances[i].resize(m_k);
      std::vector<size_t> ret_indexes(num_results);
      std::vector<double> out_dists_sqr(num_results);

      nanoflann::KNNResultSet<double> resultSet(num_results);
      resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
			kdt.index->findNeighbors(resultSet, &(*m_points)[i][0],
				                       nanoflann::SearchParams(10));
			m_neighbors[i] = std::move(ret_indexes);
			m_distances[i] = std::move(out_dists_sqr);
    }
    // m_KDTree = std::make_shared<KDTreeVectorOfVectorsAdaptor<
    //                            std::vector<std::vector<T>>, T>>(kdt);
		m_searchFlag = 0;
  }
  template<typename T>
  std::vector<size_t>
  KDTree<T>::queryNeighbors(const std::vector<T>& t_point, uint64_t t_k)
  {
    if (t_k >= m_N) {
			m_log->ERROR("KDTree " + m_name
									+ ": Attempted to query " + std::to_string(t_k)
									+ " neighbors for points in array m_points of size "
									+ std::to_string(m_N));
      m_log->INFO("KDTree " + m_name + ": Setting k to 3");
      m_k = 3;
		}
    else {
      m_k = t_k;
    }
		m_log->INFO("KDTree " + m_name + ": Querying a point in array _grid of"
	             + " size " + std::to_string(m_N) + " and dimension "
						   + std::to_string(m_dim) + " for the nearest "
							 + std::to_string(m_k) + " neighbors of a set of points");

    const size_t num_results = m_k;
    std::vector<size_t> ret_indexes(num_results);
    std::vector<double> out_dists_sqr(num_results);
    nanoflann::KNNResultSet<double> resultSet(num_results);
    resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);

    KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<T>>, T>
		kdt(m_dim, *m_points, 16);
    kdt.index->buildIndex();
		kdt.index->findNeighbors(resultSet, &t_point[0],
			                       nanoflann::SearchParams(10));
    // m_KDTree = std::make_shared<KDTreeVectorOfVectorsAdaptor<
    //                           std::vector<std::vector<T>>, T>>(kdt);
		return ret_indexes;
  }
  template<typename T>
  std::vector<double>
  KDTree<T>::queryDistances(const std::vector<T>& t_point, uint64_t t_k)
  {
    if (t_k >= m_N) {
			m_log->ERROR("KDTree " + m_name
									+ ": Attempted to query " + std::to_string(t_k)
									+ " neighbors for points in array m_points of size "
									+ std::to_string(m_N));
      m_log->INFO("KDTree " + m_name + ": Setting k to 3");
      m_k = 3;
		}
    else {
      m_k = t_k;
    }
		std::vector<size_t> neighbors(t_point.size());
		std::vector<double> distances(t_point.size());
		m_log->INFO("KDTree " + m_name + ": Querying each point in array _grid of"
	             + " size " + std::to_string(m_N) + " and dimension "
						   + std::to_string(m_dim) + " for the nearest "
							 + std::to_string(m_k) + " neighbors of a set of points");

    const size_t num_results = m_k;
		neighbors.resize(m_k);
		distances.resize(m_k);
    std::vector<size_t> ret_indexes(num_results);
    std::vector<double> out_dists_sqr(num_results);

    nanoflann::KNNResultSet<double> resultSet(num_results);
    resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
    KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<T>>, T>
		kdt(m_dim, *m_points, 16);
    kdt.index->buildIndex();
		kdt.index->findNeighbors(resultSet, &t_point[0],
			                       nanoflann::SearchParams(10));
		neighbors = std::move(ret_indexes);
		distances = std::move(out_dists_sqr);

		return distances;
  }
  template<typename T>
  std::vector<std::vector<size_t>>
  KDTree<T>::queryNeighbors(const std::vector<std::vector<T>>& t_points,
                            uint64_t t_k)
  {
    if (t_k >= m_N) {
			m_log->ERROR("KDTree " + m_name
									+ ": Attempted to query " + std::to_string(t_k)
									+ " neighbors for points in array m_points of size "
									+ std::to_string(m_N));
      m_log->INFO("KDTree " + m_name + ": Setting k to 3");
      m_k = 3;
		}
    else {
      m_k = t_k;
    }
		std::vector<std::vector<size_t>> neighbors(t_points.size());
		std::vector<std::vector<double>> distances(t_points.size());
		m_log->INFO("KDTree " + m_name + ": Querying each point in array _grid of"
	             + " size " + std::to_string(m_N) + " and dimension "
						   + std::to_string(m_dim) + " for the nearest "
							 + std::to_string(m_k) + " neighbors of a set of points");
    //  generate KDTree
    KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<T>>, T>
		kdt(m_dim, *m_points, 16);
    kdt.index->buildIndex();

    const size_t num_results = m_k;
    for (auto i = 0; i < t_points.size(); i++) {
			neighbors[i].resize(m_k);
			distances[i].resize(m_k);
      std::vector<size_t> ret_indexes(num_results);
      std::vector<double> out_dists_sqr(num_results);

      nanoflann::KNNResultSet<double> resultSet(num_results);
      resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
			kdt.index->findNeighbors(resultSet, &t_points[i][0],
				                       nanoflann::SearchParams(10));
			neighbors[i] = std::move(ret_indexes);
			distances[i] = std::move(out_dists_sqr);
    }
		// m_KDTree = std::make_shared<KDTreeVectorOfVectorsAdaptor<
    //                  std::vector<std::vector<T>>,T>>(kdt);
		return neighbors;
  }
  template<typename T>
  void KDTree<T>::queryRadius(double t_radius)
  {
    //	check if anything has changed since last query
		if (m_searchFlag == 0) {
			return;
		}
    m_neighbors_radius.resize(m_N);
		m_distances_radius.resize(m_N);
		m_radius = t_radius;
		//	the algorithm looks for points that satisfy the squared
		//	distance rather than the square root.
		m_log->INFO("KDTree " + m_name + ": Querying each point in array _grid of"
	             + " size " + std::to_string(m_N) + " and dimension "
						   + std::to_string(m_dim) + " for the nearest points within radius"
							 + std::to_string(m_radius));
    //  Nanoflann uses r^2 rather than r, since it requires one fewer
    //  operations to calculate.
		t_radius *= t_radius;
		//  generate KDTree
    KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<T>>, T>
		kdt(m_dim, *m_points, 16);
    kdt.index->buildIndex();

    for (auto i = 0; i < m_N; i++) {
      std::vector<std::pair<size_t,double> > ret_matches;

			kdt.index->radiusSearch(&(*m_points)[i][0], t_radius, ret_matches,
				                      nanoflann::SearchParams(10));
			std::vector<size_t> indices(ret_matches.size());
			std::vector<double> distances(ret_matches.size());
			for (auto j = 0; j < ret_matches.size(); j++) {
				indices[j] = ret_matches[j].first;
				distances[j] = ret_matches[j].second;
			}
			m_neighbors_radius[i] = std::move(indices);
			m_distances_radius[i] = std::move(distances);
    }
		m_searchFlag = 0;
		// m_KDTree = std::make_shared<KDTreeVectorOfVectorsAdaptor<
    //                  std::vector<std::vector<T>>,T>>(kdt);
  }
  //----------------------------------------------------------------------------
}
