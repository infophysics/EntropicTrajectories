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
  //  map for a string to enum of KDTreeBackend
  std::map<std::string, KDTreeBackend> KDTreeBackendMap =
  {
    { "NANOFLANN",  KDTreeBackend::NANOFLANN },
    { "FLANN",      KDTreeBackend::FLANN }
  };
  //  map for a enum of KDTreeBackend to a string
  std::map<KDTreeBackend, std::string> KDTreeBackendNameMap =
  {
    { KDTreeBackend::NANOFLANN, "NANOFLANN" },
    { KDTreeBackend::FLANN,     "FLANN" }
  };

  template<typename T>
  KDTree<T>::KDTree() : m_N(0), m_dim(0), m_name("default")
  {
    m_log = std::make_shared<Log>();
		m_log->init("ET:KDTree:" + m_name, ".logs/kdtree_" + m_name + ".txt");
		m_log->TRACE("KDTree '" + m_name + "' created at location "
		            + address_to_string(*this));
  }
  template<typename T>
  KDTree<T>::~KDTree()
  {
		m_log->TRACE("KDTree '" + m_name
		            + "' destroyed at location " + address_to_string(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  KDTree<T>::KDTree(const std::shared_ptr<std::vector<std::vector<T>>> t_points)
  : m_N(t_points->size()), m_dim((*t_points)[0].size()), m_name("default"),
    m_points(t_points)
  {
    m_log = std::make_shared<Log>();
		m_log->init("ET:KDTree:" + m_name, ".logs/kdtree_" + m_name + ".txt");
		m_log->TRACE("KDTree '" + m_name + "' created at location "
		            + address_to_string(*this));
    setupTree();
    m_currentNeighborIndices.resize(m_N);
    m_currentNeighborDistances.resize(m_N);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  KDTree<T>::KDTree(const std::string t_name,
                    const std::shared_ptr<std::vector<std::vector<T>>> t_points)
  : m_N(t_points->size()), m_dim((*t_points)[0].size()), m_name(t_name),
    m_points(t_points)
  {
    m_log = std::make_shared<Log>();
    m_log->init("ET:KDTree:" + m_name, ".logs/kdtree_" + m_name + ".txt");
    m_log->TRACE("KDTree '" + m_name + "' created at location "
                + address_to_string(*this));
    setupTree();
    m_currentNeighborIndices.resize(m_N);
    m_currentNeighborDistances.resize(m_N);
  }
  //----------------------------------------------------------------------------
  //  Getters and Setters
  //----------------------------------------------------------------------------
  template<typename T>
  enum KDTreeBackend KDTree<T>::getBackend() const
  {
    return m_backend;
  }
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
  enum KDTreeSearchFlags KDTree<T>::getSearchFlag() const
  {
    return m_searchFlag;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void KDTree<T>::setBackend(enum KDTreeBackend t_backend)
  {
    m_backend = t_backend;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void KDTree<T>::setName(const std::string t_name)
  {
    m_name = t_name;
    m_log->INFO("Set name to '" + m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void KDTree<T>::setDim(const size_t t_dim)
  {
    m_dim = t_dim;
    m_log->INFO("Set dimension to " + std::to_string(m_dim));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void KDTree<T>::setN(const size_t t_N)
  {
    m_N = t_N;
    m_currentNeighborIndices.resize(m_N);
    m_currentNeighborDistances.resize(m_N);
    m_log->INFO("Set number of points to " + std::to_string(m_N));
    m_searchFlag = KDTreeSearchFlags::OUTDATED;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void
  KDTree<T>::setPoints(const std::shared_ptr<std::vector<std::vector<T>>> t_points)
  {
    m_points = t_points;
    setN(m_points->size());
    setDim((*m_points)[0].size());
    m_currentNeighborIndices.resize(m_N);
    m_currentNeighborDistances.resize(m_N);
    m_log->INFO("Set points to array of size (" + std::to_string(m_N)
                + " x " + std::to_string(m_dim) + ")");
    setupTree();
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void KDTree<T>::setLog(std::shared_ptr<Log> t_log)
  {
    m_log = t_log;
    m_log->INFO("Logger passed to KDTree '" + m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void KDTree<T>::setCurrentGlobalK(const size_t t_k)
  {
    m_currentGlobalK = t_k;
    m_searchFlag = KDTreeSearchFlags::NEW_K;
    m_log->INFO("Set global k value to " + std::to_string(t_k));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void KDTree<T>::setCurrentGlobalRadius(const double t_radius)
  {
    m_currentGlobalRadius = t_radius;
    m_searchFlag = KDTreeSearchFlags::NEW_RADIUS;
    m_log->INFO("Set global radius value to " + std::to_string(t_radius));
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
    m_searchFlag = KDTreeSearchFlags::BUILT;
    m_log->INFO("Initialized KDTree with " + std::to_string(m_N) + " points "
               + "in " + std::to_string(m_dim) + " dimensions");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void KDTree<T>::queryNeighbors(size_t t_k)
  {
    //	check if anything has changed since last query
    if (t_k >= m_N) {
			m_log->ERROR("KDTree " + m_name
									+ ": Attempted to query " + std::to_string(m_currentGlobalK)
									+ " neighbors for points in array m_points of size "
									+ std::to_string(m_N) + ".  Aborting search.");
			return;
		}
    if (t_k != m_currentGlobalK) {
      setCurrentGlobalK(t_k);
    }
		if (m_searchFlag == KDTreeSearchFlags::CURRENT_K) {
      m_log->INFO("KDTree is up to date with " + std::to_string(t_k)
                  + " neighbor values.  Aborting search.");
			return;
		}
		m_log->INFO("KDTree " + m_name + ": Querying each point in array m_points of"
	             + " size (" + std::to_string(m_N) + " x "
						   + std::to_string(m_dim) + ") for the nearest "
							 + std::to_string(m_currentGlobalK) + " neighbors");
    const size_t num_results = m_currentGlobalK;
    for (auto i = 0; i < m_N; i++) {
      std::vector<size_t> ret_indexes(num_results);
      std::vector<double> out_dists_sqr(num_results);

      nanoflann::KNNResultSet<double> resultSet(num_results);
      resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
			m_kdtree->index->findNeighbors(resultSet, &(*m_points)[i][0],
				                       nanoflann::SearchParams(10));
			m_currentNeighborIndices[i] = std::move(ret_indexes);
			m_currentNeighborDistances[i] = std::move(out_dists_sqr);
    }
		m_searchFlag = KDTreeSearchFlags::CURRENT_K;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void KDTree<T>::queryNeighbors(double t_radius)
  {
    //  make sure that radius is positive
    if (t_radius < 0) {
      m_log->WARN("Attempted to use a negative search radius "
                  + std::to_string(t_radius) + " in nearest neighbor search.");
      m_log->INFO("Setting radius to its absolute value.");
      t_radius = abs(t_radius);
    }
    //	check if anything has changed since last query
    if (t_radius != m_currentGlobalRadius) {
      setCurrentGlobalRadius(t_radius);
    }
		if (m_searchFlag == KDTreeSearchFlags::CURRENT_RADIUS) {
      m_log->INFO("KDTree is up to date with " + std::to_string(t_radius)
                  + " search radius.  Aborting search.");
			return;
		}
		m_log->INFO("KDTree " + m_name + ": Querying each point in array m_points of"
	             + " size (" + std::to_string(m_N) + " x "
						   + std::to_string(m_dim) + ") for neighbors within radius "
							 + std::to_string(m_currentGlobalRadius));
    //  WARNING: Nanoflann uses radius^2 for searchs since it is
    //  computationally simpler than taking the square root for every
    //  iteration.
    t_radius *= t_radius;
    for (auto i = 0; i < m_N; i++) {
      std::vector<std::pair<size_t,double>> ret_matches;

			m_kdtree->index->radiusSearch(&(*m_points)[i][0], t_radius, ret_matches,
				                      nanoflann::SearchParams(10));

			std::vector<size_t> indices(ret_matches.size());
			std::vector<double> distances(ret_matches.size());
			for (auto j = 0; j < ret_matches.size(); j++) {
				indices[j] = ret_matches[j].first;
				distances[j] = ret_matches[j].second;
			}
			m_currentNeighborIndices[i] = std::move(indices);
			m_currentNeighborDistances[i] = std::move(distances);
    }
		m_searchFlag = KDTreeSearchFlags::CURRENT_RADIUS;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<size_t>
  KDTree<T>::queryNeighbors(const std::vector<T>& t_point, size_t t_k)
  {
    //  check that t_point has the same dimension as m_points.
    if (t_point.size() != m_dim) {
      m_log->ERROR("Attempted to search a tree of dimension "
                   + std::to_string(m_dim) + " with a point of dimension "
                   + std::to_string(t_point.size()));
      m_log->INFO("Returning empty vector.");
      return std::vector<size_t>(0);
    }
    if (t_k >= m_N) {
			m_log->ERROR("KDTree " + m_name
									+ ": Attempted to query " + std::to_string(t_k)
									+ " neighbors for points in array m_points of size "
									+ std::to_string(m_N));
      m_log->INFO("Setting k to 3 for this search.");
      t_k = 3;
		}
		m_log->INFO("KDTree " + m_name + ": Querying a point for the nearest "
							  + std::to_string(t_k) + " neighbors.");

    const size_t num_results = m_currentGlobalK;
    std::vector<size_t> ret_indexes(num_results);
    std::vector<double> out_dists_sqr(num_results);
    nanoflann::KNNResultSet<double> resultSet(num_results);
    resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
		m_kdtree->index->findNeighbors(resultSet, &t_point[0],
			                       nanoflann::SearchParams(10));
		return ret_indexes;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<double>
  KDTree<T>::queryDistances(const std::vector<T>& t_point, size_t t_k)
  {
    //  check that t_point has the same dimension as m_points.
    if (t_point.size() != m_dim) {
      m_log->ERROR("Attempted to search a tree of dimension "
                   + std::to_string(m_dim) + " with a point of dimension "
                   + std::to_string(t_point.size()));
      m_log->INFO("Returning empty vector.");
      return std::vector<double>(0);
    }
    if (t_k >= m_N) {
			m_log->ERROR("KDTree " + m_name
									+ ": Attempted to query " + std::to_string(t_k)
									+ " neighbors for points in array m_points of size "
									+ std::to_string(m_N));
      m_log->INFO("Setting k to 3 for this search.");
      t_k = 3;
		}
		m_log->INFO("KDTree " + m_name + ": Querying a point for the nearest "
							  + std::to_string(t_k) + " neighbors.");
    const size_t num_results = m_currentGlobalK;
    std::vector<size_t> ret_indexes(num_results);
    std::vector<double> out_dists_sqr(num_results);
    nanoflann::KNNResultSet<double> resultSet(num_results);
    resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
    m_kdtree->index->findNeighbors(resultSet, &t_point[0],
    	                       nanoflann::SearchParams(10));
    // convert distances squared to distances
    std::vector<double> out_dists(out_dists_sqr.size());
    for (auto i = 0; i < out_dists_sqr.size(); i++) {
      out_dists[i] = sqrt(out_dists_sqr[i]);
    }
		return out_dists;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::tuple<std::vector<size_t>,std::vector<double>>
  KDTree<T>::query(const std::vector<T>& t_point, size_t t_k)
  {
    //  check that t_point has the same dimension as m_points.
    if (t_point.size() != m_dim) {
      m_log->ERROR("Attempted to search a tree of dimension "
                   + std::to_string(m_dim) + " with a point of dimension "
                   + std::to_string(t_point.size()));
      m_log->INFO("Returning empty vector.");
      return std::tuple<std::vector<size_t>,std::vector<double>>();
    }
    if (t_k >= m_N) {
      m_log->ERROR("KDTree " + m_name
                  + ": Attempted to query " + std::to_string(t_k)
                  + " neighbors for points in array m_points of size "
                  + std::to_string(m_N));
      m_log->INFO("Setting k to 3 for this search.");
      t_k = 3;
    }
    m_log->INFO("KDTree " + m_name + ": Querying a point for the nearest "
                + std::to_string(t_k) + " neighbors.");

    const size_t num_results = m_currentGlobalK;
    std::vector<size_t> ret_indexes(num_results);
    std::vector<double> out_dists_sqr(num_results);
    nanoflann::KNNResultSet<double> resultSet(num_results);
    resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
    m_kdtree->index->findNeighbors(resultSet, &t_point[0],
                             nanoflann::SearchParams(10));
    std::vector<double> out_dists(out_dists_sqr.size());
    for (auto i = 0; i < out_dists_sqr.size(); i++) {
     out_dists[i] = sqrt(out_dists_sqr[i]);
    }
    std::tuple<std::vector<size_t>,std::vector<double>>
    result {ret_indexes, out_dists};
    return result;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<std::vector<size_t>>
  KDTree<T>::queryNeighbors(const std::vector<std::vector<T>>& t_points,
                            size_t t_k)
  {
    //  check that t_point has the same dimension as m_points.
    if (t_points[0].size() != m_dim) {
      m_log->ERROR("Attempted to search a tree of dimension "
                   + std::to_string(m_dim) + " with points of dimension "
                   + std::to_string(t_points[0].size()));
      m_log->INFO("Returning empty vector.");
      return std::vector<std::vector<size_t>>(0);
    }
    if (t_k >= m_N) {
			m_log->ERROR("KDTree " + m_name
									+ ": Attempted to query " + std::to_string(t_k)
									+ " neighbors for points in array m_points of size "
									+ std::to_string(m_N));
      m_log->INFO("Setting k to 3 for this search.");
      t_k = 3;
		}
		std::vector<std::vector<size_t>> neighbors(t_points.size());
    std::vector<std::vector<double>> distances(t_points.size());
    m_log->INFO("KDTree " + m_name + ": Querying a set of points for the "
                + "nearest " + std::to_string(t_k) + " neighbors.");

    const size_t num_results = t_k;
    for (auto i = 0; i < t_points.size(); i++) {
			neighbors[i].resize(num_results);
      distances[i].resize(num_results);
      std::vector<size_t> ret_indexes(num_results);
      std::vector<double> out_dists_sqr(num_results);

      nanoflann::KNNResultSet<double> resultSet(num_results);
      resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
			m_kdtree->index->findNeighbors(resultSet, &t_points[i][0],
				                       nanoflann::SearchParams(10));
			neighbors[i] = std::move(ret_indexes);
    }
		return neighbors;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<std::vector<double>>
  KDTree<T>::queryDistances(const std::vector<std::vector<T>>& t_points,
                            size_t t_k)
  {
    //  check that t_point has the same dimension as m_points.
    if (t_points[0].size() != m_dim) {
      m_log->ERROR("Attempted to search a tree of dimension "
                   + std::to_string(m_dim) + " with points of dimension "
                   + std::to_string(t_points[0].size()));
      m_log->INFO("Returning empty vector.");
      return std::vector<std::vector<double>>(0);
    }
    if (t_k >= m_N) {
			m_log->ERROR("KDTree " + m_name
									+ ": Attempted to query " + std::to_string(t_k)
									+ " neighbors for points in array m_points of size "
									+ std::to_string(m_N));
      m_log->INFO("Setting k to 3 for this search.");
      t_k = 3;
		}
		std::vector<std::vector<double>> distances(t_points.size());
    m_log->INFO("KDTree " + m_name + ": Querying a set of points for the "
                + "nearest " + std::to_string(t_k) + " neighbors.");

    const size_t num_results = t_k;
    for (auto i = 0; i < t_points.size(); i++) {
			distances[i].resize(num_results);
      std::vector<size_t> ret_indexes(num_results);
      std::vector<double> out_dists_sqr(num_results);

      nanoflann::KNNResultSet<double> resultSet(num_results);
      resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
			m_kdtree->index->findNeighbors(resultSet, &t_points[i][0],
				                       nanoflann::SearchParams(10));
      // convert distances squared to distances
      std::vector<double> out_dists(out_dists_sqr.size());
      for (auto j = 0; j < out_dists_sqr.size(); j++) {
       out_dists[j] = sqrt(out_dists_sqr[j]);
      }
			distances[i] = std::move(out_dists);
    }
		return distances;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::tuple<std::vector<std::vector<size_t>>,std::vector<std::vector<double>>>
  KDTree<T>::query(const std::vector<std::vector<T>>& t_points,
                            size_t t_k)
  {
    //  check that t_point has the same dimension as m_points.
    if (t_points[0].size() != m_dim) {
      m_log->ERROR("Attempted to search a tree of dimension "
                   + std::to_string(m_dim) + " with points of dimension "
                   + std::to_string(t_points[0].size()));
      m_log->INFO("Returning empty vector.");
      return std::tuple<std::vector<std::vector<size_t>>,
                        std::vector<std::vector<double>>>();
    }
    if (t_k >= m_N) {
			m_log->ERROR("KDTree " + m_name
									+ ": Attempted to query " + std::to_string(t_k)
									+ " neighbors for points in array m_points of size "
									+ std::to_string(m_N));
      m_log->INFO("Setting k to 3 for this search.");
      t_k = 3;
		}
		std::vector<std::vector<size_t>> neighbors(t_points.size());
    std::vector<std::vector<double>> distances(t_points.size());
    m_log->INFO("KDTree " + m_name + ": Querying a set of points for the "
                + "nearest " + std::to_string(t_k) + " neighbors.");

    const size_t num_results = t_k;
    for (auto i = 0; i < t_points.size(); i++) {
			neighbors[i].resize(num_results);
      distances[i].resize(num_results);
      std::vector<size_t> ret_indexes(num_results);
      std::vector<double> out_dists_sqr(num_results);

      nanoflann::KNNResultSet<double> resultSet(num_results);
      resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
			m_kdtree->index->findNeighbors(resultSet, &t_points[i][0],
				                       nanoflann::SearchParams(10));
      std::vector<double> out_dists(out_dists_sqr.size());
      for (auto i = 0; i < out_dists_sqr.size(); i++) {
        out_dists[i] = sqrt(out_dists_sqr[i]);
      }
			neighbors[i] = std::move(ret_indexes);
      distances[i] = std::move(out_dists);
    }
    std::tuple<std::vector<std::vector<size_t>>,std::vector<std::vector<double>>>
    result {neighbors,distances};
		return result;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<size_t>
  KDTree<T>::queryNeighbors(const std::vector<T>& t_point, double t_radius)
  {
    //  check that t_point has the same dimension as m_points.
    if (t_point.size() != m_dim) {
      m_log->ERROR("Attempted to search a tree of dimension "
                   + std::to_string(m_dim) + " with a point of dimension "
                   + std::to_string(t_point.size()));
      m_log->INFO("Returning empty vector.");
      return std::vector<size_t>(0);
    }
    //  make sure that radius is positive
    if (t_radius < 0) {
      m_log->WARN("Attempted to use a negative search radius "
                  + std::to_string(t_radius) + " in nearest neighbor search.");
      m_log->INFO("Setting radius to its absolute value.");
      t_radius = abs(t_radius);
    }

		m_log->INFO("KDTree " + m_name + ": Querying a point for the nearest "
							  + "neighbors within a radius of " + std::to_string(t_radius));
    //  WARNING: Nanoflann uses radius^2 for searchs since it is
    //  computationally simpler than taking the square root for every
    //  iteration.
    t_radius *= t_radius;
    std::vector<std::pair<size_t,double>> ret_matches;

  	m_kdtree->index->radiusSearch(&t_point[0], t_radius, ret_matches,
  		                      nanoflann::SearchParams(10));
    std::vector<size_t> indices(ret_matches.size());
    for (auto j = 0; j < ret_matches.size(); j++) {
      indices[j] = ret_matches[j].first;
    }
		return indices;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<double>
  KDTree<T>::queryDistances(const std::vector<T>& t_point, double t_radius)
  {
    //  check that t_point has the same dimension as m_points.
    if (t_point.size() != m_dim) {
      m_log->ERROR("Attempted to search a tree of dimension "
                   + std::to_string(m_dim) + " with a point of dimension "
                   + std::to_string(t_point.size()));
      m_log->INFO("Returning empty vector.");
      return std::vector<double>(0);
    }
    //  make sure that radius is positive
    if (t_radius < 0) {
      m_log->WARN("Attempted to use a negative search radius "
                  + std::to_string(t_radius) + " in nearest neighbor search.");
      m_log->INFO("Setting radius to its absolute value.");
      t_radius = abs(t_radius);
    }

		m_log->INFO("KDTree " + m_name + ": Querying a point for the nearest "
							  + "neighbors within a radius of " + std::to_string(t_radius));
    //  WARNING: Nanoflann uses radius^2 for searchs since it is
    //  computationally simpler than taking the square root for every
    //  iteration.
    t_radius *= t_radius;
    std::vector<std::pair<size_t,double>> ret_matches;

  	m_kdtree->index->radiusSearch(&t_point[0], t_radius, ret_matches,
  		                      nanoflann::SearchParams(10));
		std::vector<double> distances(ret_matches.size());
		for (auto j = 0; j < ret_matches.size(); j++) {
			distances[j] = ret_matches[j].second;
		}
		return distances;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::tuple<std::vector<size_t>,std::vector<double>>
  KDTree<T>::query(const std::vector<T>& t_point, double t_radius)
  {
    //  check that t_point has the same dimension as m_points.
    if (t_point.size() != m_dim) {
      m_log->ERROR("Attempted to search a tree of dimension "
                   + std::to_string(m_dim) + " with a point of dimension "
                   + std::to_string(t_point.size()));
      m_log->INFO("Returning empty vector.");
      return std::tuple<std::vector<size_t>,std::vector<double>>();
    }
    //  make sure that radius is positive
    if (t_radius < 0) {
      m_log->WARN("Attempted to use a negative search radius "
                  + std::to_string(t_radius) + " in nearest neighbor search.");
      m_log->INFO("Setting radius to its absolute value.");
      t_radius = abs(t_radius);
    }

		m_log->INFO("KDTree " + m_name + ": Querying a point for the nearest "
							  + "neighbors within a radius of " + std::to_string(t_radius));
    //  WARNING: Nanoflann uses radius^2 for searchs since it is
    //  computationally simpler than taking the square root for every
    //  iteration.
    t_radius *= t_radius;
    std::vector<std::pair<size_t,double>> ret_matches;

  	m_kdtree->index->radiusSearch(&t_point[0], t_radius, ret_matches,
  		                      nanoflann::SearchParams(10));
    std::vector<size_t> indices(ret_matches.size());
    std::vector<double> distances(ret_matches.size());
    for (auto j = 0; j < ret_matches.size(); j++) {
      indices[j] = ret_matches[j].first;
      distances[j] = sqrt(ret_matches[j].second);
    }
    std::tuple<std::vector<size_t>,std::vector<double>>
    result {indices,distances};
		return result;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<std::vector<size_t>>
  KDTree<T>::queryNeighbors(const std::vector<std::vector<T>>& t_points,
                            double t_radius)
  {
    //  check that t_point has the same dimension as m_points.
    if (t_points[0].size() != m_dim) {
      m_log->ERROR("Attempted to search a tree of dimension "
                   + std::to_string(m_dim) + " with points of dimension "
                   + std::to_string(t_points[0].size()));
      m_log->INFO("Returning empty vector.");
      return std::vector<std::vector<size_t>>(0);
    }
    //  make sure that radius is positive
    if (t_radius < 0) {
      m_log->WARN("Attempted to use a negative search radius "
                  + std::to_string(t_radius) + " in nearest neighbor search.");
      m_log->INFO("Setting radius to its absolute value.");
      t_radius = abs(t_radius);
    }

		m_log->INFO("KDTree " + m_name + ": Querying a point for the nearest "
							  + "neighbors within a radius of " + std::to_string(t_radius));
    //  WARNING: Nanoflann uses radius^2 for searchs since it is
    //  computationally simpler than taking the square root for every
    //  iteration.
    t_radius *= t_radius;
    std::vector<std::pair<size_t,double>> ret_matches;
    std::vector<std::vector<size_t>> indices;
    for (auto i = 0; i < t_points.size(); i++) {
  	  m_kdtree->index->radiusSearch(&t_points[i][0], t_radius, ret_matches,
  		                        nanoflann::SearchParams(10));
      std::vector<size_t> temp_indices(ret_matches.size());
  		for (auto j = 0; j < ret_matches.size(); j++) {
  			temp_indices[j] = ret_matches[j].first;
  		}
      indices[i] = std::move(temp_indices);
    }
		return indices;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<std::vector<double>>
  KDTree<T>::queryDistances(const std::vector<std::vector<T>>& t_points,
                            double t_radius)
  {
    //  check that t_point has the same dimension as m_points.
    if (t_points[0].size() != m_dim) {
      m_log->ERROR("Attempted to search a tree of dimension "
                   + std::to_string(m_dim) + " with points of dimension "
                   + std::to_string(t_points[0].size()));
      m_log->INFO("Returning empty vector.");
      return std::vector<std::vector<double>>(0);
    }
    //  make sure that radius is positive
    if (t_radius < 0) {
      m_log->WARN("Attempted to use a negative search radius "
                  + std::to_string(t_radius) + " in nearest neighbor search.");
      m_log->INFO("Setting radius to its absolute value.");
      t_radius = abs(t_radius);
    }

		m_log->INFO("KDTree " + m_name + ": Querying a point for the nearest "
							  + "neighbors within a radius of " + std::to_string(t_radius));
    //  WARNING: Nanoflann uses radius^2 for searchs since it is
    //  computationally simpler than taking the square root for every
    //  iteration.
    t_radius *= t_radius;
    std::vector<std::pair<size_t,double>> ret_matches;
    std::vector<std::vector<double>> distances;
    for (auto i = 0; i < t_points.size(); i++) {
  	  m_kdtree->index->radiusSearch(&t_points[i][0], t_radius, ret_matches,
  		                        nanoflann::SearchParams(10));
      std::vector<double> temp_distances(ret_matches.size());
  		for (auto j = 0; j < ret_matches.size(); j++) {
  			temp_distances[j] = ret_matches[j].second;
  		}
      distances[i] = std::move(temp_distances);
    }
		return distances;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::tuple<std::vector<std::vector<size_t>>,std::vector<std::vector<double>>>
  KDTree<T>::query(const std::vector<std::vector<T>>& t_points,
                   double t_radius)
  {
    //  check that t_point has the same dimension as m_points.
    if (t_points[0].size() != m_dim) {
      m_log->ERROR("Attempted to search a tree of dimension "
                   + std::to_string(m_dim) + " with points of dimension "
                   + std::to_string(t_points[0].size()));
      m_log->INFO("Returning empty vector.");
      return std::tuple<std::vector<std::vector<size_t>>,
                        std::vector<std::vector<double>>>();
    }
    //  make sure that radius is positive
    if (t_radius < 0) {
      m_log->WARN("Attempted to use a negative search radius "
                  + std::to_string(t_radius) + " in nearest neighbor search.");
      m_log->INFO("Setting radius to its absolute value.");
      t_radius = abs(t_radius);
    }

		m_log->INFO("KDTree " + m_name + ": Querying a point for the nearest "
							  + "neighbors within a radius of " + std::to_string(t_radius));
    //  WARNING: Nanoflann uses radius^2 for searchs since it is
    //  computationally simpler than taking the square root for every
    //  iteration.
    t_radius *= t_radius;
    std::vector<std::pair<size_t,double>> ret_matches;
    std::vector<std::vector<size_t>> indices;
    std::vector<std::vector<double>> distances;
    for (auto i = 0; i < t_points.size(); i++) {
  	  m_kdtree->index->radiusSearch(&t_points[i][0], t_radius, ret_matches,
  		                        nanoflann::SearchParams(10));
      std::vector<size_t> temp_indices(ret_matches.size());
      std::vector<double> temp_distances(ret_matches.size());
  		for (auto j = 0; j < ret_matches.size(); j++) {
  			temp_indices[j] = ret_matches[j].first;
        temp_distances[j] = sqrt(ret_matches[j].second);
  		}
      indices[i] = std::move(temp_indices);
      distances[i] = std::move(temp_distances);
    }
    std::tuple<std::vector<std::vector<size_t>>,std::vector<std::vector<double>>>
    result {indices,distances};
		return result;
  }
  //----------------------------------------------------------------------------
}
