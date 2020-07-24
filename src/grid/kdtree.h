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
#pragma once

#include <vector>
#include <string>
#include <map>

#include "utilities.h"
#include "matrix.h"
#include "log.h"

#include <nanoflann.hpp>
#include "KDTreeVectorOfVectorsAdaptor.h"

namespace ET
{
  //! \enum KDTreeBackend
  /*! An enum corresponding to the choice of backend for the kdtree.
   */
  enum class KDTreeBackend
  {
    NANOFLANN,
    FLANN,
  };
  //! \enum KDTreeSearchFlags
  /*! An enum corresponding to behavior of the KDTree
   */
  enum class KDTreeSearchFlags
  {
    OUTDATED,
    BUILT,
    NEW_K,
    NEW_RADIUS,
    NEW_POINTS,
    CURRENT_K,
    CURRENT_RADIUS,
  };

  extern std::map<std::string, KDTreeBackend> KDTreeBackendMap;
  extern std::map<KDTreeBackend, std::string> KDTreeBackendNameMap;

  //! \class KDTree
  /*! This is a wrapper for various KDTree implementations such as
   *  nanoflann.  */
  template<typename T>
  class KDTree : public std::enable_shared_from_this<KDTree<T>>
  {
  public:
    //! Default Constructor
    /*! Default constructor for a KDTree.
     */
    KDTree();
    //! Destructor
    /*! Destructor for the KDTree
     */
    ~KDTree();
    //! Constructor
    /*! Constructor taking a shared pointer of the set of points.
     *  @param t_points A std::shared_ptr of a std::vector<std::vector<T>>
     */
    KDTree(const std::shared_ptr<std::vector<std::vector<T>>> t_points);
    //! Constructor
    /*! Constructor taking a name and a shared pointer of the set of points.
     *  @param t_name An std::string for the name of the object.
     *  @param t_points A std::shared_ptr of a std::vector<std::vector<T>>
     */
    KDTree(const std::string t_name,
           const std::shared_ptr<std::vector<std::vector<T>>> t_points);

    //  Getters and Setters
    //! Get Backend
    /*! Get the backend enum for the KDTree.
     *  @return A KDTreeBackend enum.
     */
    enum KDTreeBackend getBackend() const;
    //! Get name.
    /*! Get the name of the KDTree.
     *  @return m_name The name of the KDTree.
     */
    std::string getName() const;
    //! Get dimension
    /*! Get the dimension of underlying set of points.
     *  @return A size_t object denoting the dimension of the points.
     */
    size_t getDim() const;
    //! Get N
    /*! Get the number of points.
     *  @return A size_t object denoting the number of points.
     */
    size_t getN() const;
    //! Get Points
    /*! Get the shared pointer to the underlying points.
     *  @return A std::shared_ptr of the std::vector<std::vector<T>> of points.
     */
    std::shared_ptr<std::vector<std::vector<T>>>  getPoints() const;
    //! Get current neighbor indices
    /*! Get the values of the indices for the last search that was performed.
     *  @return The list of indices for the last nearest neighbor search.
     */
    std::vector<std::vector<size_t>> getCurrentNeighborIndices();
    //! Get current neighbor distances
    /*! Get the values of the distances for the last search that was performed.
     *  @return The list of distances for the last nearest neighbor search.
     */
    std::vector<std::vector<double>> getCurrentNeighborDistances();
    //! Get current neighbor indices
    /*! Get the values of the indices for the last search that was performed
     *  for a particular point in t_points.
     *  @param t_i Index to the point of interest in t_points.
     *  @return The list of indices for the last nearest neighbor search.
     */
    std::vector<size_t> getCurrentNeighborIndices(const size_t t_i);
    //! Get current neighbor distances
    /*! Get the values of the distances for the last search that was performed
     *. for a particular point in t_points.
     *  @param t_i Index to the point of interest in t_points.
     *  @return The list of distances for the last nearest neighbor search.
     */
    std::vector<double> getCurrentNeighborDistances(const size_t t_i);
    //! Get current global k
    /*! Get the current global k value.
     *  @return The current global k value.
     */
    size_t getCurrentGlobalK() const;
    //! Get current global radius
    /*! Get the current global radius value.
     *  @return The current global radius value.
     */
    double getCurrentGlobalRadius() const;
    //! Get Log
    /*! Get the shared logger for this Grid.
     *  @return m_log The shared instance of the logger.
     */
    std::shared_ptr<Log> getLog() const;
    //! Get Search Flag
    /*! Get the current search flag
     *  @return An enum for the search flag.
     */
    enum KDTreeSearchFlags getSearchFlag() const;
    //! Set Backend.
    /*! Set the backend for the KDTree.
     *  @param t_backend A KDTreeBackend enum.
     */
    void setBackend(enum KDTreeBackend t_backend);
    //! Set name.
    /*! Set the name of the Grid.
     *  @param t_name A std::string specifying the name of the Grid.
     */
    void setName(const std::string t_name);
    //! Set dim
    /*! Set the dimension of the Grid.
     *  @param t_dim A size_t object specifying the dimension of the Grid.
     */
    void setDim(const size_t t_dim);
    //! Set N
    /*! Set the number of elements in the Grid.
     *  @param t_N A size_t object specifying the number of
     *  elements in the Grid.
     */
    void setN(const size_t t_N);
    //! Set Points
    /*! Set the shared pointer for the underlying points.
     *  @param t_points An std::shared_ptr of the std::vector<std::vector<T>>
     *  of points.
     */
    void setPoints(const std::shared_ptr<std::vector<std::vector<T>>> t_points);
    //! Set Log
    /*! Set the shared logger for this Grid.
     *  @param t_log A std::shared_ptr<Log> instance of a logger.
     */
    void setLog(std::shared_ptr<Log> t_log);
    //! Set current global k
    /*! Set the current global k value.
     *  @param t_k A size_t object denoting the number of global nearest
     *  neighbors to use.
     */
    void setCurrentGlobalK(const size_t t_k);
    //! Set current global radius
    /*! Set the current global radius value.
     *  @param t_radius A double object denoting the radius to use.
     */
    void setCurrentGlobalRadius(const double t_radius);
    //  KDTree methods
    //! Setup tree
    /*! Initialize the KDTree object using m_points.
     */
    void setupTree();
    //! Query Neighbors
    /*! Query the tree for the nearest t_k neighbors of every point and
     *  store the result in m_currentNeighborIndices and
     *  m_currentNeighborDistances.
     *  @param t_k The number of nearest neighbors to search for.
     */
    void queryNeighbors(size_t t_k);
    //! Query Neighbors
    /*! Query the tree for the nearest neighbors of every point
     *  within the radius t_radius and store the result in
     *  m_currentNeighborIndices and m_currentNeighborDistances.
     *  @param t_radius The search radius to use for nearest neighbors.
     */
    void queryNeighbors(double t_radius);
    //! Query neighbors
    /*! Query the tree for t_k nearest neighbors of a particular point
     *  and return the indices of the neighbors.
     *  @param t_point The point to search near.
     *  @param t_k The number of neighbors to search for.
     *  @return An std::vector<size_t> of the neighbors indices.
     */
    std::vector<size_t>
    queryNeighbors(const std::vector<T>& t_point, size_t t_k);
    //! Query neighbors
    /*! Query the tree for t_k nearest neighbors of a particular point
     *  and return the distances to each neighbor.
     *  @param t_point The point to search near.
     *  @param t_k The number of neighbors to search for.
     *  @return An std::vector<double> of the neighbors distances.
     */
    std::vector<double>
    queryDistances(const std::vector<T>& t_point, size_t t_k);
    //! Query
    /*! Query the tree for t_k nearest neighbors of a particular point
     *  and return the indices and distances for each neighbor.
     *  @param t_point The point to search near.
     *  @param t_k The number of neighbors to search for.
     *  @return An std::tuple<std::vector<size_t>,std::vector<double>>
     *  of indices and distances.
     */
    std::tuple<std::vector<size_t>,std::vector<double>>
    query(const std::vector<T>& t_point, size_t t_k);
    //! Query neighbors
    /*! Query the tree for t_k nearest neighbors of a set of points
     *  and return the indices of the neighbors.
     *  @param t_points The points to search near.
     *  @param t_k The number of neighbors to search for.
     *  @return An std::vector<std::vector<size_t>> of the neighbors indices.
     */
    std::vector<std::vector<size_t>>
    queryNeighbors(const std::vector<std::vector<T>>& t_points, size_t t_k);
    //! Query distances
    /*! Query the tree for t_k nearest neighbors of a set of points
     *  and return the distances of the neighbors.
     *  @param t_points The points to search near.
     *  @param t_k The number of neighbors to search for.
     *  @return An std::vector<std::vector<double>> of the neighbors distances.
     */
    std::vector<std::vector<double>>
    queryDistances(const std::vector<std::vector<T>>& t_points, size_t t_k);
    //! Query
    /*! Query the tree for t_k nearest neighbors of a set of points
     *  and return the indices and distances for each neighbor.
     *  @param t_points The points to search near.
     *  @param t_k The number of neighbors to search for.
     *  @return An std::tuple<std::vector<std::vector<size_t>>,
     *  std::vector<std::vector<double>>> of indices and distances.
     */
    std::tuple<std::vector<std::vector<size_t>>,
               std::vector<std::vector<double>>>
    query(const std::vector<std::vector<T>>& t_points, size_t t_k);
    //! Query neighbors
    /*! Query the tree for nearest neighbors of a particular point
     *  within a radius t_radius and return the indices of the neighbors.
     *  @param t_point The point to search near.
     *  @param t_radius The radius to search within.
     *  @return An std::vector<size_t> of the neighbors indices.
     */
    std::vector<size_t>
    queryNeighbors(const std::vector<T>& t_point, double t_radius);
    //! Query neighbors
    /*! Query the tree for nearest neighbors of a particular point
     *  within a radius t_radius and return the distances to each neighbor.
     *  @param t_point The point to search near.
     *  @param t_radius The radius to search within.
     *  @return An std::vector<double> of the neighbors distances.
     */
    std::vector<double>
    queryDistances(const std::vector<T>& t_point, double t_radius);
    //! Query
    /*! Query the tree for nearest neighbors of a particular point
     *  within a radius t_radius and return the indices and distances
     *  for each neighbor.
     *  @param t_point The point to search near.
     *  @param t_radius The radius to search within.
     *  @return An std::tuple<std::vector<size_t>,std::vector<double>>
     *  of indices and distances.
     */
    std::tuple<std::vector<size_t>,std::vector<double>>
    query(const std::vector<T>& t_point, double t_radius);
    //! Query neighbors
    /*! Query the tree for nearest neighbors of a set of points
     *  within a radius t_radius and return the indices of the neighbors.
     *  @param t_points The points to search near.
     *  @param t_radius The radius to search within.
     *  @return An std::vector<std::vector<size_t>> of the neighbors indices.
     */
    std::vector<std::vector<size_t>>
    queryNeighbors(const std::vector<std::vector<T>>& t_points, double t_radius);
    //! Query distances
    /*! Query the tree for nearest neighbors of a set of points
     *  within a radius t_radius and return the distances of the neighbors.
     *  @param t_points The points to search near.
     *  @param t_radius The radius to search within.
     *  @return An std::vector<std::vector<double>> of the neighbors distances.
     */
    std::vector<std::vector<double>>
    queryDistances(const std::vector<std::vector<T>>& t_points, double t_radius);
    //! Query
    /*! Query the tree for nearest neighbors of a set of points
     *  within a radius t_radius and return the indices and distances
     *  for each neighbor.
     *  @param t_points The points to search near.
     *  @param t_radius The radius to search within.
     *  @return An std::tuple<std::vector<std::vector<size_t>>,
     *  std::vector<std::vector<double>>> of indices and distances.
     */
    std::tuple<std::vector<std::vector<size_t>>,
               std::vector<std::vector<double>>>
    query(const std::vector<std::vector<T>>& t_points, double t_radius);
    //--------------------------------------------------------------------------

  private:
    /*! Backend.  The backend to use for the KDTree.  Defaulted to NANOFLANN.
     */
    enum KDTreeBackend m_backend {KDTreeBackend::NANOFLANN};
    /*! Name.  Name of the KDTree.  Defaulted to empty string.
     */
    std::string m_name {""};
    /*! Dimension.  Dimension of the underlying points.  Defaulted to 0.
     */
    size_t m_dim {0};
    /*! N.  Number of elements in the underlying points.  Default to 0.
     */
    size_t m_N {0};
    /*! Points.  A shared pointer to the underlying std::vector<std::vector<T>>
     *  of points.
     */
    std::shared_ptr<std::vector<std::vector<T>>> m_points;
    /*! nanoflann tree.  A shared pointer to the underlying nanoflann
     *  implementation of the tree.
     */
    std::shared_ptr<KDTreeVectorOfVectorsAdaptor<
                     std::vector<std::vector<T>>,T>> m_kdtree;
    /*! Logger.  Shared instance of a logger.
    */
    std::shared_ptr<Log> m_log {std::make_shared<Log>()};
    /*! Current k value.  A size_t that records the current global value of
     *  k to use in nearest neighbor searches.
     */
    size_t m_currentGlobalK {1};
    /*! Current radius value.  A double that records the current global value
     *  of the radius to use in nearest neighbor searches.
     */
    double m_currentGlobalRadius {0};
    /*! Current neighbors.  A container for the last set of nearest neighbors
     *  indices that were found in the most recent seach.
     */
    std::vector<std::vector<size_t>> m_currentNeighborIndices {{{0}}};
    /*! Current distances.  A container for the last set of nearest neighbor
     *  distances that were found in the most recent search.
     */
    std::vector<std::vector<double>> m_currentNeighborDistances {{{0}}};
    /*! Search flag.  An enum that classifies the state of the tree.
     */
    enum KDTreeSearchFlags m_searchFlag {KDTreeSearchFlags::OUTDATED};
    /*! Logging system name generator.
     */
    virtual std::string NAME() const {
      return "KDTree:[" + std::to_string(m_id.id) + "]:" + m_name + ":";
    }
    /*! Unique ID for each instance.
     */
    UniqueID m_id;
    //--------------------------------------------------------------------------
  };

  template class KDTree<double>;
}
