//------------------------------------------------------------------------------
//  ugrid.h
//  The Entropic Trajectories Framework
//  -----------------------------------
//  Copyright (C) [2020] by [N. Carrara]
//  [ncarrara@albany.edu]
//
//  Permission to use, copy, modify, and/or distribute this software for any
//  purpose with or without fee t_is hereby granted.
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

#include "grid.h"
#include "utils.h"
#include "matrix.h"
#include "kdtree.h"
#include "log.h"

#include <nanoflann.hpp>
#include "KDTreeVectorOfVectorsAdaptor.h"

namespace ET
{
  //----------------------------------------------------------------------------
  //  Class for unstructured grids
  //----------------------------------------------------------------------------
  template<typename T>
  class UGrid : public Grid<T>
  {
  public:
    //  Constructors
    //! Defualt Constructor
    /*! Default constructor for Grid.
     */
    UGrid();
    //! Destructor
    /*! Destructor for UGrid.
     */
    ~UGrid();
    //! Constructor
    /*! constructor for UGrid that takes a Logger
     *  @param t_log A shared logger instance.
     */
    UGrid(std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for UGrid that takes a name
     *  @param t_name An std::string specifying this objects name.
     */
    UGrid(std::string t_name);
    //! Constructor
    /*! constructor for UGrid that takes a name and a Logger
     *  @param t_name An std::string specifying this objects name.
     *  @param t_log A shared logger instance.
     */
    UGrid(std::string t_name, std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for UGrid that takes a dimension.
     *  @param t_dim A size_t object which specifies the dimension.
     */
    UGrid(size_t t_dim);
    //! Constructor
    /*! constructor for UGrid that takes a dimension and a logger.
     *  @param t_dim A size_t object which specifies the dimension.
     *  @param t_log A shared logger instance.
     */
    UGrid(size_t t_dim, std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for UGrid that takes a name and a dimension.
     *  @param t_name An std::string specifying this objects name.
     *  @param t_dim A size_t object which specifies the dimension.
     */
    UGrid(std::string t_name, size_t t_dim);
    //! Constructor
    /*! constructor for UGrid that takes a name, a dimension and a logger.
     *  @param t_name An std::string specifying this objects name.
     *  @param t_dim A size_t object which specifies the dimension.
     *  @param t_log A shared logger instance.
     */
    UGrid(std::string t_name, size_t t_dim, std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for UGrid that takes a dimension and a number of elements.
     *  @param t_dim A size_t object which specifies the dimension.
     *  @param t_N A size_t object which specifies the number of elements.
     */
    UGrid(size_t t_dim, size_t t_N);
    //! Constructor
    /*! constructor for UGrid that takes a dimension, a number of elements
     *  and a logger.
     *  @param t_dim A size_t object which specifies the dimension.
     *  @param t_log A shared logger instance.
     */
    UGrid(size_t t_dim, size_t t_N, std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for UGrid that takes a name, a dimension and
     *  a number of elements.
     *  @param t_name An std::string specifying this objects name.
     *  @param t_dim A size_t object which specifies the dimension.
     */
    UGrid(std::string t_name, size_t t_dim, size_t t_N);
    //! Constructor
    /*! constructor for UGrid that takes a name, a dimension,
     *  a number of elementsand a logger.
     *  @param t_name An std::string specifying this objects name.
     *  @param t_dim A size_t object which specifies the dimension.
     *  @param t_log A shared logger instance.
     */
    UGrid(std::string t_name, size_t t_dim, size_t t_N,
         std::shared_ptr<Log> t_log);

    //   Getters and Setters
    //! Get UGrid
    /*! Get the ugrid array.
     *  @return The ugrid array.
     */
    std::vector<std::vector<T>> getUGrid() const;

    //! Set UGrid
    /*! Sets the UGrid array.
     *  @param t_ugrid A std::vector<std::vector<T>> of the array.
     */
    void setUGrid(std::vector<std::vector<T>> t_ugrid);

    UGrid(std::vector<T> t_ugrid);
    UGrid(std::vector<std::vector<T>> t_ugrid);
    //--------------------------------------------------------------------------
    //  Constructors with shared loggers
    //--------------------------------------------------------------------------

    UGrid(std::vector<T> t_ugrid, std::shared_ptr<Log> t_log);
    UGrid(std::vector<std::vector<T>> t_ugrid,
          std::shared_ptr<Log> t_log);
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Getters
    //--------------------------------------------------------------------------

    std::vector<std::vector<size_t>> getNeighbors();
    std::vector<std::vector<double>> getDistances();
    std::vector<std::vector<size_t>> getNeighborsRadius();
    std::vector<std::vector<double>> getDistancesRadius();
    std::vector<size_t> getNeighbors(size_t t_index);
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Access operators for t_ugrid
    //--------------------------------------------------------------------------
    T& operator()(const size_t t_i, const size_t t_j);
    const T& operator()(const size_t t_i, const size_t t_j) const;
    std::vector<T>& operator()(const size_t t_i);
    const std::vector<T>& operator()(const size_t t_i) const;
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Points and projections
    //--------------------------------------------------------------------------
    const std::vector<T>& getPoint(const size_t t_i) const;
    std::vector<T> projection(size_t t_j);
    void setPoint(size_t t_i, std::vector<T> t_p);
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  KDTree methods
    //--------------------------------------------------------------------------
    void setupTree();
    void queryNeighbors(size_t t_k);

    std::vector<size_t>
    queryNeighbors(const std::vector<T>& t_point, size_t t_k);

    std::vector<double>
    queryDistances(const std::vector<T>& t_point, size_t t_k);

    std::vector<std::vector<size_t>>
    queryNeighbors(const std::vector<std::vector<T>>& t_points, size_t t_k);

    void queryRadius(double t_radius);
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Various functions
    //--------------------------------------------------------------------------
    bool checkConsistency();
    //--------------------------------------------------------------------------

  private:
    std::vector<std::vector<T>> m_ugrid {{{0}}};//  Vector of vectors array
    std::vector<std::string> m_coords {{""}};   //  Coordinate labels
    //--------------------------------------------------------------------------
    //  KDTree
    //--------------------------------------------------------------------------
    KDTree<T> m_kdt;
    std::shared_ptr<KDTreeVectorOfVectorsAdaptor<
                     std::vector<std::vector<T>>,T>> m_KDTree;
    //--------------------------------------------------------------------------
    //  Results from KDTree searches
    //--------------------------------------------------------------------------
    int m_searchFlag;                   //  flag for changes
    size_t m_k;                       //  number of neighbors to search
    double m_radius;                    //  search radius
    std::vector<std::vector<size_t>> m_neighbors;
    std::vector<std::vector<double>> m_distances;
    std::vector<std::vector<size_t>> m_neighbors_radius;
    std::vector<std::vector<double>> m_distances_radius;
    //--------------------------------------------------------------------------
  };
  //----------------------------------------------------------------------------

  //  Explicit construction of type double
  template class UGrid<double>;
}
