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
  class UGrid
  {
  public:
    //--------------------------------------------------------------------------
    //  Constructors
    //--------------------------------------------------------------------------
    UGrid();
    ~UGrid();
    UGrid(size_t t_dim);
    UGrid(std::string t_name, size_t t_dim);
    UGrid(size_t t_dim, size_t t_N);
    UGrid(std::string t_name, size_t t_dim, size_t t_N);
    UGrid(std::vector<T> t_ugrid);
    UGrid(std::vector<std::vector<T>> t_ugrid);
    //--------------------------------------------------------------------------
    //  Constructors with shared loggers
    //--------------------------------------------------------------------------
    UGrid(std::shared_ptr<Log> t_log);
    UGrid(size_t t_dim, std::shared_ptr<Log> t_log);
    UGrid(std::string t_name, size_t t_dim, std::shared_ptr<Log> t_log);
    UGrid(size_t t_dim, size_t t_N, std::shared_ptr<Log> t_log);
    UGrid(std::string t_name, size_t t_dim, size_t t_N, std::shared_ptr<Log> t_log);
    UGrid(std::vector<T> t_ugrid, std::shared_ptr<Log> t_log);
    UGrid(std::vector<std::vector<T>> t_ugrid,
          std::shared_ptr<Log> t_log);
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Getters
    //--------------------------------------------------------------------------
    const size_t getDim();
    const size_t getN();
    std::vector<std::vector<T>>  getUGrid();
    const std::string getName();
    std::vector<std::vector<size_t>> getNeighbors();
    std::vector<std::vector<double>> getDistances();
    std::vector<std::vector<size_t>> getNeighborsRadius();
    std::vector<std::vector<double>> getDistancesRadius();
    std::vector<size_t> getNeighbors(size_t t_index);
    std::shared_ptr<Log> getLogger();
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Setters
    //--------------------------------------------------------------------------
    void setDim(size_t t_dim);
    void setN(size_t t_N);
    void setUGrid(std::vector<std::vector<T>> t_ugrid);
    void setName(std::string t_name);
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
    std::string m_name {""};                    //  Name of the grid
    size_t m_dim {0};                           //  Dimension of the grid
    size_t m_N {0};                             //  Number of points t_in the grid
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
    //  Logger
    //--------------------------------------------------------------------------
    std::shared_ptr<Log> m_log;
    //--------------------------------------------------------------------------
  };
  //----------------------------------------------------------------------------

  //  Explicit construction of type double
  template class UGrid<double>;
}
