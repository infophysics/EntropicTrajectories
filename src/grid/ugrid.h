//------------------------------------------------------------------------------
//  ugrid.h
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
#include "utils.h"
#include "matrix.h"
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
    UGrid(uint64_t dim);
    UGrid(std::string name, uint64_t dim);
    UGrid(uint64_t dim, uint64_t N);
    UGrid(std::string name, uint64_t dim, uint64_t N);
    UGrid(std::vector<T> ugrid);
    UGrid(std::vector<std::vector<T>> ugrid);
    //--------------------------------------------------------------------------
    //  Constructors with shared loggers
    //--------------------------------------------------------------------------
    UGrid(std::shared_ptr<Log> log);
    UGrid(uint64_t dim, std::shared_ptr<Log> log);
    UGrid(std::string name, uint64_t dim, std::shared_ptr<Log> log);
    UGrid(uint64_t dim, uint64_t N, std::shared_ptr<Log> log);
    UGrid(std::string name, uint64_t dim, uint64_t N, std::shared_ptr<Log> log);
    UGrid(std::vector<T> ugrid, std::shared_ptr<Log> log);
    UGrid(std::vector<std::vector<T>> ugrid, std::shared_ptr<Log> log);
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Getters
    //--------------------------------------------------------------------------
    uint64_t getDim();
    uint64_t getN();
    std::vector<std::vector<T> >  getUGrid();
    std::string getName();
    std::vector<std::vector<size_t>> getNeighbors();
    std::vector<std::vector<double>> getDistances();
    std::vector<std::vector<size_t>> getNeighborsRadius();
    std::vector<std::vector<double>> getDistancesRadius();
    std::vector<size_t> getNeighbors(uint64_t index);
    std::shared_ptr<Log> getLogger();
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Setters
    //--------------------------------------------------------------------------
    void setDim(uint64_t dim);
    void setN(uint64_t N);
    void setUGrid(std::vector<std::vector<T>> ugrid);
    void setName(std::string name);
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Access operators for ugrid
    //--------------------------------------------------------------------------
    T& operator()(const uint64_t i, const uint64_t j);
    const T& operator()(const uint64_t i, const uint64_t j) const;
    std::vector<T>& operator()(const uint64_t i);
    const std::vector<T>& operator()(const uint64_t i) const;
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Points and projections
    //--------------------------------------------------------------------------
    const std::vector<T>& getPoint(const uint64_t i) const;
    std::vector<T> projection(uint64_t j);
    void setPoint(uint64_t i, std::vector<T> p);
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  KDTree methods
    //--------------------------------------------------------------------------
    void queryNeighbors(uint64_t k);
    std::vector<size_t> queryNeighbors(const std::vector<T>& point, uint64_t k);
    std::vector<std::vector<size_t>>
    queryNeighbors(const std::vector<std::vector<T>>& points, uint64_t k);
    void queryRadius(double radius);
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Various functions
    //--------------------------------------------------------------------------
    bool checkConsistency();
    //--------------------------------------------------------------------------

  private:
    //--------------------------------------------------------------------------
    //  Basic attributes
    //--------------------------------------------------------------------------
    std::string _name;                 //  Name of the grid
    uint64_t _dim;                     //  Dimension of the grid
    uint64_t _N;                       //  Number of points in the grid
    std::vector<std::vector<T>> _ugrid;//  Vector of vectors array
    std::vector<std::string> _coords;  //  Coordinate labels
    //--------------------------------------------------------------------------
    //  KDTree
    //--------------------------------------------------------------------------
    std::shared_ptr<KDTreeVectorOfVectorsAdaptor<
                     std::vector<std::vector<T>>,T>> _kdtree;
    //--------------------------------------------------------------------------
    //  Results from KDTree searches
    //--------------------------------------------------------------------------
    int _searchFlag;                   //  flag for changes
    uint32_t _k;                       //  number of neighbors to search
    double _radius;                    //  search radius
    std::vector<std::vector<size_t>> _neighbors;
    std::vector<std::vector<double>> _distances;
    std::vector<std::vector<size_t>> _neighbors_radius;
    std::vector<std::vector<double>> _distances_radius;
    //--------------------------------------------------------------------------
    //  Logger
    //--------------------------------------------------------------------------
    std::shared_ptr<Log> _log;
    //--------------------------------------------------------------------------
  };
  //----------------------------------------------------------------------------

  //  Explicit construction of type double
  template class UGrid<double>;
}
