//------------------------------------------------------------------------------
//  ugrid.h
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
    UGrid();
    ~UGrid();
    UGrid(uint64_t dim);
    UGrid(std::string name, uint64_t dim);
    UGrid(uint64_t dim, uint64_t N);
    UGrid(std::string name, uint64_t dim, uint64_t N);
    UGrid(std::vector<T> ugrid);
    UGrid(std::vector<std::vector<T> > ugrid);

    //  Getters
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

    //  Setters
    void setDim(uint64_t dim);
    void setN(uint64_t N);
    void setUGrid(std::vector<std::vector<T> > ugrid);
    void setName(std::string name);

    //  Access operators for ugrid
    T& operator()(const uint64_t i, const uint64_t j);
    const T& operator()(const uint64_t i, const uint64_t j) const;
    //  Access operators for points
    std::vector<T>& operator()(const uint64_t i);
    const std::vector<T>& operator()(const uint64_t i) const;

    //  points and projections
    const std::vector<T>& getPoint(const uint64_t i) const;
    std::vector<T> projection(uint64_t j);
    void setPoint(uint64_t i, std::vector<T> p);

    //  find nearest neighbors
    void queryNeighbors(uint64_t k);
    void queryRadius(double radius);

    //  various functions
    bool checkConsistency();

  private:
    //  dimension of the ugrid
    uint64_t _dim;
    //  number of points in the ugrid
    uint64_t _N;
    //  Unstructured ugrid
    std::vector<std::vector<T> > _ugrid;
    std::string _name;
    std::vector<std::string> _coords;
    //  associated kdTree neighbors
    std::vector<std::vector<size_t> > _neighbors;
    std::vector<std::vector<double> > _distances;
    std::vector<std::vector<size_t> > _neighbors_radius;
    std::vector<std::vector<double> > _distances_radius;
    //KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<T> >, T> _kdt;
    std::shared_ptr<Log> _log;
  };

  template class UGrid<double>;

}
