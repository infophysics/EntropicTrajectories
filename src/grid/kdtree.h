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
  template<typename T>
  class kdTree : public std::enable_shared_from_this<kdTree<T>>
  {
  public:
    kdTree();
    ~kdTree();
    kdTree(std::shared_ptr<std::vector<std::vector<T>>> points);

    //--------------------------------------------------------------------------
    //  Getters
    //--------------------------------------------------------------------------
    uint64_t getDim();
    uint64_t getN();
    std::shared_ptr<std::vector<std::vector<T>>>  getPoints();
    std::string getName();
    std::vector<std::vector<size_t>> getNeighbors();
    std::vector<std::vector<double>> getDistances();
    std::vector<std::vector<size_t>> getNeighborsRadius();
    std::vector<std::vector<double>> getDistancesRadius();
    std::vector<size_t> getNeighbors(uint64_t index);
    std::shared_ptr<Log> getLogger();
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  KDTree methods
    //--------------------------------------------------------------------------
    void setupTree();

    void queryNeighbors(uint64_t k);

    std::vector<size_t>
    queryNeighbors(const std::vector<T>& point, uint64_t k);

    std::vector<double>
    queryDistances(const std::vector<T>& point, uint64_t k);

    std::vector<std::vector<size_t>>
    queryNeighbors(const std::vector<std::vector<T>>& points, uint64_t k);

    void queryRadius(double radius);
    //--------------------------------------------------------------------------

  private:
    std::string _name;
    uint64_t _N;
    uint32_t _dim;
    std::shared_ptr<std::vector<std::vector<T>>> _points;
    std::shared_ptr<KDTreeVectorOfVectorsAdaptor<
                     std::vector<std::vector<T>>,T>> _kdtree;
    std::shared_ptr<Log> _log;
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
  };

  template class kdTree<double>;
}
