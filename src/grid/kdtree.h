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
#include "utils.h"
#include "matrix.h"
#include "log.h"
#include <nanoflann.hpp>
#include "KDTreeVectorOfVectorsAdaptor.h"

namespace ET
{
  template<typename T>
  class KDTree : public std::enable_shared_from_this<KDTree<T>>
  {
  public:
    KDTree();
    ~KDTree();
    KDTree(std::shared_ptr<std::vector<std::vector<T>>> t_points);

    //--------------------------------------------------------------------------
    //  Getters
    //--------------------------------------------------------------------------
    size_t getDim();
    size_t getN();
    std::shared_ptr<std::vector<std::vector<T>>>  getPoints();
    std::string getName();
    std::vector<std::vector<size_t>> getNeighbors();
    std::vector<std::vector<double>> getDistances();
    std::vector<std::vector<size_t>> getNeighborsRadius();
    std::vector<std::vector<double>> getDistancesRadius();
    std::vector<size_t> getNeighbors(size_t t_index);
    std::shared_ptr<Log> getLogger();
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

  private:
    std::string m_name{""};
    size_t m_N{0};
    size_t m_dim{0};
    std::shared_ptr<std::vector<std::vector<T>>> m_points;
    std::shared_ptr<KDTreeVectorOfVectorsAdaptor<
                     std::vector<std::vector<T>>,T>> m_KDTree;
    std::shared_ptr<Log> m_log;
    //--------------------------------------------------------------------------
    //  Results from KDTree searches
    //--------------------------------------------------------------------------
    int m_searchFlag;          //  flag for changes
    size_t m_k;               //  number of neighbors to search
    double m_radius;          //  search radius
    std::vector<std::vector<size_t>> m_neighbors;
    std::vector<std::vector<double>> m_distances;
    std::vector<std::vector<size_t>> m_neighbors_radius;
    std::vector<std::vector<double>> m_distances_radius;
    //--------------------------------------------------------------------------
  };

  template class KDTree<double>;
}
