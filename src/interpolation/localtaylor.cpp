//------------------------------------------------------------------------------
//  localtaylor.cpp
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
#include "radialbasis.h"

namespace ET
{
  std::map<std::string, SearchScheme> SearchSchemeMap =
  {
    { "NEAREST_NEIGHBORS", SearchScheme::NEAREST_NEIGHBORS },
    { "RADIUS", SearchScheme::RADIUS },
  };

  template<typename T>
  LocalTaylorInterpolator<T>::LocalTaylorInterpolator()
  : Interpolator(), m_name("default")
  {
  }
  template<typename T>
  LocalTaylorInterpolator<T>::~LocalTaylorInterpolator()
  {
  }
  template<typename T>
  LocalTaylorInterpolator<T>::LocalTaylorInterpolator
  (std::shared_ptr<UGrid<T>> t_ugrid) : Interpolator(t_ugrid)
  {
  }
  template<typename T>
  LocalTaylorInterpolator<T>::LocalTaylorInterpolator
  (std::shared_ptr<Log> t_log) : Interpolator(t_log)
  {
  }
  template<typename T>
  LocalTaylorInterpolator<T>::LocalTaylorInterpolator
  (std::shared_ptr<UGrid<T>> t_ugrid, std::shared_ptr<Log> t_log)
  : Interpolator(t_ugrid,t_log)
  {
  }

  template<typename T>
  std::string LocalTaylorInterpolator<T>::getName()
  {
    return m_name;
  }
  template<typename T>
  size_t LocalTaylorInterpolator<T>::getK()
  {
    return m_k;
  }
  template<typename T>
  size_t LocalTaylorInterpolator<T>::getN()
  {
    return m_n;
  }
  template<typename T>
  double LocalTaylorInterpolator<T>::getRadius()
  {
    return m_radius;
  }
  template<typename T>
  SearchScheme LocalTaylorInterpolator<T>::getSearchScheme()
  {
    return m_searchScheme;
  }
  template<typename T>
  void LocalTaylorInterpolator<T>::setName(std::string t_name)
  {
    m_name = t_name;
  }
  template<typename T>
  void LocalTaylorInterpolator<T>::setK(size_t t_k)
  {
    m_k = t_k;
  }
  template<typename T>
  void LocalTaylorInterpolator<T>::setN(size_t t_n)
  {
    m_n = t_n;
  }
  template<typename T>
  void LocalTaylorInterpolator<T>::setRadius(double t_radius)
  {
    m_radius = t_radius;
  }
  template<typename T>
  void LocalTaylorInterpolator<T>::setSearchScheme(std::string t_searchScheme)
  {
    m_searchScheme = SearchSchemeMap[t_searchScheme];
  }

  template<typename T>
  Matrix<T>
  LocalTaylorInterpolator<T>::constructLocalTaylorMatrix(const size_t t_index)
  {
    //  Check that t_index is not outside the bounds of ugrid.
    if (t_index >= m_ugrid->getN()) {
      m_log->ERROR("Tried to construct LTI Matrix for point "
                   + std::to_string(t_index) + " which is out of bounds for "
                   + " UGrid array of size " + std::to_string(m_ugrid->getN()));
      m_log->TRACE("Returning zero matrix.");
      return Matrix<T>("zeroes",1,0.0);
    }
    //  Determine which scheme to use
    if (m_searchScheme == SearchScheme::NEAREST_NEIGHBORS) {
      m_ugrid->queryNeighbors(m_k);
    }
    else {
      m_ugrid->queryRadius(m_radius);
    }
    //  Grab the neighbors for the point at t_index.
    std::vector<size_t> neighbors = m_ugrid->getNeighbors(t_index);
    //  Get the value of the point at t_index.
    std::vector<T> p = m_ugrid->getPoint(t_index);
    //  Construct the empty (k x m) array for matrix LTI
    std::vector<std::vector<double>> LTI(m_k,m_monomial.getTotalNumOfElements);
    for (auto i = 0; i < neighbors.size(); i++) {
      auto id = neighbors[i];
      std::vector<T> x = m_ugrid->getPoint(id);
      std::vector<double> temp = m_monomial.taylorMonomialExpansion(p,x);
      LTI.push_back(temp);
    }
    return Matrix<T>("LTI",LTI);
  }

  template<typename T>
  Matrix<T>
  LocalTaylorInterpolator<T>::constructLocalTaylorMatrix(const std::vector<T> t_point)
  {

  }

  template<typename T>
  Matrix<T>
  LocalTaylorInterpolator<T>::constructLocalTaylorMatrix(const size_t t_index,
                                                         const size_t t_k)
  {

  }

  template<typename T>
  Matrix<T>
  LocalTaylorInterpolator<T>::constructLocalTaylorMatrix(const std::vector<T>& t_point
                                                         const size_t t_k)
  {

  }

  template<typename T>
  Matrix<T>
  LocalTaylorInterpolator<T>::constructLocalTaylorMatrix(const size_t t_index,
                                                         const size_t t_k,
                                                         const size_t t_n)
  {

  }

  template<typename T>
  Matrix<T>
  LocalTaylorInterpolator<T>::constructLocalTaylorMatrix(const std::vector<T>& t_point,
                                                         const size_t t_k,
                                                         const size_t t_n)
  {

  }

}
