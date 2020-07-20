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
  //----------------------------------------------------------------------------
  template<typename T>
  LocalTaylorInterpolator<T>::LocalTaylorInterpolator()
  : Interpolator(), m_name("default")
  {
  }
  //----------------------------------------------------------------------------
  template<typename T>
  LocalTaylorInterpolator<T>::~LocalTaylorInterpolator()
  {
  }
  //----------------------------------------------------------------------------
  template<typename T>
  LocalTaylorInterpolator<T>::LocalTaylorInterpolator
  (std::string t_name) : Interpolator(t_name)
  {
  }
  //----------------------------------------------------------------------------
  template<typename T>
  LocalTaylorInterpolator<T>::LocalTaylorInterpolator
  (std::shared_ptr<UGrid<T>> t_ugrid) : Interpolator(t_ugrid)
  {
  }
  //----------------------------------------------------------------------------
  template<typename T>
  LocalTaylorInterpolator<T>::LocalTaylorInterpolator
  (std::string t_name, std::shared_ptr<UGrid<T>> t_ugrid)
  : Interpolator(t_name, t_ugrid)
  {
  }
  //----------------------------------------------------------------------------
  template<typename T>
  LocalTaylorInterpolator<T>::LocalTaylorInterpolator
  (std::shared_ptr<Log> t_log) : Interpolator(t_log)
  {
  }
  //----------------------------------------------------------------------------
  template<typename T>
  LocalTaylorInterpolator<T>::LocalTaylorInterpolator
  (std::string t_name, std::shared_ptr<Log> t_log)
  : Interpolator(t_name, t_log)
  {
  }
  //----------------------------------------------------------------------------
  template<typename T>
  LocalTaylorInterpolator<T>::LocalTaylorInterpolator
  (std::shared_ptr<UGrid<T>> t_ugrid, std::shared_ptr<Log> t_log)
  : Interpolator(t_ugrid, t_log)
  {
  }
  //----------------------------------------------------------------------------
  template<typename T>
  LocalTaylorInterpolator<T>::LocalTaylorInterpolator
  (std::string t_name, std::shared_ptr<UGrid<T>> t_ugrid,
   std::shared_ptr<Log> t_log)
  : Interpolator(t_name, t_ugrid, t_log)
  {
  }
  //----------------------------------------------------------------------------
  template<typename T>
  size_t LocalTaylorInterpolator<T>::getK()
  {
    return m_k;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  size_t LocalTaylorInterpolator<T>::getN()
  {
    return m_n;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  double LocalTaylorInterpolator<T>::getRadius()
  {
    return m_radius;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  SearchScheme LocalTaylorInterpolator<T>::getSearchScheme()
  {
    return m_searchScheme;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void LocalTaylorInterpolator<T>::setK(size_t t_k)
  {
    m_k = t_k;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void LocalTaylorInterpolator<T>::setN(size_t t_n)
  {
    m_n = t_n;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void LocalTaylorInterpolator<T>::setRadius(double t_radius)
  {
    m_radius = t_radius;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void LocalTaylorInterpolator<T>::setSearchScheme(std::string t_searchScheme)
  {
    m_searchScheme = SearchSchemeMap[t_searchScheme];
  }
  //----------------------------------------------------------------------------
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
  //----------------------------------------------------------------------------
  template<typename T>
  Matrix<T>
  LocalTaylorInterpolator<T>::constructLocalTaylorMatrix(const std::vector<T> t_point)
  {
    std::vector<size_t> neighbors;
    //  Determine which scheme to use
    if (m_searchScheme == SearchScheme::NEAREST_NEIGHBORS) {
      neighbors = m_ugrid->queryNeighbors(t_point,m_k);
    }
    else {
      neighbors = m_ugrid->queryNeighbors(t_point,m_radius);
    }
    //  Construct the empty (n x m) array for matrix LTI
    std::vector<std::vector<double>> LTI(neighbors.size(),
                                         m_monomial.getTotalNumOfElements);
    for (auto i = 0; i < neighbors.size(); i++) {
      auto id = neighbors[i];
      std::vector<T> x = m_ugrid->getPoint(id);
      std::vector<double> temp = m_monomial.taylorMonomialExpansion(p,x);
      LTI.push_back(temp);
    }
    return Matrix<T>("LTI",LTI);
    }
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Matrix<T>
  LocalTaylorInterpolator<T>::constructLocalTaylorMatrix(const size_t t_index,
                                                         const size_t t_k)
  {
    //  Check that t_k is not greater than the number of points
    if (t_k >= m_ugrid->getN()) {
      m_log->ERROR("Tried to construct LTI Matrix for point "
                   + std::to_string(t_index) + " using " + std::to_string(t_k)
                   + " neighbors, however there are only "
                   + std::to_string(m_ugrid->getN()) + " total points.");
      m_log->TRACE("Setting m_k = m_ugrid->getN()");
      m_k = ugrid->getN();
    }
    return constructLocalTaylorMatrix(t_index);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Matrix<T>
  LocalTaylorInterpolator<T>::constructLocalTaylorMatrix(const std::vector<T>& t_point
                                                         const size_t t_k)
  {
    //  Check that t_k is not greater than the number of points
    if (t_k >= m_ugrid->getN()) {
      m_log->ERROR("Tried to construct LTI Matrix for point "
                   + std::to_string(t_index) + " using " + std::to_string(t_k)
                   + " neighbors, however there are only "
                   + std::to_string(m_ugrid->getN()) + " total points.");
      m_log->TRACE("Setting m_k = m_ugrid->getN()");
      m_k = ugrid->getN();
    }
    return constructLocalTaylorMatrix(t_point);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Matrix<T>
  LocalTaylorInterpolator<T>::constructLocalTaylorMatrix(const size_t t_index,
                                                         const double t_radius)
  {
    //  Check that t_radius is positive
    if (t_radius < 0.0) {
      m_log->ERROR("Tried to construct LTI Matrix for point "
                   + std::to_string(t_index) + " using a influence radius of "
                   + std::to_string(t_radius));
    }
    return constructLocalTaylorMatrix(t_index);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Matrix<T>
  LocalTaylorInterpolator<T>::constructLocalTaylorMatrix(const std::vector<T>& t_point
                                                         const double t_radius)
  {
    //  Check that t_radius is positive
    if (t_radius < 0.0) {
      m_log->ERROR("Tried to construct LTI Matrix for point "
                   + std::to_string(t_index) + " using a influence radius of "
                   + std::to_string(t_radius));
    }
    return constructLocalTaylorMatrix(t_point);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Matrix<T>
  LocalTaylorInterpolator<T>::constructLocalTaylorMatrix(const size_t t_index,
                                                         const size_t t_k,
                                                         const size_t t_n)
  {
    //  Check that t_n is at least 1
    if (t_n == 0) {
      m_log->ERROR("Tried to set m_n to zero.");
    }
    else {
      //  Set m_n to t_n
      m_n = t_n;
      m_log->TRACE("Setting m_n to " + std::to_string(t_n));
    }
    //  Create new monomial
    m_monomial.generateMonomial(m_ugrid->getDim(),m_n);
    //  Check that t_k is not greater than the number of points
    if (t_k >= m_ugrid->getN()) {
      m_log->ERROR("Tried to construct LTI Matrix for point "
                   + std::to_string(t_index) + " using " + std::to_string(t_k)
                   + " neighbors, however there are only "
                   + std::to_string(m_ugrid->getN()) + " total points.");
      m_log->TRACE("Setting m_k = m_ugrid->getN()");
      m_k = ugrid->getN();
    }
    return constructLocalTaylorMatrix(t_index);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Matrix<T>
  LocalTaylorInterpolator<T>::constructLocalTaylorMatrix(const std::vector<T>& t_point,
                                                         const size_t t_k,
                                                         const size_t t_n)
  {
    //  Check that t_n is at least 1
    if (t_n == 0) {
      m_log->ERROR("Tried to set m_n to zero.");
    }
    else {
      //  Set m_n to t_n
      m_n = t_n;
      m_log->TRACE("Setting m_n to " + std::to_string(t_n));
    }
    //  Create new monomial
    m_monomial.generateMonomial(m_ugrid->getDim(),m_n);
    //  Check that t_k is not greater than the number of points
    if (t_k >= m_ugrid->getN()) {
      m_log->ERROR("Tried to construct LTI Matrix for point "
                   + std::to_string(t_index) + " using " + std::to_string(t_k)
                   + " neighbors, however there are only "
                   + std::to_string(m_ugrid->getN()) + " total points.");
      m_log->TRACE("Setting m_k = m_ugrid->getN()");
      m_k = ugrid->getN();
    }
    return constructLocalTaylorMatrix(t_point);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Matrix<T>
  LocalTaylorInterpolator<T>::constructLocalTaylorMatrix(const size_t t_index,
                                                         const double t_radius,
                                                         const size_t t_n)
  {
    //  Check that t_n is at least 1
    if (t_n == 0) {
      m_log->ERROR("Tried to set m_n to zero.");
    }
    else {
      //  Set m_n to t_n
      m_n = t_n;
      m_log->TRACE("Setting m_n to " + std::to_string(t_n));
    }
    //  Create new monomial
    m_monomial.generateMonomial(m_ugrid->getDim(),m_n);
    //  Check that t_radius is positive
    if (t_radius < 0.0) {
      m_log->ERROR("Tried to construct LTI Matrix for point "
                   + std::to_string(t_index) + " using a influence radius of "
                   + std::to_string(t_radius));
    }
    return constructLocalTaylorMatrix(t_index);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Matrix<T>
  LocalTaylorInterpolator<T>::constructLocalTaylorMatrix(const std::vector<T>& t_point,
                                                         const double t_radius,
                                                         const size_t t_n)
  {
    //  Check that t_n is at least 1
    if (t_n == 0) {
      m_log->ERROR("Tried to set m_n to zero.");
    }
    else {
      //  Set m_n to t_n
      m_n = t_n;
      m_log->TRACE("Setting m_n to " + std::to_string(t_n));
    }
    //  Create new monomial
    m_monomial.generateMonomial(m_ugrid->getDim(),m_n);
    //  Check that t_radius is positive
    if (t_radius < 0.0) {
      m_log->ERROR("Tried to construct LTI Matrix for point "
                   + std::to_string(t_index) + " using a influence radius of "
                   + std::to_string(t_radius));
    }
    return constructLocalTaylorMatrix(t_point);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<T>
  LocalTaylorInterpolator<T>::derivative(const size_t t_index,
                                         const size_t t_degree)
  {
    Vector<T> result(m_ugrid->getDim());
    //  Construct the LTI Matrix for the point at index.
    Matrix<T> LTI = constructLocalTaylorMatrix(t_index, t_degree);
    //  Construct the field object
    auto f = m_field->constructLocalFieldValues(t_index);
    //  Solve the system
    auto s = xGELSx(LTI,f);
    //  Grab each value corresponding to the derivative
    for (size_t j = 0; j < m_ugrid->getDim(); j++)
		{
			std::vector<size_t> deriv(m_ugrid->getDim(),0);
			deriv[j] = t_degree;
			size_t l = m_monomial.getTaylorIndex(deriv);
			result(j) = coefficients(l);
		}
    return result;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  T LocalTaylorInterpolator<T>::derivative(const size_t t_index,
                                           const size_t t_degree,
                                           const size_t t_direction)
  {
    return derivative(t_index,t_degree)[t_direction];
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Vector<T>
  LocalTaylorInterpolator<T>::derivative(const size_t std::vector<T>& point,
                                         const size_t t_degree)
  {
    Vector<T> result(m_ugrid->getDim());
    //  Construct the LTI Matrix for the point at index.
    Matrix<T> LTI = constructLocalTaylorMatrix(t_point, t_degree);
    //  Construct the field object
    auto f = m_field->constructLocalFieldValues(t_index);
    //  Solve the system
    auto s = xGELSx(LTI,f);
    //  Grab each value corresponding to the derivative
    for (size_t j = 0; j < m_ugrid->getDim(); j++)
		{
			std::vector<size_t> deriv(m_ugrid->getDim(),0);
			deriv[j] = t_degree;
			size_t l = m_monomial.getTaylorIndex(deriv);
			result(j) = coefficients(l);
		}
    return result;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  T LocalTaylorInterpolator<T>::derivative(const size_t std::vector<T>& point,
                                           const size_t t_degree,
                                           const size_t t_direction)
  {
  }
  //----------------------------------------------------------------------------
}
