//------------------------------------------------------------------------------
//  local_taylor_interpolation.cpp
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
#include "local_taylor_interpolation.h"

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
  : Interpolator<T>()
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
  (std::string t_name) : Interpolator<T>(t_name)
  {
  }
  //----------------------------------------------------------------------------
  template<typename T>
  LocalTaylorInterpolator<T>::LocalTaylorInterpolator
  (std::shared_ptr<Grid<T>> t_ugrid) : Interpolator<T>(t_ugrid)
  {
    m_monomial.setDim(this->m_Grid->getDim());
  }
  //----------------------------------------------------------------------------
  template<typename T>
  LocalTaylorInterpolator<T>::LocalTaylorInterpolator
  (std::string t_name, std::shared_ptr<Grid<T>> t_ugrid)
  : Interpolator<T>(t_name, t_ugrid)
  {
    m_monomial.setDim(this->m_Grid->getDim());
  }
  //----------------------------------------------------------------------------
  template<typename T>
  LocalTaylorInterpolator<T>::LocalTaylorInterpolator
  (std::shared_ptr<Log> t_log) : Interpolator<T>(t_log)
  {
  }
  //----------------------------------------------------------------------------
  template<typename T>
  LocalTaylorInterpolator<T>::LocalTaylorInterpolator
  (std::string t_name, std::shared_ptr<Log> t_log)
  : Interpolator<T>(t_name, t_log)
  {
  }
  //----------------------------------------------------------------------------
  template<typename T>
  LocalTaylorInterpolator<T>::LocalTaylorInterpolator
  (std::shared_ptr<Grid<T>> t_ugrid, std::shared_ptr<Log> t_log)
  : Interpolator<T>(t_ugrid, t_log)
  {
    m_monomial.setDim(this->m_Grid->getDim());
  }
  //----------------------------------------------------------------------------
  template<typename T>
  LocalTaylorInterpolator<T>::LocalTaylorInterpolator
  (std::string t_name, std::shared_ptr<Grid<T>> t_ugrid,
   std::shared_ptr<Log> t_log)
  : Interpolator<T>(t_name, t_ugrid, t_log)
  {
    m_monomial.setDim(this->m_Grid->getDim());
  }
  //----------------------------------------------------------------------------
  template<typename T>
  LocalTaylorInterpolator<T>::LocalTaylorInterpolator(std::shared_ptr<Grid<T>> t_Grid,
                                std::shared_ptr<Field<T>> t_Field,
                                std::shared_ptr<Log> t_log)
  : Interpolator<T>(t_Grid,t_Field,t_log)
  {
    this->m_lsdriver = LSDriver::xGELS;
    this->m_log->TRACE(NAME() + "LocalTaylorInterpolator '" + this->m_name + "' created at location "
                + address_to_string(*this));
    this->m_log->INFO(NAME() + "Log passed to LocalTaylorInterpolator '" + this->m_name + "'");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  LocalTaylorInterpolator<T>::LocalTaylorInterpolator(std::string t_name,
                                std::shared_ptr<Grid<T>> t_Grid,
                                std::shared_ptr<Field<T>> t_Field,
                                std::shared_ptr<Log> t_log)
  : Interpolator<T>(t_name,t_Grid,t_Field,t_log)
  {
    this->m_lsdriver = LSDriver::xGELS;
    this->m_log->TRACE(NAME() + "LocalTaylorInterpolator '" + this->m_name + "' created at location "
                + address_to_string(*this));
    this->m_log->INFO(NAME() + "Log passed to LocalTaylorInterpolator '" + this->m_name + "'");
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
    this->m_log->TRACE(NAME() + "Setting m_k to " + std::to_string(t_k));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void LocalTaylorInterpolator<T>::setN(size_t t_n)
  {
    m_n = t_n;
    this->m_log->TRACE(NAME() + "Setting m_n to " + std::to_string(t_n));
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
    if (t_index >= this->m_Grid->getN()) {
      this->m_log->ERROR(NAME() + "Tried to construct LTI Matrix for point "
                   + std::to_string(t_index) + " which is out of bounds for "
                   + " Grid array of size " + std::to_string(this->m_Grid->getN()));
      this->m_log->TRACE(NAME() + "Returning zero matrix.");
      return Matrix<T>("zeroes",1,0.0);
    }
    //  Determine which scheme to use
    if (m_searchScheme == SearchScheme::NEAREST_NEIGHBORS) {
      this->m_Grid->getKDTree()->queryNeighbors(m_k);
    }
    else {
      this->m_Grid->getKDTree()->queryNeighbors(m_radius);
    }
    //  Reset the monomial to use the correct m_k and m_n values
    m_monomial.setDeg(m_n);
    //  Grab the neighbors for the point at t_index.
    std::vector<size_t>
    neighbors = this->m_Grid->getKDTree()->getCurrentNeighborIndices(t_index);
    //  Get the value of the point at t_index.
    std::vector<T> p = this->m_Grid->getPoint(t_index);
    //  Construct the empty (k x m) array for matrix LTI
    std::vector<std::vector<double>> LTI(neighbors.size());
    for (auto i = 0; i < neighbors.size(); i++) {
      auto id = neighbors[i];
      std::vector<T> x = this->m_Grid->getPoint(id);
      std::vector<double> temp = m_monomial.taylorMonomialExpansion(p,x);
      LTI[i] = std::move(temp);
    }
    return Matrix<T>("LTI",LTI);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Matrix<T>
  LocalTaylorInterpolator<T>::constructLocalTaylorMatrix(
                              const std::vector<T>& t_point)
  {
    std::vector<size_t> neighbors;
    //  Determine which scheme to use
    if (m_searchScheme == SearchScheme::NEAREST_NEIGHBORS) {
      neighbors = this->m_Grid->getKDTree()->queryNeighbors(t_point,m_k);
    }
    else {
      neighbors = this->m_Grid->getKDTree()->queryNeighbors(t_point,m_radius);
    }
    //  Reset the monomial to use the correct m_k and m_n values
    m_monomial.setDeg(m_n);
    //  Construct the empty (n x m) array for matrix LTI
    std::vector<std::vector<double>>
    LTI(neighbors.size());
    for (auto i = 0; i < neighbors.size(); i++) {
      auto id = neighbors[i];
      std::vector<T> x = this->m_Grid->getPoint(id);
      std::vector<double> temp = m_monomial.taylorMonomialExpansion(t_point,x);
      LTI[i] = std::move(temp);
    }
    return Matrix<T>("LTI",LTI);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Matrix<T>
  LocalTaylorInterpolator<T>::constructLocalTaylorMatrix(const size_t t_index,
                                                         const size_t t_k)
  {
    //  Check that t_k is not greater than the number of points
    if (t_k > this->m_Grid->getN()) {
      this->m_log->ERROR(NAME() + "Tried to construct LTI Matrix for point "
                   + std::to_string(t_index) + " using " + std::to_string(t_k)
                   + " neighbors, however there are only "
                   + std::to_string(this->m_Grid->getN()) + " total points.");
      this->m_log->TRACE(NAME() + "Setting m_k = this->m_Grid->getN()");
      m_k = this->m_Grid->getN();
    }
    else {
      m_k = t_k;
      this->m_log->TRACE(NAME() + "Setting m_k to " + std::to_string(t_k));
    }
    return constructLocalTaylorMatrix(t_index);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Matrix<T>
  LocalTaylorInterpolator<T>::constructLocalTaylorMatrix(
                              const std::vector<T>& t_point,
                              const size_t t_k)
  {
    //  Check that t_k is not greater than the number of points
    if (t_k > this->m_Grid->getN()) {
      this->m_log->ERROR(NAME() + "Tried to construct LTI Matrix for a point using "
                         + std::to_string(t_k) + " neighbors,"
                         + " however there are only "
                         + std::to_string(this->m_Grid->getN())
                         + " total points.");
      this->m_log->TRACE(NAME() + "Setting m_k = this->m_Grid->getN()");
      m_k = this->m_Grid->getN();
    }
    else {
      m_k = t_k;
      this->m_log->TRACE(NAME() + "Setting m_k to " + std::to_string(t_k));
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
      this->m_log->ERROR(NAME() + "Tried to construct LTI Matrix for point "
                   + std::to_string(t_index) + " using a influence radius of "
                   + std::to_string(t_radius));
    }
    return constructLocalTaylorMatrix(t_index);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Matrix<T>
  LocalTaylorInterpolator<T>::constructLocalTaylorMatrix(
                              const std::vector<T>& t_point,
                              const double t_radius)
  {
    //  Check that t_radius is positive
    if (t_radius < 0.0) {
      this->m_log->ERROR(NAME() + "Tried to construct LTI Matrix for point using an influence radius of " + std::to_string(t_radius));
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
      this->m_log->ERROR(NAME() + "Tried to set m_n to zero.");
    }
    else if (m_n < t_n) {
      //  Set m_n to t_n
      m_n = t_n;
      this->m_log->TRACE(NAME() + "Setting m_n to " + std::to_string(t_n));
    }
    //  Create new monomial
    m_monomial.generateMonomial(m_n);
    //  Check that t_k is not greater than the number of points
    if (t_k > this->m_Grid->getN()) {
      this->m_log->ERROR(NAME() + "Tried to construct LTI Matrix for point "
                   + std::to_string(t_index) + " using " + std::to_string(t_k)
                   + " neighbors, however there are only "
                   + std::to_string(this->m_Grid->getN()) + " total points.");
      this->m_log->TRACE(NAME() + "Setting m_k = this->m_Grid->getN()");
      m_k = this->m_Grid->getN();
    }
    else {
      m_k = t_k;
      this->m_log->TRACE(NAME() + "Setting m_k to " + std::to_string(t_k));
    }
    return constructLocalTaylorMatrix(t_index);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Matrix<T>
  LocalTaylorInterpolator<T>::constructLocalTaylorMatrix(
                              const std::vector<T>& t_point,
                              const size_t t_k,
                              const size_t t_n)
  {
    //  Check that t_n is at least 1
    if (t_n == 0) {
      this->m_log->ERROR(NAME() + "Tried to set m_n to zero.");
    }
    else if (m_n < t_n) {
      //  Set m_n to t_n
      m_n = t_n;
      this->m_log->TRACE(NAME() + "Setting m_n to " + std::to_string(t_n));
    }
    //  Create new monomial
    m_monomial.generateMonomial(m_n);
    //  Check that t_k is not greater than the number of points
    if (t_k > this->m_Grid->getN()) {
      this->m_log->ERROR(NAME() + "Tried to construct LTI Matrix for a point using "
                         + std::to_string(t_k) + " neighbors,"
                         + " however there are only "
                         + std::to_string(this->m_Grid->getN())
                         + " total points.");
      this->m_log->TRACE(NAME() + "Setting m_k = this->m_Grid->getN()");
      m_k = this->m_Grid->getN();
    }
    else {
      m_k = t_k;
      this->m_log->TRACE(NAME() + "Setting m_k to " + std::to_string(t_k));
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
      this->m_log->ERROR(NAME() + "Tried to set m_n to zero.");
    }
    else if (m_n < t_n) {
      //  Set m_n to t_n
      m_n = t_n;
      this->m_log->TRACE(NAME() + "Setting m_n to " + std::to_string(t_n));
    }
    //  Create new monomial
    m_monomial.generateMonomial(m_n);
    //  Check that t_radius is positive
    if (t_radius < 0.0) {
      this->m_log->ERROR(NAME() + "Tried to construct LTI Matrix for point "
                   + std::to_string(t_index) + " using a influence radius of "
                   + std::to_string(t_radius));
    }
    else {
      m_radius = t_radius;
      this->m_log->TRACE(NAME() + "Setting m_radius to " + std::to_string(t_radius));
    }
    return constructLocalTaylorMatrix(t_index);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Matrix<T>
  LocalTaylorInterpolator<T>::constructLocalTaylorMatrix(
                              const std::vector<T>& t_point,
                              const double t_radius,
                              const size_t t_n)
  {
    //  Check that t_n is at least 1
    if (t_n == 0) {
      this->m_log->ERROR(NAME() + "Tried to set m_n to zero.");
    }
    else if (m_n < t_n) {
      //  Set m_n to t_n
      m_n = t_n;
      this->m_log->TRACE(NAME() + "Setting m_n to " + std::to_string(t_n));
    }
    //  Create new monomial
    m_monomial.generateMonomial(m_n);
    //  Check that t_radius is positive
    if (t_radius < 0.0) {
      this->m_log->ERROR(NAME() + "Tried to construct LTI Matrix for point using an influence radius of " + std::to_string(t_radius));
    }
    else {
      m_radius = t_radius;
      this->m_log->TRACE(NAME() + "Setting m_radius to " + std::to_string(t_radius));
    }
    return constructLocalTaylorMatrix(t_point);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Vector<T>
  LocalTaylorInterpolator<T>::derivative(const size_t t_index,
                                         const size_t t_degree)
  {
    //  Check the type of the field
    if (this->m_Field->getType() == FieldType::SCALAR) {
      return scalarFieldDerivative(t_index, t_degree);
    }
  }
  //----------------------------------------------------------------------------
  template<typename T>
  T LocalTaylorInterpolator<T>::derivative(const size_t t_index,
                                           const size_t t_degree,
                                           const size_t t_direction)
  {
    //  Check the type of the field
    if (this->m_Field->getType() == FieldType::SCALAR) {
      return scalarFieldDerivative(t_index, t_degree, t_direction);
    }
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Vector<T>
  LocalTaylorInterpolator<T>::derivative(const std::vector<T>& t_point,
                                         const size_t t_degree)
  {
    //  Check the type of the field
    if (this->m_Field->getType() == FieldType::SCALAR) {
      return scalarFieldDerivative(t_point, t_degree);
    }
  }
  //----------------------------------------------------------------------------
  template<typename T>
  T LocalTaylorInterpolator<T>::derivative(const std::vector<T>& t_point,
                                           const size_t t_degree,
                                           const size_t t_direction)
  {
    //  Check the type of the field
    if (this->m_Field->getType() == FieldType::SCALAR) {
      return scalarFieldDerivative(t_point, t_degree, t_direction);
    }
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<Vector<T>>
  LocalTaylorInterpolator<T>::fieldDerivative(const size_t t_degree)
  {
    std::vector<Vector<T>> result(this->m_Grid->getN());
    for (auto i = 0; i < result.size(); i++) {
      result[i] = derivative(i,t_degree);
    }
    return result;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<T>
  LocalTaylorInterpolator<T>::fieldDerivative(const size_t t_degree,
                                              const size_t t_direction)
  {
    std::vector<T> result(this->m_Grid->getN());
    for (auto i = 0; i < result.size(); i++) {
      result[i] = derivative(i,t_degree,t_direction);
    }
    return result;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Vector<T>
  LocalTaylorInterpolator<T>::scalarFieldDerivative(const size_t t_index,
                                                    const size_t t_degree)
  {
    Vector<T> result(this->m_Grid->getDim());
    //  Construct the LTI Matrix for the point at index.
    Matrix<T> LTI = constructLocalTaylorMatrix(t_index, m_k, t_degree);
    //  Construct the field object
    auto f = this->m_Field->constructLocalFieldValues(t_index);
    //  Solve the system
    auto s = this->xGELSx(LTI,f);
    //  Grab each value corresponding to the derivative
    for (size_t j = 0; j < this->m_Grid->getDim(); j++)
		{
			std::vector<size_t> deriv(this->m_Grid->getDim(),0);
			deriv[j] = t_degree;
			size_t l = m_monomial.getTaylorIndex(deriv);
			result(j) = factorial(t_degree)*s(l);
		}
    return result;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  T LocalTaylorInterpolator<T>::scalarFieldDerivative(const size_t t_index,
                                                      const size_t t_degree,
                                                      const size_t t_direction)
  {
    return scalarFieldDerivative(t_index,t_degree)(t_direction);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Vector<T>
  LocalTaylorInterpolator<T>::scalarFieldDerivative(
                              const std::vector<T>& t_point,
                              const size_t t_degree)
  {
    Vector<T> result(this->m_Grid->getDim());
    //  Construct the LTI Matrix for the point at index.
    Matrix<T> LTI = constructLocalTaylorMatrix(t_point, m_k, t_degree);
    //  Construct the field object
    auto f = this->m_Field->constructLocalFieldValues(t_point,m_k);
    //  Solve the system
    auto s = this->xGELSx(LTI,f);
    //  Grab each value corresponding to the derivative
    for (size_t j = 0; j < this->m_Grid->getDim(); j++)
		{
			std::vector<size_t> deriv(this->m_Grid->getDim(),0);
			deriv[j] = t_degree;
			size_t l = m_monomial.getTaylorIndex(deriv);
			result(j) = factorial(t_degree)*s(l);
		}
    return result;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  T LocalTaylorInterpolator<T>::scalarFieldDerivative(
                                const std::vector<T>& t_point,
                                const size_t t_degree,
                                const size_t t_direction)
  {
    return scalarFieldDerivative(t_point,t_degree)(t_direction);
  }
  //----------------------------------------------------------------------------
}
