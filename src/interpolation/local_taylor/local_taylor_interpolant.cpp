//------------------------------------------------------------------------------
//  local_taylor_interpolant.cpp
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
#include "local_taylor_interpolant.h"

namespace ET
{
  //----------------------------------------------------------------------------
  template<typename T>
  LocalTaylorInterpolant<T>::LocalTaylorInterpolant()
  : Interpolant<T>()
  {
  }
  //----------------------------------------------------------------------------
  template<typename T>
  LocalTaylorInterpolant<T>::~LocalTaylorInterpolant()
  {
  }
  //----------------------------------------------------------------------------
  template<typename T>
  size_t LocalTaylorInterpolant<T>::get_n() const
  {
    return m_n;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<T> LocalTaylorInterpolant<T>::getExpansionPoint() const
  {
    return m_expansion_point;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<T> LocalTaylorInterpolant<T>::getExpansionCoefficients() const
  {
    return m_expansion_coefficients;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Monomial LocalTaylorInterpolant<T>::getMonomial() const
  {
    return m_monomial;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void LocalTaylorInterpolant<T>::set_n(const size_t t_n)
  {
    m_n = t_n;
    m_monomial.setDeg(m_n);
    m_monomial.generateMonomial();
    this->m_log->INFO(NAME() + "Set 'n' to " + std::to_string(t_n));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void
  LocalTaylorInterpolant<T>::setExpansionPoint(const std::vector<T>& t_point)
  {
    m_expansion_point = t_point;
    this->m_log->INFO(NAME() + "Set expansion point.");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void
  LocalTaylorInterpolant<T>::setExpansionCoefficients
  (const std::vector<T>& t_coefficients)
  {
    m_expansion_coefficients = t_coefficients;
    this->m_log->INFO(NAME() + "Set expansion coefficients.");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void LocalTaylorInterpolant<T>::setMonomial(const Monomial t_monomial)
  {
    m_monomial = t_monomial;
    this->m_log->INFO(NAME() + "Set Monomial object.");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  T LocalTaylorInterpolant<T>::operator()(const std::vector<T>& t_point)
  {
    //  Construct the monomial expansion for the point with respect to the
    //  expansion point.
    std::vector<double>
    temp = m_monomial.taylorMonomialExpansion(m_expansion_point,t_point);
    T result = 0;
    for (auto i = 0; i < temp.size(); i++) {
      result += temp[i] * m_expansion_coefficients[i];
    }
    return result;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<T>
  LocalTaylorInterpolant<T>::operator()(const std::vector<std::vector<T>>& t_points)
  {
    //  Construct the monomial expansion for the point with respect to the
    //  expansion point.
    std::vector<T> result(t_points.size(),0.0);
    for (auto i = 0; i < result.size(); i++) {
      std::vector<double>
      temp = m_monomial.taylorMonomialExpansion(m_expansion_point,t_points[i]);
      for (auto j = 0; j < temp.size(); j++) {
        result[i] += temp[j] * m_expansion_coefficients[j];
      }
    }
    return result;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Global interpolation functions
  //----------------------------------------------------------------------------
  template<typename T>
  LocalTaylorInterpolant<T>
  createLocalTaylorInterpolant(const std::shared_ptr<Field<T>> t_field,
                               const std::vector<T>& t_expansion_point,
                               const size_t& t_n)
  {
    const size_t num_points = t_field->getGrid()->getN();
    //  Construct a local interpolator
    LocalTaylorInterpolator<T> interp(t_field->getGrid());
    interp.setField(t_field);

    //  Create the local interpolation matrix using every point
    Matrix<T>
    LTI = interp.constructLocalTaylorMatrix(t_expansion_point,
                                            num_points,
                                            t_n);
    auto f = interp.getField()->constructLocalFieldValues(t_expansion_point,
                                                       num_points);
    //  Solve the system
    auto s = interp.xGELSx(LTI,f);
    //  Create the interpolant
    LocalTaylorInterpolant<T> interpolant;
    interpolant.set_n(t_n);
    interpolant.setExpansionPoint(t_expansion_point);
    interpolant.setExpansionCoefficients(s.getVec());
    return interpolant;
  }
  //----------------------------------------------------------------------------

  template LocalTaylorInterpolant<double>
  createLocalTaylorInterpolant<double>
  (const std::shared_ptr<Field<double>> t_field,
   const std::vector<double>& t_expansion_point, const size_t& t_n);
}
