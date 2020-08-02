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
  template<typename T>
  T LocalTaylorInterpolant<T>::d(const std::vector<T>& t_point)
  {
    //  Construct the monomial expansion for the point with respect to the
    //  expansion point.
    std::vector<double>
    temp = m_monomial.taylorMonomialExpansion(m_expansion_point,t_point);
    T result = 0;
    for (auto i = 0; i < temp.size()-1; i++) {
      result += T(i+1) * temp[i] * m_expansion_coefficients[i+1];
    }
    return result;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<T>
  LocalTaylorInterpolant<T>::d(const std::vector<std::vector<T>>& t_points)
  {
    //  Construct the monomial expansion for the point with respect to the
    //  expansion point.
    std::vector<T> result(t_points.size(),0.0);
    for (auto i = 0; i < result.size(); i++) {
      std::vector<double>
      temp = m_monomial.taylorMonomialExpansion(m_expansion_point,t_points[i]);
      for (auto j = 0; j < temp.size()-1; j++) {
        result[i] += T(j+1) * temp[j] * m_expansion_coefficients[j+1];
      }
    }
    return result;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  T LocalTaylorInterpolant<T>::dd(const std::vector<T>& t_point)
  {
    //  Construct the monomial expansion for the point with respect to the
    //  expansion point.
    std::vector<double>
    temp = m_monomial.taylorMonomialExpansion(m_expansion_point,t_point);
    T result = 0;
    for (auto i = 0; i < temp.size()-2; i++) {
      result += T(i+2)*(i+1) * temp[i] * m_expansion_coefficients[i+2];
    }
    return result;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<T>
  LocalTaylorInterpolant<T>::dd(const std::vector<std::vector<T>>& t_points)
  {
    //  Construct the monomial expansion for the point with respect to the
    //  expansion point.
    std::vector<T> result(t_points.size(),0.0);
    for (auto i = 0; i < result.size(); i++) {
      std::vector<double>
      temp = m_monomial.taylorMonomialExpansion(m_expansion_point,t_points[i]);
      for (auto j = 0; j < temp.size()-2; j++) {
        result[i] += T(j+2)*(j+1) * temp[j] * m_expansion_coefficients[j+2];
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
  createLocalTaylorInterpolant(Field<T>& t_Field,
                               const std::vector<T>& t_expansion_point,
                               const size_t& t_n)
  {
    const size_t num_points = t_Field.getGrid()->getN();
    //  Construct a local interpolator
    LocalTaylorInterpolator<T> interp(t_Field.getGrid());

    //  Create the local interpolation matrix using every point
    Matrix<T>
    LTI = interp.constructLocalTaylorMatrix(t_expansion_point,
                                            num_points,
                                            t_n);
    auto f = t_Field.constructLocalFieldValues(t_expansion_point,
                                               num_points);
    //  Incorporate boundary conditions
    std::vector<BoundaryCondition<T>>
    conditions = t_Field.getBoundaryConditions();
    //  TODO: What strategy should we use here?
    //        Either we remove the boundary points from LTI
    //        Or we have boundary points be separate from the Grid.
    //        Either way, we need a local taylor matrix for the boundary
    //        points, A_ij.
    //  NOTE: Going with removing rows for now to construct a new matrix.
    //        The goal here is to remove the rows corresponding to the
    //        boundary points, and then construct the constraint matrix.
    Matrix<T>
    C = interp.constructBoundaryConditionMatrix(t_Field,t_expansion_point);
    Vector<T> s;
    Vector<T> f_C;
    if (conditions.size() == 0) {
      //  Solve the system
      s = interp.xGELSx(LTI,f);
    }
    else {
      //  Create larger matrix [[A^TA C^T],[C 0]]
      Matrix<T> LTI_T = LTI.transpose();
      Matrix<T> A_TA = LTI_T * LTI;
      Matrix<T> C_T = C.transpose();
      Matrix<T> B = A_TA;
      Matrix<T> B2 = C_T;
      for (auto i = 0; i < C.getNumRows(); i++) {
        B.addRow(B.getNumRows(),C.getRow(i));
        std::vector<T> z(C.getNumCols(),0.0);
        B2.addRow(B.getNumRows(),z);
      }
      for (auto i = 0; i < B2.getNumCols(); i++) {
        B.addCol(B.getNumCols(),B2.getCol(i));
      }
      //  Create solution vector [A^Tf f_C]
      auto new_f = LTI_T * f;
      for (auto i = 0; i < f_C.getDim(); i++) {
        new_f.addVal(new_f.getDim(),f_C(i));
      }
      //  find the solution vector
      s = interp.xGELSx(B,new_f);
      //  reduce solution to the coefficients of the expansion
      for (auto i = 0; i < f_C.getDim(); i++) {
        s.removeVal(s.getDim());
      }
    }
    //  Solve the system
    //auto s = interp.xGELSx(LTI,f);
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
  (Field<double>& t_Field,
   const std::vector<double>& t_expansion_point, const size_t& t_n);
}
