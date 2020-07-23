//------------------------------------------------------------------------------
//  local_taylor_interpolant.h
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

namespace ET
{
  //! \class Local Taylor Interpolant
  /*! Local Taylor Interpolant class.
   */
  template<typename T>
  class LocalTaylorInterpolant : public Interpolant<T>
  {
  public:
    LocalTaylorInterpolant();
    ~LocalTaylorInterpolant();

    //  Operator overloads
    /*! Access operator.  Operator for interpolating at a point.
     *  @param t_point.  A std::vector<T> of the point to interpolate at.
     *  @return The interpolation at t_point.
    */
    T& operator()(const std::vector<T>& t_point);


  private:
    /*! Expansion point.  The point that the LTI is expanded around.
     */
    std::vector<T> m_expansion_point {0};
    /*! Expansion coefficients.  The coefficients of the Taylor expansion.
     */
    std::vector<T> m_expansion_coeffs {0};

  };

  //----------------------------------------------------------------------------
  //  Global interpolation functions
  //----------------------------------------------------------------------------
}
