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

#include <vector>
#include <iostream>
#include <map>
#include <memory>

#include "ugrid.h"
#include "geometry.h"
#include "params.h"
#include "utilities.h"
#include "vector.h"
#include "matrix.h"
#include "field.h"
#include "interpolator.h"
#include "local_taylor_interpolation.h"

#include "interpolant.h"

namespace ET
{
  //! \class Local Taylor Interpolant
  /*! Local Taylor Interpolant class.
   */
  template<typename T>
  class LocalTaylorInterpolant : public Interpolant<T>
  {
  public:
    //! Default Constructor
    /*! Default constructor for an instance of the base class.
     */
    LocalTaylorInterpolant();
    //! Default Destructor
    /*! Default destructor for an instance of the base class.
     */
    ~LocalTaylorInterpolant();

    //  Getters and Setters
    //! Get n.
    /*! Get n.  Get the order of the Taylor expansion.
     *  @return A size_t object denoting the order of the Taylor expansion.
     */
    size_t get_n() const;
    //! Get expansion point
    /*! Get the point that the Taylor expansion is defined around.
     *  @return A std::vector<T> quantifying the expansion point.
     */
    std::vector<T> getExpansionPoint() const;
    //! Get expansion coefficients
    /*! Get the expansion coefficients for the expansion.
     *  @return A std:vector<T> of the expansion coefficients.
     */
    std::vector<T> getExpansionCoefficients() const;
    //! Get monomial
    /*! Get the monomial object.
     *  @return A monomial object.
     */
    Monomial getMonomial() const;
    //! Set n.
    /*! Set n.  Set the order of the Taylor expansion.
     *  @param t_n A size_t object denoting the order of the expansion.
     */
    void set_n(const size_t t_n);
    //! Set expansion point
    /*! Set the point that the expansion is defined around.
     *  @param t_point A std::vector<T> quantifying the expansion point.
     */
    void setExpansionPoint(const std::vector<T>& t_point);
    //! Set expansion coefficients
    /*! Set the expansion coefficients for the Taylor expansion.
     *  @param t_coefficients A std::vector<T> quantifying the expansion.
     */
    void setExpansionCoefficients(const std::vector<T>& t_coefficients);
    //! Set monomial
    /*! Sets the monomial object.
     *  @param t_monomial A Monomial object.
     */
    void setMonomial(const Monomial t_monomial);
    //  Operator overloads
    /*! Access operator.  Operator for interpolating at a point.
     *  @param t_point.  A std::vector<T> of the point to interpolate at.
     *  @return The interpolation at t_point.
    */
    virtual T operator()(const std::vector<T>& t_point);


  private:
    /*! Order of the expansion.  Default is set to n=3.
     */
    size_t m_n {3};
    /*! Expansion point.  The point that the LTI is expanded around.
     */
    std::vector<T> m_expansion_point {0};
    /*! Expansion coefficients.  The coefficients of the Taylor expansion.
     */
    std::vector<T> m_expansion_coefficients {0};
    /*! Monomial.  An instance of the monomial class for constructing
     *  Taylor monomials.  Defaulted to use dim=1 and n=3.
     */
    Monomial m_monomial {1,3};
    /*! Logging system name generator.
     */
    virtual std::string NAME() const {
      return "LTInterpolant:[" + std::to_string(m_id.id) + "]:"
             + this->m_name + ":";
    }
    /*! Unique ID for each instance.
     */
    UniqueID m_id;

  };

  //----------------------------------------------------------------------------
  //  Global interpolation functions
  //----------------------------------------------------------------------------
  template<typename T>
  LocalTaylorInterpolant<T>
  createLocalTaylorInterpolant(const std::shared_ptr<Field<T>> t_field,
                               const std::vector<T>& t_expansion_point,
                               const size_t& t_n);
  //----------------------------------------------------------------------------

  template class LocalTaylorInterpolant<double>;

}
