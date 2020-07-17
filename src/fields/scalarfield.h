//------------------------------------------------------------------------------
//  scalarfield.h
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
#include <memory>

#include "ugrid.h"
#include "log.h"

//------------------------------------------------------------------------------
//  Forward declaration of Interpolator, Integrator and DiffEQ
//------------------------------------------------------------------------------
namespace ET
{
  //template<typename T> class Field;
  template<typename T> class Interpolator;
  template<typename T> class DiffEQ;
  template<typename T> class Integrator;
}

#include "field.h"
#include "interpolator.h"
#include "diffeq.h"
#include "integrator.h"

namespace ET
{
  //----------------------------------------------------------------------------
  //  Scalar fields class.
  //----------------------------------------------------------------------------
  template<typename T>
  class ScalarField: public Field<T>
  {
  public:
    //! Defualt Constructor
    /*! Default constructor for Field.
     */
    ScalarField();
    //! Destructor
    /*! Destructor for Field.
     */
    virtual ~ScalarField();
    //! Constructor
    /*! constructor for Field that takes a Logger
     */
    ScalarField(std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for Field that takes a UGrid
     */
    ScalarField(std::shared_ptr<UGrid<T>> t_ugrid);
    //! Constructor
    /*! constructor for Field that takes a Interpolator
     */
    ScalarField(std::shared_ptr<Interpolator<T>> t_interpolator);
    //! Constructor
    /*! constructor for Field that takes a UGrid and a logger
     */
    ScalarField(std::shared_ptr<UGrid<T>> t_ugrid, std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for Field that takes a field
     */
    ScalarField(std::vector<T> t_field);
    //! Constructor
    /*! constructor for Field that takes a field and aLogger
     */
    ScalarField(std::vector<T> t_field, std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for Field that takes a field and a UGrid
     */
    ScalarField(std::vector<T> t_field, std::shared_ptr<UGrid<T>> t_ugrid);
    //! Constructor
    /*! constructor for Field that takes a field and a Interpolator
     */
    ScalarField(std::vector<T> t_field,
                std::shared_ptr<Interpolator<T>> t_interpolator);
    //! Constructor
    /*! constructor for Field that takes a field, UGrid and a logger
     */
    ScalarField(std::vector<T> t_field,
                std::shared_ptr<UGrid<T>> t_ugrid,
                std::shared_ptr<Log> t_log);

    //! Get field.
    /*! get field.  Returns a const copy of the field values vector.
     *  @return m_field An std::vector<T> of field values.
     */
    std::vector<T> getField() const;
    //! Access field.
    /*! access field.  Returns a pointer to the field values vector.
     *  @return *m_field A pointer to m_field.
     */
    std::vector<T>* accessField();
    //! Data
    /*! data.  Returns a pointer to the beginning of the values vector.
     *  @return m_field.data() A pointer to the beginning of m_field.
     */
    T* data();
    //! Set field
    /*! set field.  Sets the field values vector.
     *  @param t_field An std::vector<T> of field values.
     */
    void setField(std::vector<T> t_field);
    //  Operator overloads
    //! Plus.
    /*! plus operator.  Add two scalar fields together
     *  @param t_scalar A const ScalarField<T>& reference.
     *  @return A new ScalarField<T> that is the sum of this and t_scalar.
     */
    ScalarField<T> operator+(const ScalarField<T>& t_scalar) const;
    //! Minus.
    /*! subtraction operator.  Subtract two scalar fields from each other.
     *  @param t_scalar A const ScalarField<T>& reference.
     *  @return A new ScalarField<T> that is the difference of
     *  t_scalar from this.
     */
    ScalarField<T> operator-(const ScalarField<T>& t_scalar) const;
    //! Mulitply.
    /*! multiplication operator.  Multiply two scalar fields together
     *  @param t_scalar A const ScalarField<T>& reference.
     *  @return A new ScalarField<T> that is the product of this and t_scalar.
     */
    ScalarField<T> operator*(const ScalarField<T>& t_scalar) const;
    //! Divide.
    /*! divide operator.  Divide two scalar fields.
     *  @param t_scalar A const ScalarField<T>& reference.
     *  @return A new ScalarField<T> that is the ratio of this and t_scalar.
     */
    ScalarField<T> operator/(const ScalarField<T>& t_scalar) const;
    //! Plus equals.
    /*! plus equals operator.  Add sa scalar field to this.
     *  @param t_scalar A const ScalarField<T>& reference.
     *  @return This field plus t_scalar.
     */
    ScalarField<T>& operator+=(const ScalarField<T>& t_scalar);
    //! Minus equals.
    /*! minus equals operator.  Subtract a scalar field from this.
     *  @param t_scalar A const ScalarField<T>& reference.
     *  @return This field minus t_scalar.
     */
    ScalarField<T>& operator-=(const ScalarField<T>& t_scalar);
    //! Mulitply equals.
    /*! multiply equals operator.  Multiply this scalar field by another.
     *  @param t_scalar A const ScalarField<T>& reference.
     *  @return This field times t_scalar.
     */
    ScalarField<T>& operator*=(const ScalarField<T>& t_scalar);
    //! Divide equals.
    /*! divide equals operator.  Divide this scalar field by another.
     *  @param t_scalar A const ScalarField<T>& reference.
     *  @return This field divided by t_scalar.
     */
    ScalarField<T>& operator/=(const ScalarField<T>& t_scalar);
    //! Access operator
    /*! access operator.  Get the value of the scalar field at a point.
     *  @param t_i The index of the desired point.
     *  @return The value of the scalar field at t_i.
     */
    T& operator()(const size_t& t_i);
    //! Access operator
    /*! access operator.  Get the value of the scalar field at a point.
     *  @param t_i The index of the desired point.
     *  @return The value of the scalar field at t_i.
     */
    const T& operator()(const size_t& t_i) const;
    //! Gradient
    /*! Computes the gradient of the entire field.
     *  @return An std::vector<std::vector<T>> of derivative values
     *  along each direction of the underlying microstates.
     */
    std::vector<std::vector<T>> gradient();
    //! Laplacian
    /*! Computes the Laplacian of the entire field.
     *  @return An std::vector<T> of Laplacian values at each
     *  point of the scalar field.
     */
    std::vector<T> laplacian();
    //! Construct local field values
    /*! Function for generating a vector
     *  of local field values to use for interpolation.
     *  @param t_index The index of the point to construct around.
     *  @return A vector of field values.
     */
    Vector<T> constructLocalFieldValues(size_t t_index);
    //! Construct local field values
    /*! Function for generating vectors and matrices
     *  of local field values to use for interpolation.
     *  @param t_point The point to construct around.
     *  @param t_k The number of neighbors to use.
     *  @return Either vectors or matrices.
     */
    Vector<T> constructLocalFieldValues(const std::vector<T>& t_point,
                                        size_t t_k);

    // T laplacian(size_t index);
    // //--------------------------------------------------------------------------
    // //  Derivatives along the entire grid
    // //--------------------------------------------------------------------------
    // std::vector<std::vector<T>> derivative(size_t n);
    // std::vector<T> derivative(size_t dir, size_t n);
    // std::vector<T> derivative(std::vector<size_t> deriv);
    // //--------------------------------------------------------------------------
    // //  Derivatives at points on the grid
    // //--------------------------------------------------------------------------
    // std::vector<T> derivativePoint(size_t index, size_t n);
    // T derivativePoint(size_t index, size_t dir, size_t n);
    // T derivativePoint(size_t index, std::vector<size_t> deriv);
    // //--------------------------------------------------------------------------
    // //  Derivatives at arbitrary points
    // //--------------------------------------------------------------------------
    // std::vector<T> derivativePoint(std::vector<T> point, size_t n);
    // T derivativePoint(std::vector<T> point, size_t dir, size_t n);
    // T derivativePoint(std::vector<T> point, std::vector<size_t> deriv);
    // //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Methods for calculating integrals
    //--------------------------------------------------------------------------
    void queryNeighborCache(size_t index, size_t k);
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Various functions
    //--------------------------------------------------------------------------
    const std::string summary();
    //--------------------------------------------------------------------------

  private:
    /*! field.  std::vector<T> of field values associated to ugrid. */
    std::vector<T> m_field;
    //--------------------------------------------------------------------------
    //  Objects for differential equation solver
    //--------------------------------------------------------------------------
    size_t _point_cache;
    std::vector<size_t> _neighbor_cache;
    std::vector<double> _distance_cache;

    int _flag;
    std::string _info;
    //--------------------------------------------------------------------------

  };
  //----------------------------------------------------------------------------

  //  Explicit instantiation of double
  template class ScalarField<double>;

}
