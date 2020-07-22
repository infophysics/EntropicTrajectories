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

#include "grid.h"
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
    /*! constructor for Field that takes a name
     */
    ScalarField(std::string t_name);
    //! Constructor
    /*! constructor for Field that takes a Logger
     */
    ScalarField(std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for Field that takes a name and a Logger
     */
    ScalarField(std::string t_name, std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for Field that takes a Grid
     */
    ScalarField(std::shared_ptr<Grid<T>> t_grid);
    //! Constructor
    /*! constructor for Field that takes a name and a Grid
     */
    ScalarField(std::string t_name, std::shared_ptr<Grid<T>> t_grid);
    //! Constructor
    /*! constructor for Field that takes a Interpolator
     */
    ScalarField(std::shared_ptr<Interpolator<T>> t_interpolator);
    //! Constructor
    /*! constructor for Field that takes a name and a Interpolator
     */
    ScalarField(std::string t_name,
                std::shared_ptr<Interpolator<T>> t_interpolator);
    //! Constructor
    /*! constructor for Field that takes a Grid and a logger
     */
    ScalarField(std::shared_ptr<Grid<T>> t_grid, std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for Field that takes a name, Grid and a logger
     */
    ScalarField(std::string t_name, std::shared_ptr<Grid<T>> t_grid,
                std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for Field that takes a field
     */
    ScalarField(std::vector<T> t_field);
    //! Constructor
    /*! constructor for Field that takes a name and a field
     */
    ScalarField(std::string t_name, std::vector<T> t_field);
    //! Constructor
    /*! constructor for Field that takes a field and a Logger
     */
    ScalarField(std::vector<T> t_field, std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for Field that takes a name, field and a Logger
     */
    ScalarField(std::string t_name, std::vector<T> t_field,
                std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for Field that takes a field and a Grid
     */
    ScalarField(std::vector<T> t_field, std::shared_ptr<Grid<T>> t_grid);
    //! Constructor
    /*! constructor for Field that takes a name, field and a Grid
     */
    ScalarField(std::string t_name, std::vector<T> t_field,
                std::shared_ptr<Grid<T>> t_grid);
    //! Constructor
    /*! constructor for Field that takes a field and a Interpolator
     */
    ScalarField(std::vector<T> t_field,
                std::shared_ptr<Interpolator<T>> t_interpolator);
    //! Constructor
    /*! constructor for Field that takes a name, field and a Interpolator
     */
    ScalarField(std::string t_name, std::vector<T> t_field,
                std::shared_ptr<Interpolator<T>> t_interpolator);
    //! Constructor
    /*! constructor for Field that takes a field, Grid and a logger
     */
    ScalarField(std::vector<T> t_field,
                std::shared_ptr<Grid<T>> t_grid,
                std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for Field that takes a name, field, Grid and a logger
     */
    ScalarField(std::string t_name, std::vector<T> t_field,
                std::shared_ptr<Grid<T>> t_grid,
                std::shared_ptr<Log> t_log);

    //! Get type
    /*! Get type.  Get the type of field.
     *  @return An FieldType enum
     */
    enum FieldType getType() const;
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
    //! Set Interpolator
    /*! set Interpolator.  Sets the shared pointer for the
     *  associated Interpolator.
     *  @param t_interpolator A shared pointer for a Interpolator instance.
     */
    virtual void setInterpolator(std::shared_ptr<Interpolator<T>> t_interpolator);
    //! Set DiffEQ
    /*! set DiffEQ.  Sets the shared pointer for the associated DiffEQ.
     *  @param t_diffeq A shared pointer for a DiffEQ instance.
     */
    virtual void setDiffEQ(std::shared_ptr<DiffEQ<T>> t_diffeq);
    //! Set Integrator
    /*! set Integrator.  Sets the shared pointer for the associated Integrator.
     *  @param t_integrator A shared pointer for a Integrator instance.
     */
    virtual void setIntegrator(std::shared_ptr<Integrator<T>> t_integrator);
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
    virtual Vector<T> constructLocalFieldValues(size_t t_index);
    //! Construct local field values
    /*! Function for generating vectors and matrices
     *  of local field values to use for interpolation.
     *  @param t_point The point to construct around.
     *  @param t_k The number of neighbors to use.
     *  @return Either vectors or matrices.
     */
    virtual Vector<T>
    constructLocalFieldValues(const std::vector<T>& t_point,
                                    size_t t_k);

    //  Derivative functions
    //! Derivative
    /*! derivative.  Derivative for a point in the Grid given by index,
     *  of degree t_degree.
     *  @return The nth-derivative at the point given
     *  by the index.
     */
    virtual Vector<T> derivative(const size_t t_index,
                         const size_t t_degree);
    //! Derivative
    /*! derivative.  Derivative for a point in the Grid given by index,
     *  of degree t_degree and in the direction t_direction.
     *  @return The nth-derivative in the lth-direction at the point given
     *  by the index.
     */
    virtual T derivative(const size_t t_index,
                 const size_t t_degree,
                 const size_t t_direction);
   //! Derivative
   /*! derivative.  Derivative for an arbitrary point,
    *  of degree t_degree.
    *  @return The nth-derivative at the point given
    *  by the index.
    */
   virtual Vector<T> derivative(const std::vector<T>& point,
                        const size_t t_degree);
   //! Derivative
   /*! derivative.  Derivative for an arbitrary point,
    *  of degree t_degree and in the direction t_direction.
    *  @return The nth-derivative in the lth-direction at the point given
    *  by the index.
    */
   virtual T derivative(const std::vector<T>& point,
                const size_t t_degree,
                const size_t t_direction);

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
    //! Assignment operator
    /*! Virtual function for the assignment operator which prevents slicing.
     */
    virtual ScalarField<T>& operator=(const Field<T>& t_field) {
      if (const ScalarField<T>* scalarfield =
          dynamic_cast<const ScalarField<T>*>(&t_field)) {
            assign(*scalarfield);
      }
      else {
        throw;
      }
      return *this;
    }


  private:
    /*! field.  std::vector<T> of field values associated to grid. */
    std::vector<T> m_field;
    /*! Field Type.  A const enum declaring the type of the Field.
     *  Defaulted to default.
     */
    const enum FieldType m_type {FieldType::SCALAR};
    //--------------------------------------------------------------------------
    //  Objects for differential equation solver
    //--------------------------------------------------------------------------
    size_t _point_cache;
    std::vector<size_t> _neighbor_cache;
    std::vector<double> _distance_cache;

  protected:
    /*! Assignment function.  */
    void assign(const ScalarField<T>& t_scalarfield) {
      Field<T>::assign(t_scalarfield);
      m_field = t_scalarfield.getField();
    }

  };
  //----------------------------------------------------------------------------

  //  Explicit instantiation of double
  template class ScalarField<double>;

}
