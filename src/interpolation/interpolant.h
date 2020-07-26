//------------------------------------------------------------------------------
//  interpolant.h
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
#include <memory>
#include <iostream>
#include <ostream>
#include <stdio.h>
#include <lapacke.h>
#include <cblas.h>
#include <stdint.h>

#include "log.h"
#include "utilities.h"

namespace ET
{
  //! \class Interpolant
  /*! Interpolant class.  Class for storing interpolated functions
   *  and fields.  It is assumed in most Interpolants that the function
   *  to be interpolated is C^0 continuous, except when using piece wise
   *  Interpolants.
   */
  template<typename T>
  class Interpolant
  {
  public:
    //! Default Constructor
    /*! Default constructor for an instance of the base class.
     */
    Interpolant();
    //! Default Destructor
    /*! Default destructor for an instance of the base class.
     */
    ~Interpolant();

    // Getters and Setters
    //! Get Name.
    /*! Get name. Returns the name of the Interpolant.
     *  @return An std::string of the name.
     */
    std::string getName() const;
    //! Get Dim.
    /*! Get dimension.  Returns the dimension of the underlying space.
     *  @return A size_t object specifying the dimension.
     */
    size_t getDim() const;
    //! Get Ranges.
    /*! Get ranges.  Returns the ranges of each dimension.
     *  @return An std::vector<std::vector<T>> of the ranges.
     */
    std::vector<std::vector<T>> getRanges() const;
    //! Get Range.
    /*! Get range.  Returns the range for a given dimension.
     *  @param t_i A size_t denoting the dimension.
     *  @return An std::vector<T> of the Ranges for the dimension t_i.
     */
    std::vector<T> getRange(const size_t t_i) const;
    //! Get Range Min
    /*! Get range min.  Returns the left hand side of the range for a dimension.
     *  @param t_i A size_t object denoting the dimension.
     *  @return A T object specifying the left boundary.
     */
    T getRangeMin(const size_t t_i) const;
    //! Get Range Max
    /*! Get range max.  Returns the right hand side of the range for a dimension.
     *  @param t_i A size_t object denoting the dimension.
     *  @return A T object specifying the right boundary.
     */
    T getRangeMax(const size_t t_i) const;
    //! Set Name
    /*! Set name.  Sets the name of the Interpolant.
     *  @param t_name The name to set.
     */
    void setName(const std::string t_name);
    //! Set Dim
    /*! Set dimension.  Sets the dimension of the Interpolant.
     *  @param t_dim The dimension to set.
     */
    void setDim(const size_t t_dim);
    //! Set Ranges.
    /*! Set ranges.  Sets the ranges for all dimensions.
     *  @param t_ranges A std::vector<std::vector<T>> of the ranges to be set.
     */
    void setRanges(const std::vector<std::vector<T>> t_ranges);
    //! Set range
    /*! Set range.  Sets the range for a given dimension.
     *  @param t_i The index of the dimension to be set.
     *  @param t_range An std::vector<T> of the range.
     */
    void setRange(const size_t t_i, const std::vector<T> t_range);
    //! Set Range Min
    /*! Set range min.  Sets the left hand side of the range for a dimension.
     *  @param t_i The index of the dimension to be set.
     *  @param t_min The left boundary.
     */
    void setRangeMin(const size_t t_i, const T t_min);
    //! Set Range Max
    /*! Set range max.  Sets the right hand side of the range for a dimension.
     *  @param t_i The index of the dimension to be set.
     *  @param t_max The right boundary.
     */
    void setRangeMax(const size_t t_i, const T t_max);
    //  Access operators
    /*! Access operator.  This operator acts as the function call for
     *  the interpolant.
     *  @param t_point A std::vector<T> object denoting the point
     *  to interpolate at.
     *  @return The return value of the interpolant.
     */
    virtual T operator()(const std::vector<T>& t_point);
    //! Range operator.
    /*! range operator.  This operator allows one to access the
     *  ranges in each dimension.
     *  @param t_i The index of the range.
     *  @return An std::vector<T> of the range.
     */
    std::vector<T>& operator[](const size_t& t_i);
    /*! range operator.  This operator allows one to access the
     *  ranges in each dimension.
     *  @param t_i The index of the range.
     *  @return An std::vector<T> of the range.
     */
    const std::vector<T>& operator[](const size_t& t_i) const;

  protected:
    /*! Name.  Name of the interpolant.
     */
    std::string m_name {""};
    /*! Dimension.  The dimension of the underlying microstates.
     */
    size_t m_dim {0};
    /*! Ranges.  Container storing the ranges for each dimension.
     */
    std::vector<std::vector<T>> m_ranges {{{0,0}}};
    /*! Logger.  A shared instance of a Logger.
     */
    std::shared_ptr<Log> m_log;
    /*! Logging system name generator.
     */
    virtual std::string NAME() const {
      return "Interpolant:[" + std::to_string(m_id.id) + "]:" + m_name + ":";
    }
    /*! Unique ID for each instance.
     */
    UniqueID m_id;
  };

  template class Interpolant<double>;

}
