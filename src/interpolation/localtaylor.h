//------------------------------------------------------------------------------
//  localtaylor.h
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
#include "utils.h"
#include "matrix.h"

#include "interpolator.h"

namespace ET
{
  //! Search Scheme enum
  /*! Enum for determining the influence domain.
   *
   */
  enum class SearchScheme
  {
    /*! Enum value: ET::SearchScheme::NEAREST_NEIGHBORS. */
    NEAREST_NEIGHBORS,
    /*! Enum value: ET::SearchScheme::RADIUS. */
    RADIUS
  };

  //! \class LocalTaylorInterpolator Class
  /*! A derived class of Interpolator<T> which uses taylor expansions
   *  for interpolation.
   */
  template<typename T>
  class LocalTaylorInterpolator : public Interpolator<T>
  {
  public:
    //! Defualt Constructor
    /*! Default constructor for LTI.
     */
    LocalTaylorInterpolator();
    //! Destructor
    /*! Destructor for LTI
     */
    ~LocalTaylorInterpolator();
    //! Constructor
    /*! constructor for LTI that takes a UGrid
     */
    LocalTaylorInterpolator(std::shared_ptr<UGrid<T>> t_ugrid);
    //! Constructor
    /*! constructor for LTI that takes a Logger
     */
    LocalTaylorInterpolator(std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for LTI that takes a UGrid and a logger
     */
     LocalTaylorInterpolator(std::shared_ptr<UGrid<T>> t_ugrid,
                  std::shared_ptr<Log> t_log);
    /*! Get name.  Get the name of the LTI.
     *  @return The name of the LTI.
     */
    std::string getName();
    /*! Get k.  Get the nearest neighbor value.
     *  @return The nearest neighbor value k.
     */
    size_t getK();
    /*! Get n.  Get the Taylor expansion exponent n.
     *  @return The order of the expansion n.
     */
    size_t getN();
    /*! Get radius.  Get the radius for the influence domain.
     *  @return The radius of the influence domain.
     */
    double getRadius();
    /*! Get search scheme.  Get the search scheme type.
     *  @return The type of search scheme used.
     */
    SearchScheme getSearchScheme();
    /*! Set name.  Set the name of the LTI.
     *  @param t_name The name to assign to the LTI.
     *  @return void
     */
    void setName(std::string t_name);
    /*! Set k.  Set the nearest neighbor value.
     *  @param t_k The nearest neighbor value to set to m_k.
     *  @return void
     */
    void setK(size_t t_k);
    /*! Set n.  Set the Taylor expansion exponent n.
     *  @param t_n The order of the expansion to set to m_n.
     *  @return void
     */
    void setN(size_t t_n);
    /*! Set radius.  Set the radius for the influence domain.
     *  @param t_radius The radius of the influence domain to set to m_radius.
     *  @return void
     */
    void setRadius(double t_radius);
    /*! Set search scheme.  Set the search scheme type.
     *  @param t_searchScheme The search scheme to be set to m_searchScheme.
     *  @return The type of search scheme used.
     */
    void setSearchScheme(std::string t_searchScheme);

    //! LTI Matrix
    /*! Local Taylor Interpolator Matrix constructor that takes in a
     *  point from an associated grid.
     *  @param t_index The index of the point in the associated grid.
     *  @return The LTI Matrix for the point referenced by index.
     */
    Matrix<T> constructLocalTaylorMatrix(const size_t t_index);
    //! LTI Matrix
    /*! Local Taylor Interpolator Matrix constructor that takes in a
     *  arbitrary point.
     *  @param t_point An arbitrary point to build a Taylor Matrix around.
     *  @return The LTI Matrix for the point.
     */
    Matrix<T> constructLocalTaylorMatrix(const std::vector<T> t_point);
    //! LTI Matrix
    /*! Local Taylor Interpolator Matrix constructor that takes in a
     *  point from an associated grid and a value for k.
     *  @param t_index The index of the point in the associated grid.
     *  @param t_k The number of neighbors to use for the influence domain.
     *  @return The LTI Matrix for the point referenced by index.
     */
    Matrix<T> constructLocalTaylorMatrix(const size_t t_index,
                                         const size_t t_k);
    //! LTI Matrix
    /*! Local Taylor Interpolator Matrix constructor that takes in a
     *  arbitrary point and a value for k.
     *  @param t_point An arbitrary point to build a Taylor Matrix around.
     *  @param t_k The number of neighbors to use for the influence domain.
     *  @return The LTI Matrix for the point.
     */
    Matrix<T> constructLocalTaylorMatrix(const std::vector<T> t_point,
                                         const size_t t_k);
    //! LTI Matrix
    /*! Local Taylor Interpolator Matrix constructor that takes in a
    *  point from an associated grid and a value for radius.
    *  @param t_index The index of the point in the associated grid.
    *  @param t_radius The radius to use for the influence domain.
    *  @return The LTI Matrix for the point referenced by index.
    */
    Matrix<T> constructLocalTaylorMatrix(const size_t t_index,
                                         const double t_radius);
    //! LTI Matrix
    /*! Local Taylor Interpolator Matrix constructor that takes in a
    *  arbitrary point and a value for radius.
    *  @param t_point An arbitrary point to build a Taylor Matrix around.
    *  @param t_radius The radius to use for the influence domain.
    *  @return The LTI Matrix for the point.
    */
    Matrix<T> constructLocalTaylorMatrix(const std::vector<T> t_point,
                                         const double t_radius);
    //! LTI Matrix
    /*! Local Taylor Interpolator Matrix constructor that takes in a
    *  point from an associated grid and values for k and n.
    *  @param t_index The index of the point in the associated grid.
    *  @param t_k The number of neighbors to use for the influence domain.
    *  @param t_n The order of the Taylor expansion to be used.
    *  @return The LTI Matrix for the point referenced by index.
    */
    Matrix<T> constructLocalTaylorMatrix(const size_t t_index,
                                         const size_t t_k,
                                         const size_t t_n);
    //! LTI Matrix
    /*! Local Taylor Interpolator Matrix constructor that takes in a
    *  arbitrary point and values for k and n.
    *  @param t_point An arbitrary point to build a Taylor Matrix around.
    *  @param t_k The number of neighbors to use for the influence domain.
    *  @param t_n The order of the Taylor expansion to be used.
    *  @return The LTI Matrix for the point.
    */
    Matrix<T> constructLocalTaylorMatrix(const std::vector<T> t_point,
                                         const size_t t_k,
                                         const size_t t_n);
    //! LTI Matrix
    /*! Local Taylor Interpolator Matrix constructor that takes in a
    *  point from an associated grid and values for radius and n.
    *  @param t_index The index of the point in the associated grid.
    *  @param t_radius The radius to use for the influence domain.
    *  @param t_n The order of the Taylor expansion to be used.
    *  @return The LTI Matrix for the point referenced by index.
    */
    Matrix<T> constructLocalTaylorMatrix(const size_t t_index,
                                         const double t_radius,
                                         const size_t t_n);
    //! LTI Matrix
    /*! Local Taylor Interpolator Matrix constructor that takes in a
    *  arbitrary point and values for radius and n.
    *  @param t_point An arbitrary point to build a Taylor Matrix around.
    *  @param t_radius The radius to use for the influence domain.
    *  @param t_n The order of the Taylor expansion to be used.
    *  @return The LTI Matrix for the point.
    */
    Matrix<T> constructLocalTaylorMatrix(const std::vector<T> t_point,
                                         const double t_radius,
                                         const size_t t_n);

    //  Overloaded derivative functions

    //! Derivative
    /*! derivative.  Derivative for a point in the UGrid given by index,
    *  of degree t_degree.
    *  @return The nth-derivative at the point given
    *  by the index.
    */
    std::vector<T> derivative(const size_t t_index,
                              const size_t t_degree);
    //! Derivative
    /*! derivative.  Derivative for a point in the UGrid given by index,
    *  of degree t_degree and in the direction t_direction.
    *  @return The nth-derivative in the lth-direction at the point given
    *  by the index.
    */
    T derivative(const size_t t_index,
                 const size_t t_degree,
                 const size_t t_direction);
    //! Derivative
    /*! derivative.  Derivative for an arbitrary point,
    *  of degree t_degree.
    *  @return The nth-derivative at the point given
    *  by the index.
    */
    std::vector<T> derivative(const size_t std::vector<T>& point,
                              const size_t t_degree);
    //! Derivative
    /*! derivative.  Derivative for an arbitrary point,
    *  of degree t_degree and in the direction t_direction.
    *  @return The nth-derivative in the lth-direction at the point given
    *  by the index.
    */
    T derivative(const size_t std::vector<T>& point,
                 const size_t t_degree,
                 const size_t t_direction);
  private:
    /*! Name.  Name of the LTI.  Defaulted to empty string.
     */
    std::string m_name {""};
    /*! k.  Number of nearest neighbors to use in interpolation.
        Defaulted to k = 3
     */
    size_t m_k {3};
    /*! n.  Number of terms to use in the Taylor expansion.
        Defaulted to n = 3.
     */
    size_t m_n {3};
    /*! radius.  Search radius to use for finding neighbors.
        Defaulted to 0.
     */
    double m_radius {0.0};
    /*! search scheme.  Scheme to use when constructing influence domain.
        Defaulted to using k nearest neighbors.
     */
    enum SearchScheme m_searchScheme {SearchScheme::NEAREST_NEIGHBORS};
    /*! Monomial.  An instance of the monomial class for constructing
     *  Taylor monomials.  Defaulted to use k=3 and n=3.
     */
    Monomial m_monomial {3,3};
  };

}
