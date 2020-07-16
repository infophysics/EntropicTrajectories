//------------------------------------------------------------------------------
//  Interpolator.h
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

//------------------------------------------------------------------------------
//  Forward declaration of ScalarField
//------------------------------------------------------------------------------
namespace ET
{
  template<typename T> class Field;
  template<typename T> class ScalarField;
  template<typename T> class VectorField;
}
#include "field.h"
#include "scalarfield.h"
#include "vectorfield.h"

namespace ET
{
  //! LSDriver enum
  /*! Enum for determining the type of least squares driver to use.
   *
   */
  enum class LSDriver
  {
    /*! Enum value: ET::LSDriver::xGELS.*/
    xGELS,
    /*! Enum value: ET::LSDriver::xGELSY.
     *  This uses a complete orthogonal transformation.
     */
    xGELSY,
    /*! Enum value: ET::LSDriver::xGELSD.
     *  This uses the SVD algorithm with divide and conquer.
     */
    xGELSD,
    /*! Enum value: ET::LSDriver::xGELSS.
     *  This uses standard SVD to solve the system.
     */
    xGELSS,
  };

  //! \class Interpolator Class
  /*! A Base class for interpolating on unstructured grids.
   *  Each instance has an associated shared pointer to a UGrid which
   *  gets passed to each derived class.
   */
  template<typename T>
  class Interpolator
  {
  public:
    //! Defualt Constructor
    /*! Default constructor for Interpolator.
     */
    Interpolator();
    //! Destructor
    /*! Destructor for Interpolator.
     */
    ~Interpolator();
    //! Constructor
    /*! constructor for Interpolator that takes a UGrid
     */
    Interpolator(std::shared_ptr<UGrid<T>> t_ugrid);
    //! Constructor
    /*! constructor for Interpolator that takes a Logger
     */
    Interpolator(std::shared_ptr<Log> t_log);
    //! Constructor
    /*! constructor for Interpolator that takes a UGrid and a logger
     */
    Interpolator(std::shared_ptr<UGrid<T>> t_ugrid,
                 std::shared_ptr<Log> t_log);

    /*! Get name.  Get the name of the Interpolator.
     *  @return The name of the Interpolator.
     */
    std::string getName() const;
    /*! Get UGrid.  Get the shared instance of the UGrid.
     *  @return A shared pointer to the UGrid.
     */
    std::shared_ptr<UGrid<T>> getUGrid() const;
    /*! Get Field.  Get the shared instance of the Field.
     *  @return A shared pointer to the Field.
     */
    std::shared_ptr<Field<T>> getField() const;
    /*! Get log.  Get the shared instance of the logger.
     *  @return A shared pointer to the logger
     */
    std::shared_ptr<Log> getLogger() const;
    //! Get flag
    /*! get flag.
     *  @return An int that classifies the type of information stored in
     *  info.
     */
    int getFlag() const;
    //! Get info
    /*! get info.
     *  @return The std::string info that contains relevant information.
     */
    std::string getInfo() const;
    //! Get LSDriver
    /*! get least squares driver.
     *  @return The LSDriver enum that is currently set to be used.
     */
    LSDriver getLSDriver() const;
    //! Set name
    /*! set name.  Sets the name of the Vector.
        @param t_name a std::string for the name of the Vector.
    */
    void setName(std::string t_name);
    //! Set UGrid
    /*! set UGrid.  Sets the shared pointer for the associated UGrid.
     *  @param t_ugrid A shared pointer for a UGrid instance.
     */
    void setUGrid(std::shared_ptr<UGrid<T>> t_ugrid);
    //! Set Field
    /*! set Field.  Sets the shared pointer for the associated Field.
     *  @param t_field A shared pointer for a Field instance.
     */
    void setField(std::shared_ptr<Field<T>> t_field);
    //! Set Logger
    /*! set Logger.  Sets the shared pointer for the associated logger.
     *  @param t_log A shared pointer for a Log instance.
     */
    void setLogger(std::shared_ptr<Log> t_log);
    //! Set flag
    /*! set flag.  Sets the flag pertaining to info.
        @param t_flag an int that classifies the type of information stored
        in m_info.
    */
    void setFlag(int t_flag);
    //! Set info
    /*! set info.  Sets useful information pertaining to Vector.
        @param t_info an std::string containing useful messages.
    */
    void setInfo(std::string t_info);
    //! Set LSDriver
    /*! set least squares driver.  Sets the type of the least squares
     *  method to use.
     *  @param t_type A string denoting the type of LS driver to use.
     */
    void setLSDriver(std::string type);

    //  Derivative functions that must be overloaded in derived classes.

    //! Derivative
    /*! derivative.  Derivative for a point in the UGrid given by index,
     *  of degree t_degree.
     *  @return The nth-derivative at the point given
     *  by the index.
     */
    virtual std::vector<T> derivative(const size_t t_index,
                                      const size_t t_degree);
    //! Derivative
    /*! derivative.  Derivative for a point in the UGrid given by index,
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
   virtual std::vector<T> derivative(const size_t std::vector<T>& point,
                                     const size_t t_degree);
   //! Derivative
   /*! derivative.  Derivative for an arbitrary point,
    *  of degree t_degree and in the direction t_direction.
    *  @return The nth-derivative in the lth-direction at the point given
    *  by the index.
    */
   virtual T derivative(const size_t std::vector<T>& point,
                        const size_t t_degree,
                        const size_t t_direction);

    //--------------------------------------------------------------------------
    //  Gradient functions
    //--------------------------------------------------------------------------

    // //--------------------------------------------------------------------------
    // //  Scalar fields
    // //--------------------------------------------------------------------------
    // //--------------------------------------------------------------------------
    // //  scalarGradientPoint - approximate the gradient at a point
    // //  Arguments:  ugrid   - UGrid<T> pointer
    // //              field   - ScalarField<T> pointer
    // //              index   - index of the point
    // //
    // //  Returns:    std::vector<T> of the gradient.
    // //--------------------------------------------------------------------------
    // std::vector<T>
    // scalarGradientPoint(const std::shared_ptr<UGrid<T>> ugrid,
    //                     const std::shared_ptr<ScalarField<T>> field,
    //                     uint64_t index);
    // //--------------------------------------------------------------------------
    // //  scalarGradientLSPoint - approximate the gradient at a point
    // //                           using the vanilla LS method
    // //  Arguments:  ugrid      - UGrid<T> pointer
    // //              field      - ScalarField<T> pointer
    // //              index      - index of the point
    // //
    // //  Returns:    std::vector<T> of the gradient.
    // //--------------------------------------------------------------------------
    // std::vector<T>
    // scalarGradientLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
    //                        const std::shared_ptr<ScalarField<T>> field,
    //                        uint64_t index);
    // //--------------------------------------------------------------------------
    // //  scalarGradient      - approximate the gradient for an entire field
    // //  Arguments:  ugrid   - UGrid<T> pointer
    // //              field   - ScalarField<T> pointer
    // //
    // //  Returns:    std::vector<std::vector<T>> of the gradient.
    // //--------------------------------------------------------------------------
    // std::vector<std::vector<T>>
    // scalarGradient(const std::shared_ptr<UGrid<T>> ugrid,
    //                const std::shared_ptr<ScalarField<T>> field);
    // //--------------------------------------------------------------------------
    // //  scalarGradientLS      - approximate the gradient for an entire field
    // //                           using the vanilla LS method
    // //  Arguments:  ugrid      - UGrid<T> pointer
    // //              field      - ScalarField<T> pointer
    // //
    // //  Returns:    std::vector<std::vector<T>> of the gradient.
    // //--------------------------------------------------------------------------
    // std::vector<std::vector<T>>
    // scalarGradientLS(const std::shared_ptr<UGrid<T>> ugrid,
    //                   const std::shared_ptr<ScalarField<T>> field);
    // //--------------------------------------------------------------------------
    // //--------------------------------------------------------------------------
    // //  Passing field as a const refernce
    // //--------------------------------------------------------------------------
    // //--------------------------------------------------------------------------
    // //  scalarGradientPoint - approximate the gradient at a point
    // //  Arguments:  ugrid   - UGrid<T> pointer
    // //              field   - const ScalarField<T>& reference
    // //              index   - index of the point
    // //
    // //  Returns:    std::vector<T> of the gradient.
    // //--------------------------------------------------------------------------
    // std::vector<T>
    // scalarGradientPoint(const std::shared_ptr<UGrid<T>> ugrid,
    //                     const ScalarField<T>& field,
    //                     uint64_t index);
    // //--------------------------------------------------------------------------
    // //  scalarGradientLSPoint - approximate the gradient at a point
    // //                           using the vanilla LS method
    // //  Arguments:  ugrid      - UGrid<T> pointer
    // //              field      - const ScalarField<T>& reference
    // //              index      - index of the point
    // //
    // //  Returns:    std::vector<T> of the gradient.
    // //--------------------------------------------------------------------------
    // std::vector<T>
    // scalarGradientLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
    //                        const ScalarField<T>& field,
    //                        uint64_t index);
    // //--------------------------------------------------------------------------
    // //  scalarGradient      - approximate the gradient for an entire field
    // //  Arguments:  ugrid   - UGrid<T> pointer
    // //              field   - const ScalarField<T>& reference
    // //              index   - index of the point
    // //
    // //  Returns:    std::vector<std::vector<T>> of the gradient.
    // //--------------------------------------------------------------------------
    // std::vector<std::vector<T>>
    // scalarGradient(const std::shared_ptr<UGrid<T>> ugrid,
    //                const ScalarField<T>& field);
    // //--------------------------------------------------------------------------
    // //  scalarGradientLS      - approximate the gradient for an entire field
    // //                           using the vanilla LS method
    // //  Arguments:  ugrid      - UGrid<T> pointer
    // //              field      - const ScalarField<T>& reference
    // //              index      - index of the point
    // //
    // //  Returns:    std::vector<std::vector<T>> of the gradient.
    // //--------------------------------------------------------------------------
    // std::vector<std::vector<T>>
    // scalarGradientLS(const std::shared_ptr<UGrid<T>> ugrid,
    //                   const ScalarField<T>& field);
    // //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //  nth-derivatives of scalar field
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //  scalarDerivative       - approximate the derivative for an entire field
    //  Arguments:  ugrid      - UGrid<T> pointer
    //              field      - ScalarField<T> pointer
    //              n          - order of derivative
    //
    //  Returns:    std::vector<std::vector<T>> of the derivatives
    //--------------------------------------------------------------------------
    std::vector<std::vector<T>>
    scalarDerivative(const std::shared_ptr<UGrid<T>> ugrid,
                     const std::shared_ptr<ScalarField<T>> field,
                     uint32_t n);
    //--------------------------------------------------------------------------
    //  scalarDerivative       - approximate the derivative for an entire field
    //                           of order n in the direction dir
    //  Arguments:  ugrid      - UGrid<T> pointer
    //              field      - ScalarField<T> pointer
    //              dir        - direction of the derivative
    //              n          - order of the derivative
    //
    //  Returns:    std::vector<T> of the derivative along dir.
    //--------------------------------------------------------------------------
    std::vector<T>
    scalarDerivative(const std::shared_ptr<UGrid<T>> ugrid,
                          const std::shared_ptr<ScalarField<T>> field,
                          uint32_t dir,
                          uint32_t n);
    //--------------------------------------------------------------------------
    //  scalarDerivative       - approximate the derivative for an entire field
    //                           of order n in the direction dir
    //  Arguments:  ugrid      - UGrid<T> pointer
    //              field      - ScalarField<T> pointer
    //              deriv      - vector of ints denoting direction and order
    //
    //  Returns:    std::vector<T> of the gradient.
    //--------------------------------------------------------------------------
    std::vector<T>
    scalarDerivative(const std::shared_ptr<UGrid<T>> ugrid,
                     const std::shared_ptr<ScalarField<T>> field,
                     std::vector<uint32_t> deriv);
    //--------------------------------------------------------------------------
    //  scalarDerivativePoint  - approximate the derivative for a point
    //                           of order n
    //  Arguments:  ugrid      - UGrid<T> pointer
    //              field      - ScalarField<T> pointer
    //              index      - index of the point
    //              n          - order of the derivative
    //
    //  Returns:    std::vector<T> of the gradient.
    //--------------------------------------------------------------------------
    std::vector<T>
    scalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
                          const std::shared_ptr<ScalarField<T>> field,
                          uint64_t index,
                          uint32_t n);
    //--------------------------------------------------------------------------
    //  scalarDerivativePoint  - approximate the derivative for a point
    //                           of order n
    //  Arguments:  ugrid      - UGrid<T> pointer
    //              field      - ScalarField<T> pointer
    //              point      - std::vector<T> of the point
    //              n          - order of the derivative
    //
    //  Returns:    std::vector<T> of the gradient.
    //--------------------------------------------------------------------------
    std::vector<T>
    scalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
                          const std::shared_ptr<ScalarField<T>> field,
                          std::vector<T> point,
                          uint32_t n);
    //--------------------------------------------------------------------------
    //  scalarDerivativePoint  - approximate the derivative for a point
    //                           of order n in direction dir
    //  Arguments:  ugrid      - UGrid<T> pointer
    //              field      - ScalarField<T> pointer
    //              index      - index of the point
    //              dir        - direction of the derivative
    //              n          - order of the derivative
    //
    //  Returns:    T          - the gradient in direction dir and order n.
    //--------------------------------------------------------------------------
    T scalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
                            const std::shared_ptr<ScalarField<T>> field,
                            uint64_t index,
                            uint32_t dir,
                            uint32_t n);
    //--------------------------------------------------------------------------
    //  scalarDerivativePoint  - approximate the derivative for a point
    //                           of order n in direction dir
    //  Arguments:  ugrid      - UGrid<T> pointer
    //              field      - ScalarField<T> pointer
    //              point      - std::vector<T> of the point
    //              dir        - direction of the derivative
    //              n          - order of the derivative
    //
    //  Returns:    T          - the gradient in direction dir and order n.
    //--------------------------------------------------------------------------
    T scalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
                            const std::shared_ptr<ScalarField<T>> field,
                            std::vector<T> point,
                            uint32_t dir,
                            uint32_t n);
    //--------------------------------------------------------------------------
    //  scalarDerivativePoint  - approximate the derivative for a point
    //                           of order n in direction dir
    //  Arguments:  ugrid      - UGrid<T> pointer
    //              field      - ScalarField<T> pointer
    //              index      - index of the point
    //              deriv      - vector denoting the direction and order
    //
    //  Returns:    T          - the gradient in direction dir and order n.
    //--------------------------------------------------------------------------
    T scalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
                            const std::shared_ptr<ScalarField<T>> field,
                            uint64_t index,
                            std::vector<uint32_t> deriv);
    //--------------------------------------------------------------------------
    //  scalarDerivativePoint  - approximate the derivative for a point
    //                           of order n in direction dir
    //  Arguments:  ugrid      - UGrid<T> pointer
    //              field      - ScalarField<T> pointer
    //              point      - std::vector<T> of the point
    //              deriv      - vector denoting the direction and order
    //
    //  Returns:    T          - the gradient in direction dir and order n.
    //--------------------------------------------------------------------------
    T scalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
                            const std::shared_ptr<ScalarField<T>> field,
                            std::vector<T> point,
                            std::vector<uint32_t> deriv);
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //  Passing field as a const reference
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //  scalarDerivative       - approximate the derivative for an entire field
    //  Arguments:  ugrid      - UGrid<T> pointer
    //              field      - const ScalarField<T>& reference
    //              n          - order of derivative
    //
    //  Returns:    std::vector<std::vector<T>> of the derivatives
    //--------------------------------------------------------------------------
    std::vector<std::vector<T>>
    scalarDerivative(const std::shared_ptr<UGrid<T>> ugrid,
                     const ScalarField<T>& field,
                     uint32_t n);
    //--------------------------------------------------------------------------
    //  scalarDerivative       - approximate the derivative for an entire field
    //                           of order n in the direction dir
    //  Arguments:  ugrid      - UGrid<T> pointer
    //              field      - const ScalarField<T>& reference
    //              dir        - direction of the derivative
    //              n          - order of the derivative
    //
    //  Returns:    std::vector<T> of the derivative along dir.
    //--------------------------------------------------------------------------
    std::vector<T> scalarDerivative(const std::shared_ptr<UGrid<T>> ugrid,
                                    const ScalarField<T>& field,
                                    uint32_t dir,
                                    uint32_t n);
    //--------------------------------------------------------------------------
    //  scalarDerivative       - approximate the derivative for an entire field
    //                           of order n in the direction dir
    //  Arguments:  ugrid      - UGrid<T> pointer
    //              field      - const ScalarField<T>& reference
    //              deriv      - vector of ints denoting direction and order
    //
    //  Returns:    std::vector<T> of the gradient.
    //--------------------------------------------------------------------------
    std::vector<T> scalarDerivative(const std::shared_ptr<UGrid<T>> ugrid,
                                    const ScalarField<T>& field,
                                    std::vector<uint32_t> deriv);
    //--------------------------------------------------------------------------
    //  scalarDerivativePoint  - approximate the derivative for a point
    //                           of order n
    //  Arguments:  ugrid      - UGrid<T> pointer
    //              field      - const ScalarField<T>& reference
    //              index      - index of the point
    //              n          - order of the derivative
    //
    //  Returns:    std::vector<T> of the gradient.
    //--------------------------------------------------------------------------
    std::vector<T> scalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
                                         const ScalarField<T>& field,
                                         uint64_t index,
                                         uint32_t n);
    //--------------------------------------------------------------------------
    //  scalarDerivativePoint  - approximate the derivative for a point
    //                           of order n
    //  Arguments:  ugrid      - UGrid<T> pointer
    //              field      - const ScalarField<T>& reference
    //              point      - std::vector<T> of the point
    //              n          - order of the derivative
    //
    //  Returns:    std::vector<T> of the gradient.
    //--------------------------------------------------------------------------
    std::vector<T> scalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
                                        const ScalarField<T>& field,
                                        std::vector<T> point,
                                        uint32_t n);
    //--------------------------------------------------------------------------
    //  scalarDerivativePoint  - approximate the derivative for a point
    //                           of order n in direction dir
    //  Arguments:  ugrid      - UGrid<T> pointer
    //              field      - const ScalarField<T>& pointer
    //              index      - index of the point
    //              dir        - direction of the derivative
    //              n          - order of the derivative
    //
    //  Returns:    T          - the gradient in direction dir and order n.
    //--------------------------------------------------------------------------
    T scalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
                            const ScalarField<T>& field,
                            uint64_t index,
                            uint32_t dir,
                            uint32_t n);
    //--------------------------------------------------------------------------
    //  scalarDerivativePoint  - approximate the derivative for a point
    //                           of order n in direction dir
    //  Arguments:  ugrid      - UGrid<T> pointer
    //              field      - const ScalarField<T>& pointer
    //              point      - std::vector<T> of the point
    //              dir        - direction of the derivative
    //              n          - order of the derivative
    //
    //  Returns:    T          - the gradient in direction dir and order n.
    //--------------------------------------------------------------------------
    T scalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
                            const ScalarField<T>& field,
                            std::vector<T> point,
                            uint32_t dir,
                            uint32_t n);
    //--------------------------------------------------------------------------
    //  scalarDerivativePoint  - approximate the derivative for a point
    //                           of order n in direction dir
    //  Arguments:  ugrid      - UGrid<T> pointer
    //              field      - const ScalarField<T>& reference
    //              index      - index of the point
    //              deriv      - vector denoting the direction and order
    //
    //  Returns:    T          - the gradient in direction dir and order n.
    //--------------------------------------------------------------------------
    T scalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
                            const ScalarField<T>& field,
                            uint64_t index,
                            std::vector<uint32_t> deriv);
    //--------------------------------------------------------------------------
    //  scalarDerivativePoint  - approximate the derivative for a point
    //                           of order n in direction dir
    //  Arguments:  ugrid      - UGrid<T> pointer
    //              field      - const ScalarField<T>& reference
    //              point      - std::vector<T> of the point
    //              deriv      - vector denoting the direction and order
    //
    //  Returns:    T          - the gradient in direction dir and order n.
    //--------------------------------------------------------------------------
    T scalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
                            const ScalarField<T>& field,
                            std::vector<T> point,
                            std::vector<uint32_t> deriv);
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //  Driver routines for derivatives
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //  Vanilla least squares
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //  scalarDerivativeLSPoint - approximate the derivative for a point
    //                             using the vanilla LS method
    //  Arguments:  ugrid        - UGrid<T> pointer
    //              field        - ScalarField<T> pointer
    //              index        - uint64_t
    //              n            - order of the derivative
    //
    //  Returns:    Vector<T> of the derivatives up to order n.
    //--------------------------------------------------------------------------
    Vector<T>
    scalarDerivativeLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
                            const std::shared_ptr<ScalarField<T>> field,
                            uint64_t index,
                            uint32_t n);
    //--------------------------------------------------------------------------
    //  scalarDerivativeLSPoint - approximate the derivative for a point
    //                             using the vanilla LS method
    //  Arguments:  ugrid        - UGrid<T> pointer
    //              field        - ScalarField<T> pointer
    //              point        - std::vector<T> of the point
    //              n            - order of the derivative
    //
    //  Returns:    Vector<T> of the derivatives up to order n.
    //--------------------------------------------------------------------------
    Vector<T>
    scalarDerivativeLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
                            const std::shared_ptr<ScalarField<T>> field,
                            std::vector<T> point,
                            uint32_t n);
    //--------------------------------------------------------------------------
    //  scalarDerivativeLS      - approximate the derivative for an entire
    //                             field using the vanilla LS method
    //  Arguments:  ugrid        - UGrid<T> pointer
    //              field        - ScalarField<T> pointer
    //              n            - order of the derivative
    //
    //  Returns:    std::vector<T> of the gradient.
    //--------------------------------------------------------------------------
    std::vector<Vector<T>>
    scalarDerivativeLS(const std::shared_ptr<UGrid<T>> ugrid,
                       const std::shared_ptr<ScalarField<T>> field,
                       uint32_t n);
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //  Passing field as a const reference
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //  scalarDerivativeLSPoint - approximate the derivative for a point
    //                             using the vanilla LS method
    //  Arguments:  ugrid        - UGrid<T> pointer
    //              field        - ScalarField<T> pointer
    //              index        - uint64_t
    //              n            - order of the derivative
    //
    //  Returns:    std::vector<T> of the gradient.
    //--------------------------------------------------------------------------
    Vector<T>
    scalarDerivativeLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
                            const ScalarField<T>& field,
                            uint64_t index,
                            uint32_t n);
    //--------------------------------------------------------------------------
    //  scalarDerivativeLSPoint - approximate the derivative for a point
    //                             using the vanilla LS method
    //  Arguments:  ugrid        - UGrid<T> pointer
    //              field        - const ScalarField<T>& reference
    //              point        - std::vector<T> of the point
    //              n            - order of the derivative
    //
    //  Returns:    Vector<T> of the derivatives up to order n.
    //--------------------------------------------------------------------------
    Vector<T>
    scalarDerivativeLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
                            const ScalarField<T>& field,
                            std::vector<T> point,
                            uint32_t n);
    //--------------------------------------------------------------------------
    //  scalarDerivativeLS      - approximate the derivative for an entire
    //                             field using the vanilla LS method
    //  Arguments:  ugrid        - UGrid<T> pointer
    //              field        - ScalarField<T> pointer
    //              n            - order of the derivative
    //
    //  Returns:    std::vector<T> of the gradient.
    //--------------------------------------------------------------------------
    std::vector<Vector<T>>
    scalarDerivativeLS(const std::shared_ptr<UGrid<T>> ugrid,
                       const ScalarField<T>& field,
                       uint32_t n);
    //--------------------------------------------------------------------------
    //  Moving least squares
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //  scalarDerivativeMLSPoint - approximate the derivative for a point
    //                             using the MLS method
    //  Arguments:  ugrid        - UGrid<T> pointer
    //              field        - ScalarField<T> pointer
    //              index        - uint64_t
    //              n            - order of the derivative
    //
    //  Returns:    Vector<T> of the derivatives up to order n.
    //--------------------------------------------------------------------------
    Vector<T>
    scalarDerivativeMLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
                             const std::shared_ptr<ScalarField<T>> field,
                             uint64_t index,
                             uint32_t n);
    //--------------------------------------------------------------------------
    //  scalarDerivativeMLSPoint - approximate the derivative for a point
    //                             using the MLS method
    //  Arguments:  ugrid        - UGrid<T> pointer
    //              field        - ScalarField<T> pointer
    //              point        - coordinates of the point of interest
    //              n            - order of the derivative
    //
    //  Returns:    Vector<T> of the derivatives up to order n.
    //--------------------------------------------------------------------------
    Vector<T>
    scalarDerivativeMLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
                            const std::shared_ptr<ScalarField<T>> field,
                            std::vector<T> point,
                            uint32_t n);
    //--------------------------------------------------------------------------
    //  scalarDerivativeMLS      - approximate the derivative for all points
    //                             using the  MLS method
    //  Arguments:  ugrid        - UGrid<T> pointer
    //              field        - ScalarField<T> pointer
    //              index        - uint64_t
    //              n            - order of the derivative
    //
    //  Returns:    Vector<T> of the derivatives up to order n.
    //--------------------------------------------------------------------------
    std::vector<Vector<T>>
    scalarDerivativeMLS(const std::shared_ptr<UGrid<T>> ugrid,
                        const std::shared_ptr<ScalarField<T>> field,
                        uint32_t n);
    //--------------------------------------------------------------------------
    //  Passing field as a const reference
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //  scalarDerivativeMLSPoint - approximate the derivative for a point
    //                             using the MLS method
    //  Arguments:  ugrid        - UGrid<T> pointer
    //              field        - ScalarField<T> pointer
    //              index        - uint64_t
    //              n            - order of the derivative
    //
    //  Returns:    Vector<T> of the derivatives up to order n.
    //--------------------------------------------------------------------------
    Vector<T>
    scalarDerivativeMLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
                             const ScalarField<T>& field,
                             uint64_t index,
                             uint32_t n);
    //--------------------------------------------------------------------------
    //  scalarDerivativeMLSPoint - approximate the derivative for a point
    //                             using the MLS method
    //  Arguments:  ugrid        - UGrid<T> pointer
    //              field        - const ScalarField<T>& referencer
    //              point        - coordinates of the point of interest
    //              n            - order of the derivative
    //
    //  Returns:    Vector<T> of the derivatives up to order n.
    //--------------------------------------------------------------------------
    Vector<T>
    scalarDerivativeMLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
                           const ScalarField<T>& field,
                           std::vector<T> point,
                           uint32_t n);
    //--------------------------------------------------------------------------
    //  scalarDerivativeMLS      - approximate the derivative for all points
    //                             using the MLS method
    //  Arguments:  ugrid        - UGrid<T> pointer
    //              field        - ScalarField<T> pointer
    //              index        - uint64_t
    //              n            - order of the derivative
    //
    //  Returns:    Vector<T> of the derivatives up to order n.
    //--------------------------------------------------------------------------
    std::vector<Vector<T>>
    scalarDerivativeMLS(const std::shared_ptr<UGrid<T>> ugrid,
                        const ScalarField<T>& field,
                        uint32_t n);
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //  Weighted Moving least squares
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //  scalarDerivativeWMLSPoint - approximate the derivative for a point
    //                             using the WMLS method
    //  Arguments:  ugrid        - UGrid<T> pointer
    //              field        - ScalarField<T> pointer
    //              index        - uint64_t
    //              n            - order of the derivative
    //
    //  Returns:    Vector<T> of the derivatives up to order n.
    //--------------------------------------------------------------------------
    Vector<T>
    scalarDerivativeWMLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
                              const std::shared_ptr<ScalarField<T>> field,
                              uint64_t index,
                              uint32_t n);
    //--------------------------------------------------------------------------
    //  scalarDerivativeWMLSPoint - approximate the derivative for a point
    //                             using the WMLS method
    //  Arguments:  ugrid        - UGrid<T> pointer
    //              field        - ScalarField<T> pointer
    //              point        - std::vector<T> of the point
    //              n            - order of the derivative
    //
    //  Returns:    Vector<T> of the derivatives up to order n.
    //--------------------------------------------------------------------------
    Vector<T>
    scalarDerivativeWMLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
                              const std::shared_ptr<ScalarField<T>> field,
                              std::vector<T> point,
                              uint32_t n);
    //--------------------------------------------------------------------------
    //  scalarDerivativeWMLSPoint - approximate the derivative for a point
    //                             using the  WMLS method
    //  Arguments:  ugrid        - UGrid<T> pointer
    //              field        - ScalarField<T> pointer
    //              index        - uint64_tstd::vector<T> neighbors(points.size());
    //              n            - order of the derivative
    //
    //  Returns:    Vector<T> of the derivatives up to order n.
    //--------------------------------------------------------------------------
    std::vector<Vector<T>>
    scalarDerivativeWMLS(const std::shared_ptr<UGrid<T>> ugrid,
                         const std::shared_ptr<ScalarField<T>> field,
                         uint32_t n);
    //--------------------------------------------------------------------------
    //  Passing field as a const reference
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //  scalarDerivativeWMLSPoint - approximate the derivative for a point
    //                             using the WMLS method
    //  Arguments:  ugrid        - UGrid<T> pointer
    //              field        - ScalarField<T> pointer
    //              index        - uint64_t
    //              n            - order of the derivative
    //
    //  Returns:    Vector<T> of the derivatives up to order n.
    //--------------------------------------------------------------------------
    Vector<T>
    scalarDerivativeWMLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
                              const ScalarField<T>& field,
                              uint64_t index,
                              uint32_t n);
    //--------------------------------------------------------------------------
    //  scalarDerivativeWMLSPoint - approximate the derivative for a point
    //                             using the WMLS method
    //  Arguments:  ugrid        - UGrid<T> pointer
    //              field        - const ScalarField<T>& reference
    //              point        - std::vector<T> of the point
    //              n            - order of the derivative
    //
    //  Returns:    Vector<T> of the derivatives up to order n.
    //--------------------------------------------------------------------------
    Vector<T>
    scalarDerivativeWMLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
                              const ScalarField<T>& field,
                              std::vector<T> point,
                              uint32_t n);
    //--------------------------------------------------------------------------
    //  scalarDerivativeWMLSPoint - approximate the derivative for a point
    //                             using the WMLS method
    //  Arguments:  ugrid        - UGrid<T> pointer
    //              field        - ScalarField<T> pointer
    //              index        - uint64_t
    //              n            - order of the derivative
    //
    //  Returns:    Vector<T> of the derivatives up to order n.
    //--------------------------------------------------------------------------
    std::vector<Vector<T>>
    scalarDerivativeWMLS(const std::shared_ptr<UGrid<T>> ugrid,
                         const ScalarField<T>& field,
                         uint32_t n);
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Helper functions
    //--------------------------------------------------------------------------
    std::vector<Vector<T>>
    xScalarDerivative(const std::shared_ptr<UGrid<T>> ugrid,
                      const std::shared_ptr<ScalarField<T>> field,
                      uint32_t n);
    std::vector<Vector<T>>
    xScalarDerivative(const std::shared_ptr<UGrid<T>> ugrid,
                      const ScalarField<T>& field,
                      uint32_t n);
    Vector<T> xScalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
                                     const std::shared_ptr<ScalarField<T>> field,
                                     uint64_t index, uint32_t n);
    Vector<T> xScalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
                                    const ScalarField<T>& field,
                                    uint64_t index, uint32_t n);
    Vector<T> xScalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
                                     const std::shared_ptr<ScalarField<T>> field,
                                     std::vector<T> point, uint32_t n);
    Vector<T> xScalarDerivativePoint(const std::shared_ptr<UGrid<T>> ugrid,
                                    const ScalarField<T>& field,
                                    std::vector<T> point, uint32_t n);
    Vector<T> xGELSx(Matrix<T> B, Vector<T> u);
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Vector fields
    //--------------------------------------------------------------------------
    std::vector<std::vector<T>>
    vectorDerivative(const std::shared_ptr<UGrid<T>> ugrid,
                     const VectorField<T>& field,
                     uint32_t dir, uint32_t n);
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  various functions
    //--------------------------------------------------------------------------
    const std::string summary();
    //--------------------------------------------------------------------------

  protected:
    /*! Name.  The name of the Interpolator. */
    std::string m_name;
    /*! Logger.  A shared instance of a Logger.*/
    std::shared_ptr<Log> m_log;
    /*! UGrid.  A shared instance of a UGrid.*/
    std::shared_ptr<UGrid<T>> m_ugrid;
    /*! Field.  A shared instance of a Field.*/
    std::shared_ptr<Field<T>> m_field;
    /*! LSDriver.  An enum that denotes the type of least squares
     *  method to use.
     */
    enum LSDriver m_lsdriver;
    /*! Flag.  An in that denotes the type of message stored in info.*/
    int m_flag;
    /*! Info.  A container for storing important information.*/
    std::string m_info;
  };

  template class Interpolator<double>;
}
