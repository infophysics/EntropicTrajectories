//------------------------------------------------------------------------------
//  approximator.h
//  The Entropic Trajectories Framework
//  -----------------------------------
//  Copyright (C) [2020] by [N. Carrara, F. Costa, P. Pessoa]
//  [ncarrara@albany.edu,felipecosta.physics@gmail.com,
//    pedroh.pessoa100@gmail.com]
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
#include "params.h"
#include "utils.h"
#include "matrix.h"

//------------------------------------------------------------------------------
//  Forward declaration of ScalarField
//------------------------------------------------------------------------------
namespace ET
{
  template<typename T> class ScalarField;
  template<typename T> class VectorField;
}
#include "scalarfield.h"
#include "vectorfield.h"
namespace ET
{
  //----------------------------------------------------------------------------
  //  Approximator class.
  //----------------------------------------------------------------------------
  template<typename T>
  class Approximator
  {
  public:
    //--------------------------------------------------------------------------
    //  Constructors
    //--------------------------------------------------------------------------
    Approximator();
    ~Approximator();
    Approximator(int type);
    Approximator(std::string type);
    //--------------------------------------------------------------------------
    //  Constructors with shared loggers
    //--------------------------------------------------------------------------
    Approximator(std::shared_ptr<Log> log);
    Approximator(int type, std::shared_ptr<Log> log);
    Approximator(std::string type, std::shared_ptr<Log> log);
    //--------------------------------------------------------------------------
    //  Getters
    //--------------------------------------------------------------------------
    int getApproxType() const;
    ApproxParams getApproxParams() const;
    int getLSDriver() const;
    int getFlag() const;
    std::string getInfo() const;
    std::shared_ptr<Log> getLogger();
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Setters
    //--------------------------------------------------------------------------
    void setApproxType(std::string type);
    void setApproxParams(ApproxParams params);
    void setLSDriver(std::string type);
    //  Set parameters
    void set_k(uint64_t k);
    void set_n(uint64_t n);
    void setFlag(int flag);
    void setInfo(std::string info);
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Gradient functions
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Scalar fields
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //  scalarGradientPoint - approximate the gradient at a point
    //  Arguments:  ugrid   - UGrid<T> pointer
    //              field   - ScalarField<T> pointer
    //              index   - index of the point
    //
    //  Returns:    std::vector<T> of the gradient.
    //--------------------------------------------------------------------------
    std::vector<T>
    scalarGradientPoint(const std::shared_ptr<UGrid<T>> ugrid,
                        const std::shared_ptr<ScalarField<T>> field,
                        uint64_t index);
    //--------------------------------------------------------------------------
    //  scalarGradientMLSPoint - approximate the gradient at a point
    //                           using the vanilla MLS method
    //  Arguments:  ugrid      - UGrid<T> pointer
    //              field      - ScalarField<T> pointer
    //              index      - index of the point
    //
    //  Returns:    std::vector<T> of the gradient.
    //--------------------------------------------------------------------------
    std::vector<T>
    scalarGradientMLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
                           const std::shared_ptr<ScalarField<T>> field,
                           uint64_t index);
    //--------------------------------------------------------------------------
    //  scalarGradient      - approximate the gradient for an entire field
    //  Arguments:  ugrid   - UGrid<T> pointer
    //              field   - ScalarField<T> pointer
    //
    //  Returns:    std::vector<std::vector<T>> of the gradient.
    //--------------------------------------------------------------------------
    std::vector<std::vector<T>>
    scalarGradient(const std::shared_ptr<UGrid<T>> ugrid,
                   const std::shared_ptr<ScalarField<T>> field);
    //--------------------------------------------------------------------------
    //  scalarGradientMLS      - approximate the gradient for an entire field
    //                           using the vanilla MLS method
    //  Arguments:  ugrid      - UGrid<T> pointer
    //              field      - ScalarField<T> pointer
    //
    //  Returns:    std::vector<std::vector<T>> of the gradient.
    //--------------------------------------------------------------------------
    std::vector<std::vector<T>>
    scalarGradientMLS(const std::shared_ptr<UGrid<T>> ugrid,
                      const std::shared_ptr<ScalarField<T>> field);
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //  Passing field as a const refernce
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //  scalarGradientPoint - approximate the gradient at a point
    //  Arguments:  ugrid   - UGrid<T> pointer
    //              field   - const ScalarField<T>& reference
    //              index   - index of the point
    //
    //  Returns:    std::vector<T> of the gradient.
    //--------------------------------------------------------------------------
    std::vector<T>
    scalarGradientPoint(const std::shared_ptr<UGrid<T>> ugrid,
                        const ScalarField<T>& field,
                        uint64_t index);
    //--------------------------------------------------------------------------
    //  scalarGradientMLSPoint - approximate the gradient at a point
    //                           using the vanilla MLS method
    //  Arguments:  ugrid      - UGrid<T> pointer
    //              field      - const ScalarField<T>& reference
    //              index      - index of the point
    //
    //  Returns:    std::vector<T> of the gradient.
    //--------------------------------------------------------------------------
    std::vector<T>
    scalarGradientMLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
                           const ScalarField<T>& field,
                           uint64_t index);
    //--------------------------------------------------------------------------
    //  scalarGradient      - approximate the gradient for an entire field
    //  Arguments:  ugrid   - UGrid<T> pointer
    //              field   - const ScalarField<T>& reference
    //              index   - index of the point
    //
    //  Returns:    std::vector<std::vector<T>> of the gradient.
    //--------------------------------------------------------------------------
    std::vector<std::vector<T>>
    scalarGradient(const std::shared_ptr<UGrid<T>> ugrid,
                   const ScalarField<T>& field);
    //--------------------------------------------------------------------------
    //  scalarGradientMLS      - approximate the gradient for an entire field
    //                           using the vanilla MLS method
    //  Arguments:  ugrid      - UGrid<T> pointer
    //              field      - const ScalarField<T>& reference
    //              index      - index of the point
    //
    //  Returns:    std::vector<std::vector<T>> of the gradient.
    //--------------------------------------------------------------------------
    std::vector<std::vector<T>>
    scalarGradientMLS(const std::shared_ptr<UGrid<T>> ugrid,
                      const ScalarField<T>& field);
    //--------------------------------------------------------------------------
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
    //--------------------------------------------------------------------------
    //  Driver routines for derivatives
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //  scalarDerivativeMLSPoint - approximate the derivative for a point
    //                             using the vanilla MLS method
    //  Arguments:  ugrid        - UGrid<T> pointer
    //              field        - ScalarField<T> pointer
    //              index        - uint64_t
    //              n            - order of the derivative
    //
    //  Returns:    std::vector<T> of the gradient.
    //--------------------------------------------------------------------------
    std::vector<T>
    scalarDerivativeMLSPoint(const std::shared_ptr<UGrid<T>> ugrid,
                             const std::shared_ptr<ScalarField<T>> field,
                             uint64_t index,
                             uint32_t n);
    //--------------------------------------------------------------------------
    //  scalarDerivativeMLS      - approximate the derivative for an entire
    //                             field using the vanilla MLS method
    //  Arguments:  ugrid        - UGrid<T> pointer
    //              field        - ScalarField<T> pointer
    //              n            - order of the derivative
    //
    //  Returns:    std::vector<T> of the gradient.
    //--------------------------------------------------------------------------
    std::vector<std::vector<T>>
    scalarDerivativeMLS(const std::shared_ptr<UGrid<T>> ugrid,
                        const std::shared_ptr<ScalarField<T>> field,
                        uint32_t n);
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //  Helper functions
    //--------------------------------------------------------------------------
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
    //  MLS functions
    //--------------------------------------------------------------------------
    Matrix<T> constructTaylorMatrix(const std::shared_ptr<UGrid<T>> ugrid,
                                    const std::vector<uint64_t> neighbors,
                                    uint64_t index,
                                    uint64_t order);
    Matrix<T> constructTaylorMatrix(const std::shared_ptr<UGrid<T>> ugrid,
                                    const std::vector<uint64_t> neighbors,
                                    uint64_t index,
                                    Monomial& mono);
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  various functions
    //--------------------------------------------------------------------------
    std::string summary();
    //--------------------------------------------------------------------------

  private:
    //--------------------------------------------------------------------------
    //  Basic attributes
    //--------------------------------------------------------------------------
    std::string _name;
    enum ApproxType _type;
    struct ApproxParams _params;
    enum LSDriver _lsdriver;
    //--------------------------------------------------------------------------
    //  Shared objects
    //--------------------------------------------------------------------------
    std::shared_ptr<Log> _log;

    //  conatiner for message status
    int _flag;
    //  container for messages
    std::string _info;
  };


  template class Approximator<double>;
}
