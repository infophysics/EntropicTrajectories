//------------------------------------------------------------------------------
//  rbf.h
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

namespace ET
{
  //----------------------------------------------------------------------------
  //  Class for radial basis functions and algorithms
  //  such as the RBF-FD algorithm
  //----------------------------------------------------------------------------
  template<typename T>
  class RadialBasisFunction
  {
  public:
    RadialBasisFunction();
    ~RadialBasisFunction();

    std::shared_ptr<RBFParams> getRBFParams();

    //--------------------------------------------------------------------------
    //  Various functions for constructing RBF Matrices
    //--------------------------------------------------------------------------
    Matrix<T> constructRBFMatrix(const std::shared_ptr<UGrid<T>> ugrid,
                                 const std::vector<size_t> neighbors,
                                 size_t index);
    Matrix<T> constructRBFdMatrix(const std::shared_ptr<UGrid<T>> ugrid,
                                 const std::vector<size_t> neighbors,
                                 size_t index, size_t dir);
    Matrix<T> constructRBFMatrix(const std::shared_ptr<UGrid<T>> ugrid);
    Vector<T> constructRBFVector(const std::shared_ptr<UGrid<T>> ugrid,
                                 std::vector<T> point);
    Matrix<T> constructRBFdMatrix(const std::shared_ptr<UGrid<T>> ugrid,
                                  size_t dir);
    Vector<T> constructRBFdVector(const std::shared_ptr<UGrid<T>> ugrid,
                                  std::vector<T> point,
                                  size_t dir);
    //--------------------------------------------------------------------------
    //  Various helper functions
    //--------------------------------------------------------------------------
    T xRBFx(const std::vector<T>& p1, const std::vector<T>& p2);
    std::vector<T> xRBFdx(const std::vector<T>& p1, const std::vector<T>& p2);
    //--------------------------------------------------------------------------
    //  Algorithms
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //  The Radial Basis Function - Finite Difference algorithm
    //
    //--------------------------------------------------------------------------
    void RBF_FD();

  private:
    enum RBFKernelType m_type{RBFKernelType::GAUSSIAN};
    struct RBFParams m_params;
  };
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  RBF Kernels
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //  Gaussian Kernels
  //----------------------------------------------------------------------------
  double gaussianKernel(const std::vector<double>& p1,
                        const std::vector<double>& p2,
                        double shape=3.05048);
  //----------------------------------------------------------------------------
  //  Gaussian first derivative
  //----------------------------------------------------------------------------
  double gaussianKernelFirstDerivative(const std::vector<double>& p1,
                                       const std::vector<double>& p2,
                                       size_t dir, double shape=3.05048);
  //----------------------------------------------------------------------------
  //  Gaussian first derivative
  //----------------------------------------------------------------------------
  std::vector<double>
  gaussianKernelFirstDerivative(const std::vector<double>& p1,
                                const std::vector<double>& p2,
                                double shape=3.05048);
  //----------------------------------------------------------------------------
  //  Gaussian second derivative
  //----------------------------------------------------------------------------
  double gaussianKernelSecondDerivative(const std::vector<double>& p1,
                                        const std::vector<double>& p2,
                                        double shape=3.05048);
  //----------------------------------------------------------------------------
  //  Gaussian nth-derivative
  //  This formula was taken from the paper
  //  https://pdfs.semanticscholar.org/88a4/6c09acf598170e98d250a86ccf00f69b1544.pdf
  //----------------------------------------------------------------------------
  double gaussianKernelNthDerivative(const std::vector<double>& p1,
                                     const std::vector<double>& p2,
                                     double shape=3.05048, size_t n=0);
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Multiquadric Kernels
  //----------------------------------------------------------------------------
  double multiquadricKernel(const std::vector<double>& p1,
                            const std::vector<double>& p2,
                            double shape=3.05048);
  //----------------------------------------------------------------------------
  //  Multiquadric first derivative
  //----------------------------------------------------------------------------
  double multiquadricKernelFirstDerivative(const std::vector<double>& p1,
                                           const std::vector<double>& p2,
                                           double shape=3.05048);
  //----------------------------------------------------------------------------
  //  Multiquadric second derivative
  //----------------------------------------------------------------------------
  double multiquadricKernelSecondDerivative(const std::vector<double>& p1,
                                            const std::vector<double>& p2,
                                            double shape=3.05048);
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Inverse Multiquadric Kernels
  //----------------------------------------------------------------------------
  double
  inverseMultiquadricKernel(const std::vector<double>& p1,
                            const std::vector<double>& p2,
                            double shape=3.05048);
  //----------------------------------------------------------------------------
  //  Inverse Multiquadric first derivative
  //----------------------------------------------------------------------------
  double
  inverseMultiquadricKernelFirstDerivative(const std::vector<double>& p1,
                                           const std::vector<double>& p2,
                                           double shape=3.05048);
  //----------------------------------------------------------------------------
  //  Inverse Multiquadric second derivative
  //----------------------------------------------------------------------------
  double
  inverseMultiquadricKernelSecondDerivative(const std::vector<double>& p1,
                                            const std::vector<double>& p2,
                                            double shape=3.05048);
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Inverse Quadratic Kernels
  //----------------------------------------------------------------------------
  double
  inverseQuadraticKernel(const std::vector<double>& p1,
                         const std::vector<double>& p2,
                         double shape=3.05048);
  //----------------------------------------------------------------------------
  //  Inverse Quadratic first derivative
  //----------------------------------------------------------------------------
  double
  inverseQuadraticKernelFirstDerivative(const std::vector<double>& p1,
                                        const std::vector<double>& p2,
                                        double shape=3.05048);
  //----------------------------------------------------------------------------
  //  Inverse Quadratic second derivative
  //----------------------------------------------------------------------------
  double
  inverseQuadraticKernelSecondDerivative(const std::vector<double>& p1,
                                         const std::vector<double>& p2,
                                         double shape=3.05048);
  //----------------------------------------------------------------------------
  template class RadialBasisFunction<double>;
}
