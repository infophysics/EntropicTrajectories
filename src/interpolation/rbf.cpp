//------------------------------------------------------------------------------
//  rbf.cpp
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
#include "rbf.h"

namespace ET
{
  template<typename T>
  RadialBasisFunction<T>::RadialBasisFunction()
  {
  }
  template<typename T>
  RadialBasisFunction<T>::~RadialBasisFunction()
  {
  }

  template<typename T>
  std::shared_ptr<RBFParams> RadialBasisFunction<T>::getRBFParams()
  {
    return std::make_shared<RBFParams>(m_params);
  }
  //----------------------------------------------------------------------------
  //	Construct local RBF matrix
  //----------------------------------------------------------------------------
  template<typename T>
  Matrix<T>
  RadialBasisFunction<T>::constructRBFMatrix(const std::shared_ptr<UGrid<T>> ugrid,
    const std::vector<size_t> neighbors, size_t index)
  {
    //	For a simple gaussian weight, find the distances from index
    //	to each point in neighbors.
    Matrix<T> RBF(neighbors.size());
    for (auto i = 0; i < neighbors.size(); i++)
    {
      for (auto j = 0; j < neighbors.size(); j++)
      {
        RBF(i,j) = xRBFx(ugrid->getPoint(neighbors[i]),
                         ugrid->getPoint(neighbors[j]));
      }
    }
    return RBF;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //	Construct local RBF derivative matrix
  //----------------------------------------------------------------------------
  template<typename T>
  Matrix<T>
  RadialBasisFunction<T>::constructRBFdMatrix(const std::shared_ptr<UGrid<T>> ugrid,
    const std::vector<size_t> neighbors, size_t index, size_t dir)
  {
    //	For a simple gaussian weight, find the distances from index
    //	to each point in neighbors.
    Matrix<T> RBFd(neighbors.size());
    for (auto i = 0; i < neighbors.size(); i++)
    {
      for (auto j = 0; j < neighbors.size(); j++)
      {
        RBFd(i,j) = xRBFdx(ugrid->getPoint(neighbors[i]),
                           ugrid->getPoint(neighbors[j]))[dir];
      }
    }
    return RBFd;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //	Construct RBF matrix
  //----------------------------------------------------------------------------
  template<typename T>
  Matrix<T>
  RadialBasisFunction<T>::constructRBFMatrix(const std::shared_ptr<UGrid<T>> ugrid)
  {
    //	For a simple gaussian weight, find the distances between each point
    Matrix<T> RBF(ugrid->getN());
    for (auto i = 0; i < ugrid->getN(); i++)
    {
      for (auto j = 0; j < ugrid->getN(); j++)
      {
        RBF(i,j) = xRBFx(ugrid->getPoint(i),
                         ugrid->getPoint(j));
      }
    }
    return RBF;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //	Construct RBF vector
  //----------------------------------------------------------------------------
  template<typename T>
  Vector<T>
  RadialBasisFunction<T>::constructRBFVector(const std::shared_ptr<UGrid<T>> ugrid,
                                      std::vector<T> point)
  {
    //	For a simple gaussian weight, find the distances between each point
    Vector<T> RBF(ugrid->getN());
    for (auto i = 0; i < ugrid->getN(); i++) {
      RBF(i) = xRBFx(point,ugrid->getPoint(i));
    }
    return RBF;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //	Construct RBF derivative matrix
  //----------------------------------------------------------------------------
  template<typename T>
  Matrix<T>
  RadialBasisFunction<T>::constructRBFdMatrix(const std::shared_ptr<UGrid<T>> ugrid,
                                              size_t dir)
  {
    //	For a simple gaussian weight, find the distances between each point
    Matrix<T> RBFd(ugrid->getN());
    for (auto i = 0; i < ugrid->getN(); i++)
    {
      for (auto j = 0; j < ugrid->getN(); j++)
      {
        RBFd(i,j) = xRBFdx(ugrid->getPoint(i),
                           ugrid->getPoint(j))[dir];
      }
    }
    return RBFd;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //	Construct RBF derivative vector
  //----------------------------------------------------------------------------
  template<typename T>
  Vector<T>
  RadialBasisFunction<T>::constructRBFdVector(const std::shared_ptr<UGrid<T>> ugrid,
                                              std::vector<T> point, size_t dir)
  {
    //	For a simple gaussian weight, find the distances between each point
    Vector<T> RBFd(ugrid->getN());
    for (auto i = 0; i < ugrid->getN(); i++) {
      RBFd(i) = xRBFdx(ugrid->getPoint(i),point)[dir];
    }
    return RBFd;
  }
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  template<typename T>
  T RadialBasisFunction<T>::xRBFx(const std::vector<T>& p1,
                                  const std::vector<T>& p2)
  {
    if (m_type == RBFKernelType::GAUSSIAN)
    {
      return gaussianKernel(p1,p2,m_params.m_rbfshape);
    }
    else
    {
      return gaussianKernel(p1,p2,m_params.m_rbfshape);
    }
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<T> RadialBasisFunction<T>::xRBFdx(const std::vector<T>& p1,
                                                const std::vector<T>& p2)
  {
    if (m_type == RBFKernelType::GAUSSIAN)
    {
      return gaussianKernelFirstDerivative(p1,p2,m_params.m_rbfshape);
    }
    else
    {
      return gaussianKernelFirstDerivative(p1,p2,m_params.m_rbfshape);
    }
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  RBF functions
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //  Gaussian
  //----------------------------------------------------------------------------
  double gaussianKernel(const std::vector<double>& p1,
                        const std::vector<double>& p2,
                        double shape)
  {
    double dist = L2(p1,p2);
    return exp(-pow((shape*dist),2));
  }
  //----------------------------------------------------------------------------
  //  Gaussian first derivative
  //----------------------------------------------------------------------------
  double gaussianKernelFirstDerivative(const std::vector<double>& p1,
                                       const std::vector<double>& p2,
                                       size_t dir, double shape)
  {
    double diff = p1[dir] - p2[dir];
    double dist = L2(p1,p2);
    return -2*shape*diff*gaussianKernel(p1,p2,shape)/dist;
  }
  //----------------------------------------------------------------------------
  //  Gaussian first derivative
  //----------------------------------------------------------------------------
  std::vector<double>
  gaussianKernelFirstDerivative(const std::vector<double>& p1,
                                const std::vector<double>& p2,
                                double shape)
  {
    std::vector<double> result(p1.size());
    double dist = L2(p1,p2);
    for (auto i = 0; i < p1.size(); i++) {
      result[i] = -2*shape*(p1[i] - p2[i])*gaussianKernel(p1,p2,shape)/dist;
    }
    return result;
  }
  //----------------------------------------------------------------------------
  //  Gaussian second derivative
  //----------------------------------------------------------------------------
  double gaussianKernelSecondDerivative(const std::vector<double>& p1,
                                        const std::vector<double>& p2,
                                        double shape)
  {
    return -2*shape;//*gaussianRBF(val,shape)
           //+ 4*pow((shape*val),2)*gaussianRBF(val,shape);
  }
  //----------------------------------------------------------------------------
  //  Gaussian nth-derivative
  //  This formula was taken from the paper
  //  https://pdfs.semanticscholar.org/88a4/6c09acf598170e98d250a86ccf00f69b1544.pdf
  //----------------------------------------------------------------------------
  double gaussianKernelNthDerivative(const std::vector<double>& p1,
                                     const std::vector<double>& p2,
                                     double shape, size_t n)
  {
    return 0;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Multiquadric Kernels
  //----------------------------------------------------------------------------
  double multiquadricKernel(const std::vector<double>& p1,
                            const std::vector<double>& p2,
                            double shape)
  {
    double dist = L2(p1,p1);
    return std::sqrt(1 + std::pow(shape*dist,2));
  }
  //----------------------------------------------------------------------------
  //  Multiquadric first derivative
  //----------------------------------------------------------------------------
  double multiquadricKernelFirstDerivative(const std::vector<double>& p1,
                                           const std::vector<double>& p2,
                                           double shape)
  {
  }
  //----------------------------------------------------------------------------
  //  Multiquadric second derivative
  //----------------------------------------------------------------------------
  double multiquadricKernelSecondDerivative(const std::vector<double>& p1,
                                            const std::vector<double>& p2,
                                            double shape)
  {
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Inverse Multiquadric Kernels
  //----------------------------------------------------------------------------
  double
  inverseMultiquadricKernel(const std::vector<double>& p1,
                            const std::vector<double>& p2,
                            double shape)
  {
  }
  //----------------------------------------------------------------------------
  //  Inverse Multiquadric first derivative
  //----------------------------------------------------------------------------
  double
  inverseMultiquadricKernelFirstDerivative(const std::vector<double>& p1,
                                           const std::vector<double>& p2,
                                           double shape)
  {
  }
  //----------------------------------------------------------------------------
  //  Inverse Multiquadric second derivative
  //----------------------------------------------------------------------------
  double
  inverseMultiquadricKernelSecondDerivative(const std::vector<double>& p1,
                                            const std::vector<double>& p2,
                                            double shape)
  {
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Inverse Quadratic Kernels
  //----------------------------------------------------------------------------
  double
  inverseQuadraticKernel(const std::vector<double>& p1,
                         const std::vector<double>& p2,
                         double shape)
  {
  }
  //----------------------------------------------------------------------------
  //  Inverse Quadratic first derivative
  //----------------------------------------------------------------------------
  double
  inverseQuadraticKernelFirstDerivative(const std::vector<double>& p1,
                                        const std::vector<double>& p2,
                                        double shape)
  {
  }
  //----------------------------------------------------------------------------
  //  Inverse Quadratic second derivative
  //----------------------------------------------------------------------------
  double
  inverseQuadraticKernelSecondDerivative(const std::vector<double>& p1,
                                         const std::vector<double>& p2,
                                         double shape)
   {
   }
  //----------------------------------------------------------------------------

}
