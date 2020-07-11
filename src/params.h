//------------------------------------------------------------------------------
//  params.h
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

#include <map>

namespace ET
{
  enum MatrixType
  {
    SQUARE,
    SYMMETRIC,
    PERSYMMETRIC,
    CENTROSYMMETRIC,
    ANTI_SYMMETRIC,
    UPPER_TRIANGULAR,
    LOWER_TRIANGULAR,
    DIAGONAL,
    BIDIAGONAL,
    TRIDIAGONAL,
    ANTI_DIAGONAL,
    BAND,
  };

  //----------------------------------------------------------------------------
  //  Least squares driver routine
  //----------------------------------------------------------------------------
  enum LSDriver
  {
    xGELS,  //  default driver
    xGELSY, //  complete orthogonal factorization
    xGELSD, //  SVD with divide and conquer
    xGELSS, //  SVD
  };
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Approximator type
  //----------------------------------------------------------------------------
  enum ApproxType
  {
    LS,
    MLS,
    WMLS,
    RBF,
  };
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Weight function type
  //----------------------------------------------------------------------------
  enum WeightFunctionType
  {
    GAUSSIAN,
  };
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  RBF function type
  //----------------------------------------------------------------------------
  enum RBFType
  {
    GAUSS,
    MULTIQUADRIC,
    INVERSE_QUADRATIC,
    INVERSE_MULTIQUADRIC
  };
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Struct for approximator parameters
  //----------------------------------------------------------------------------
  struct ApproxParams
  {
    //  TODO:: implement a set of parameters for each type.
    //  Basic parameters
    uint64_t k = 3;             //  number of nearest neighbors
    uint64_t n = 3;             //  order of polynomial expansion
    //  type of weight matrix
    enum WeightFunctionType _weight;
    //  type of RBF kernel
    enum RBFType _rbf;
    //  RBF params
    double _rbfshape = 3.05048;

    ApproxParams() {}
  };
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Integrator type
  //----------------------------------------------------------------------------
  enum IntegratorType
  {
    RG4,
    VORONOI,
  };
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Struct for integrator parameters
  //----------------------------------------------------------------------------
  struct IntegratorParams
  {
    //  TODO: implement a set of parameters for each type
  };
  //----------------------------------------------------------------------------

}
