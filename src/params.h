//------------------------------------------------------------------------------
//  params.h
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

#include <map>

namespace ET
{
  enum class MatrixType
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
  //  Interpolator type
  //----------------------------------------------------------------------------
  enum class InterpolatorType
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
  enum class WeightFunctionType
  {
    GAUSSIAN,
  };
  //----------------------------------------------------------------------------


  //----------------------------------------------------------------------------
  //  Integrator type
  //----------------------------------------------------------------------------
  enum class IntegratorType
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
