//------------------------------------------------------------------------------
//  geometry.cpp
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
#include "geometry.h"

namespace ET
{
  //----------------------------------------------------------------------------
  //  Binomial coefficient
  //----------------------------------------------------------------------------
  double binomialCoeff(double n, double m)
  {
    return (double)i4_choose((int)n,(int)m);
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Bernoulli number
  //----------------------------------------------------------------------------
  double bernoulli(double n, double m)
  {
    double result = 0;
    for (uint32_t i = 0; i < n; i++)
    {
      for (uint32_t j = 0; j < i; j++)
      {
        double sign = pow(-1,j);
        double coeff = binomialCoeff(j,i);
        double num = pow((1 + j),n);
        double den = (m + 1);
        result += sign*coeff*num/den;
      }
    }
    return result;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Distance functions
  //----------------------------------------------------------------------------
  double L2(const std::vector<double>& p1, const std::vector<double>& p2)
  {
    double dist = 0;
    for (uint32_t i = 0; i < p1.size(); i++)
    {
      dist += pow((p1[i] - p2[i]),2);
    }
    return sqrt(dist);
  }
  //----------------------------------------------------------------------------

}
