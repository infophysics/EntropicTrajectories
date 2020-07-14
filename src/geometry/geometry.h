//------------------------------------------------------------------------------
//  geometry.h
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
#include <memory>
#include <vector>
#include <string>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include <iterator>
#include <math.h>
#include <iostream>
#include <numeric>

#include "monomial.hpp"

namespace ET
{
  //----------------------------------------------------------------------------
  //  Binomial coefficient
  //----------------------------------------------------------------------------
  double binomialCoeff(double n, double m);
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Bernoulli number
  //----------------------------------------------------------------------------
  double bernoulli(double n, double m);
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  RBF functions
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //  Gaussian
  //----------------------------------------------------------------------------
  double gaussianRBF(const std::vector<double>& p1,
                     const std::vector<double>& p2,
                     double shape=3.05048);
  //----------------------------------------------------------------------------
  //  Gaussian first derivative
  //----------------------------------------------------------------------------
  double gaussianRBFd(const std::vector<double>& p1,
                     const std::vector<double>& p2,
                     double shape=3.05048);
  //----------------------------------------------------------------------------
  //  Gaussian second derivative
  //----------------------------------------------------------------------------
  double gaussianRBFdd(double val, double shape=3.05048);
  //----------------------------------------------------------------------------
  //  Gaussian nth-derivative
  //  This formula was taken from the paper
  //  https://pdfs.semanticscholar.org/88a4/6c09acf598170e98d250a86ccf00f69b1544.pdf
  //----------------------------------------------------------------------------
  double gaussianRBFdn(double val, double shape, uint32_t n);

  //----------------------------------------------------------------------------
  //  Distance functions
  //----------------------------------------------------------------------------
  double L2(const std::vector<double>& p1, const std::vector<double>& p2);

}
