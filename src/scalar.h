//------------------------------------------------------------------------------
//  scalar.h
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

#include "ugrid.h"

namespace ET
{
  template<typename T> class Approximator;
}
#include "approximator.h"
namespace ET
{
  template<typename T>
  class ScalarField
  {
  public:
    ScalarField();
    ~ScalarField();
    ScalarField(UGrid<T>* micro, std::vector<T> field);

    //  Getters
    Approximator<T>* getApproximator();

    //  Setters
    void setDerivative(std::string type);
    //  Methods for calculating derivatives
    Matrix<T> constructBMatrix();


  private:
    UGrid<T>* _micro;
    std::vector<T> _field;
    Approximator<T>* _approx;

  };

  template class ScalarField<double>;
}
