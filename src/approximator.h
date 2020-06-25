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
#include "grid.h"
#include "params.h"
#include "utils.h"
#include "matrix.h"


namespace ET
{
  template<typename T> class ScalarField;
}
#include "scalar.h"
namespace ET
{
  template<typename T>
  class Approximator
  {
  public:
    Approximator();
    ~Approximator();

    //  Setters
    void setDerivative(std::string type);
    void set_k(uint64_t k);
    void set_n(uint64_t n);

    //  Gradient functions
    std::vector<T> gradient(Grid<T>* grid, ScalarField<T>* field, uint64_t index);
    std::vector<T> gradientMLS(Grid<T>* grid, ScalarField<T>* field, uint64_t index);

    Matrix<T> constructBMatrix(Grid<T>* grid, std::vector<uint64_t>* neighbors,
      uint64_t index, uint64_t order);
  private:
    enum ApproxType _type;
    struct ApproxParams _params;
  };


  template class Approximator<double>;
}
