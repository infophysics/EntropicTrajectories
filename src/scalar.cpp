//------------------------------------------------------------------------------
//  scalar.cpp
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
#include "scalar.h"


namespace ET
{
  template<typename T>
  ScalarField<T>::ScalarField()
  {

  }
  template<typename T>
  ScalarField<T>::~ScalarField()
  {

  }

  template<typename T>
  Approximator<T>* ScalarField<T>::getApproximator()
  {
    return _approx;
  }

  template<typename T>
  void ScalarField<T>::setDerivative(std::string type)
  {
    _approx->setDerivative(type);
  }

  template<typename T>
  ScalarField<T>::ScalarField(UGrid<T>* micro, std::vector<T> field)
  : _micro(micro), _field(field)
  {
    _approx = new Approximator<T>();
  }
}
