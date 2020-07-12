//------------------------------------------------------------------------------
//  diffeq.cpp
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
#include "diffeq.h"

namespace ET
{
  //----------------------------------------------------------------------------
  //  Base class for differential equations
  //----------------------------------------------------------------------------
  template<typename T>
  diffEQ<T>::diffEQ()
  {
  }
  template<typename T>
  diffEQ<T>::~diffEQ()
  {
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Base class for first order differential equations
  //----------------------------------------------------------------------------
  template<typename T>
  firstOrderODE<T>::firstOrderODE()
  {
  }
  template<typename T>
  firstOrderODE<T>::~firstOrderODE()
  {
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Setters and getters
  //----------------------------------------------------------------------------
  template<typename T>
  void firstOrderODE<T>::setScalarField(std::shared_ptr<ScalarField<T>> field)
  {
    _sfield = field;
  }
  template<typename T>
  void firstOrderODE<T>::setVectorField(std::shared_ptr<VectorField<T>> field)
  {
    _vfield = field;
  }
  template<typename T>
  std::shared_ptr<ScalarField<T>> firstOrderODE<T>::getScalarField()
  {
    return _sfield;
  }
  template<typename T>
  std::shared_ptr<VectorField<T>> firstOrderODE<T>::getVectorField()
  {
    return _vfield;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Base class for second order differential equations
  //----------------------------------------------------------------------------
  template<typename T>
  secondOrderODE<T>::secondOrderODE()
  {
  }
  template<typename T>
  secondOrderODE<T>::~secondOrderODE()
  {
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Setters and getters
  //----------------------------------------------------------------------------
  template<typename T>
  void secondOrderODE<T>::setScalarField(std::shared_ptr<ScalarField<T>> field)
  {
    _sfield = field;
  }
  template<typename T>
  void secondOrderODE<T>::setVectorField(std::shared_ptr<VectorField<T>> field)
  {
    _vfield = field;
  }
  template<typename T>
  std::shared_ptr<ScalarField<T>> secondOrderODE<T>::getScalarField()
  {
    return _sfield;
  }
  template<typename T>
  std::shared_ptr<VectorField<T>> secondOrderODE<T>::getVectorField()
  {
    return _vfield;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Various first order differential equations
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Heat equation
  //----------------------------------------------------------------------------
  template<typename T>
  heatEquation<T>::heatEquation()
  {
  }
  template<typename T>
  heatEquation<T>::~heatEquation()
  {
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Setters and getters
  //----------------------------------------------------------------------------
  template<typename T>
  void heatEquation<T>::setAlpha(T alpha)
  {
    _alpha = alpha;
  }
  template<typename T>
  T heatEquation<T>::getAlpha()
  {
    return _alpha;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Evaluation function
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<T> heatEquation<T>::dt(std::vector<T> point, double t, T k)
  {
    std::vector<T> result(this->getScalarField()->getN());
    for(uint64_t i = 0; i < this->getScalarField()->getN(); i++)
    {
      std::vector<T> d2f_dx2 = this->getScalarField()->derivativePoint(point,2);
      result[i] = _alpha * std::accumulate(d2f_dx2.begin(),d2f_dx2.end(),0);
    }
    return result;
  }
  //----------------------------------------------------------------------------
}
