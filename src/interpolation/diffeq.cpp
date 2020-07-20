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
  DiffEQ<T>::DiffEQ()
  {
  }
  template<typename T>
  DiffEQ<T>::~DiffEQ()
  {
  }
  template<typename T>
  DiffEQ<T>::DiffEQ(std::shared_ptr<Log> log)
  {
  }
  template<typename T>
  DiffEQ<T>::DiffEQ(std::shared_ptr<Grid<T>> ugrid, std::shared_ptr<Log> log)
  {
  }
  //----------------------------------------------------------------------------

  // //----------------------------------------------------------------------------
  // //  Base class for first order differential equations
  // //----------------------------------------------------------------------------
  // template<typename T>
  // FirstOrderODE<T>::FirstOrderODE()
  // {
  // }
  // template<typename T>
  // FirstOrderODE<T>::~FirstOrderODE()
  // {
  // }
  // //----------------------------------------------------------------------------
  //
  // //----------------------------------------------------------------------------
  // //  Setters and getters
  // //----------------------------------------------------------------------------
  // template<typename T>
  // void FirstOrderODE<T>::setScalarField(std::shared_ptr<ScalarField<T>> field)
  // {
  //   _sfield = field;
  // }
  // template<typename T>
  // void FirstOrderODE<T>::setVectorField(std::shared_ptr<VectorField<T>> field)
  // {
  //   _vfield = field;
  // }
  // template<typename T>
  // std::shared_ptr<ScalarField<T>> FirstOrderODE<T>::getScalarField()
  // {
  //   return _sfield;
  // }
  // template<typename T>
  // std::shared_ptr<VectorField<T>> FirstOrderODE<T>::getVectorField()
  // {
  //   return _vfield;
  // }
  // //----------------------------------------------------------------------------
  //
  // //----------------------------------------------------------------------------
  // //  Base class for second order differential equations
  // //----------------------------------------------------------------------------
  // template<typename T>
  // SecondOrderODE<T>::SecondOrderODE()
  // {
  // }
  // template<typename T>
  // SecondOrderODE<T>::~SecondOrderODE()
  // {
  // }
  // //----------------------------------------------------------------------------
  //
  // //----------------------------------------------------------------------------
  // //  Setters and getters
  // //----------------------------------------------------------------------------
  // template<typename T>
  // void SecondOrderODE<T>::setScalarField(std::shared_ptr<ScalarField<T>> field)
  // {
  //   _sfield = field;
  // }
  // template<typename T>
  // void SecondOrderODE<T>::setVectorField(std::shared_ptr<VectorField<T>> field)
  // {
  //   _vfield = field;
  // }
  // template<typename T>
  // std::shared_ptr<ScalarField<T>> SecondOrderODE<T>::getScalarField()
  // {
  //   return _sfield;
  // }
  // template<typename T>
  // std::shared_ptr<VectorField<T>> SecondOrderODE<T>::getVectorField()
  // {
  //   return _vfield;
  // }
  // //----------------------------------------------------------------------------
  //
  // //----------------------------------------------------------------------------
  // //  Various first order differential equations
  // //----------------------------------------------------------------------------
  //
  // //----------------------------------------------------------------------------
  // //  Heat equation
  // //----------------------------------------------------------------------------
  // template<typename T>
  // HeatEquation<T>::HeatEquation()
  // {
  // }
  // template<typename T>
  // HeatEquation<T>::~HeatEquation()
  // {
  // }
  // //----------------------------------------------------------------------------
  //
  // //----------------------------------------------------------------------------
  // //  Setters and getters
  // //----------------------------------------------------------------------------
  // template<typename T>
  // void HeatEquation<T>::setAlpha(T alpha)
  // {
  //   _alpha = alpha;
  // }
  // template<typename T>
  // T HeatEquation<T>::getAlpha()
  // {
  //   return _alpha;
  // }
  // //----------------------------------------------------------------------------
  //
  // //----------------------------------------------------------------------------
  // //  Evaluation functions
  // //----------------------------------------------------------------------------
  // template<typename T>
  // T HeatEquation<T>::dt(std::vector<T> point, double t, T k)
  // {
  //   std::vector<T> d2f_dx2 = this->getScalarField()->derivativePoint(point,2);
  //   return _alpha * std::accumulate(d2f_dx2.begin(),d2f_dx2.end(),0);
  // }
  // //----------------------------------------------------------------------------
  // template<typename T>
  // std::vector<T> HeatEquation<T>::dt(std::vector<std::vector<T>> points,
  //                                    double t, std::vector<T> k)
  // {
  //   std::vector<T> result(this->getScalarField()->getN());
  //   for(uint64_t i = 0; i < this->getScalarField()->getN(); i++)
  //   {
  //     std::vector<T> d2f_dx2 = this->getScalarField()->derivativePoint(points[i],2);
  //     T val = d2f_dx2[0];
  //     result[i] = _alpha * val;
  //   }
  //   return result;
  // }
  // //----------------------------------------------------------------------------
}
