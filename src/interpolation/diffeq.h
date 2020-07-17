//------------------------------------------------------------------------------
//  diffeq.h
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

#include<vector>
#include<string>
#include <iostream>
#include <map>
#include <memory>

#include "ugrid.h"
#include "log.h"

//------------------------------------------------------------------------------
//  Forward declaration of Interpolator, Integrator and ScalarField
//------------------------------------------------------------------------------
// namespace ET
// {
//   template<typename T> class Field;
//   template<typename T> class Interpolator;
//   //template<typename T> class ScalarField;
//   template<typename T> class Integrator;
// }
// #include "field.h"
// #include "interpolator.h"
// //#include "scalarfield.h"
// #include "integrator.h"

namespace ET
{
  //----------------------------------------------------------------------------
  //  Base class for differential equations
  //----------------------------------------------------------------------------
  template<typename T>
  class DiffEQ
  {
  public:
    DiffEQ();
    ~DiffEQ();
    DiffEQ(std::shared_ptr<Log> log);
    DiffEQ(std::shared_ptr<UGrid<T>> ugrid, std::shared_ptr<Log> log);
  private:
  };
  //----------------------------------------------------------------------------

  // //----------------------------------------------------------------------------
  // //  Base class for first order differential equations
  // //----------------------------------------------------------------------------
  // template<typename T>
  // class FirstOrderODE : public DiffEQ<T>
  // {
  // public:
  //   FirstOrderODE();
  //   ~FirstOrderODE();
  //
  //   //--------------------------------------------------------------------------
  //   //  Setters and getters
  //   //--------------------------------------------------------------------------
  //   // void setScalarField(std::shared_ptr<ScalarField<T>> field);
  //   // void setVectorField(std::shared_ptr<VectorField<T>> field);
  //   //
  //   // std::shared_ptr<ScalarField<T>> getScalarField();
  //   // std::shared_ptr<VectorField<T>> getVectorField();
  //
  //   //--------------------------------------------------------------------------
  //   //  evaluation functions
  //   //--------------------------------------------------------------------------
  //   T  dt(std::vector<T> point, double t, T k);
  //   std::vector<T>  dt(std::vector<std::vector<T>> points, double t,
  //                      std::vector<T> k);
  //   //--------------------------------------------------------------------------
  // private:
  //   // std::shared_ptr<ScalarField<T>> _sfield;
  //   // std::shared_ptr<VectorField<T>> _vfield;
  // };
  // //----------------------------------------------------------------------------
  //
  // //----------------------------------------------------------------------------
  // //  Base class for second order differential equations
  // //----------------------------------------------------------------------------
  // template<typename T>
  // class SecondOrderODE : public DiffEQ<T>
  // {
  // public:
  //   SecondOrderODE();
  //   ~SecondOrderODE();
  //
  //   // //--------------------------------------------------------------------------
  //   // //  Setters and getters
  //   // //--------------------------------------------------------------------------
  //   // void setScalarField(std::shared_ptr<ScalarField<T>> field);
  //   // void setVectorField(std::shared_ptr<VectorField<T>> field);
  //   //
  //   // std::shared_ptr<ScalarField<T>> getScalarField();
  //   // std::shared_ptr<VectorField<T>> getVectorField();
  //
  //   //--------------------------------------------------------------------------
  //   //  evaluation functions
  //   //--------------------------------------------------------------------------
  //   T  dt(std::vector<T> point, double t, T k);
  //   std::vector<T>  dt(std::vector<std::vector<T>> points, double t,
  //                      std::vector<T> k);
  //   T  d2t(std::vector<T> point, double t, T k);
  //   std::vector<T>  d2t(std::vector<std::vector<T>> points, double t,
  //                       std::vector<T> k);
  //   //--------------------------------------------------------------------------
  // private:
  //   // std::shared_ptr<ScalarField<T>> _sfield;
  //   // std::shared_ptr<VectorField<T>> _vfield;
  // };
  // //----------------------------------------------------------------------------
  //
  //
  // //----------------------------------------------------------------------------
  // //  Various first order differential equations
  // //----------------------------------------------------------------------------
  //
  // //----------------------------------------------------------------------------
  // //  Heat equation
  // //----------------------------------------------------------------------------
  // template<typename T>
  // class HeatEquation : public FirstOrderODE<T>
  // {
  // public:
  //   HeatEquation();
  //   ~HeatEquation();
  //   //--------------------------------------------------------------------------
  //   //  Setters and getters
  //   //--------------------------------------------------------------------------
  //   void setAlpha(T alpha);
  //   T getAlpha();
  //
  //   //--------------------------------------------------------------------------
  //   //  evaluation functions
  //   //--------------------------------------------------------------------------
  //   T  dt(std::vector<T> point, double t, T k);
  //   std::vector<T>  dt(std::vector<std::vector<T>> points, double t,
  //                      std::vector<T> k);
  //   //--------------------------------------------------------------------------
  // private:
  //   T _alpha;
  // };
  // //----------------------------------------------------------------------------

  template class DiffEQ<double>;
  // template class FirstOrderODE<double>;
  // template class SecondOrderODE<double>;
  // template class HeatEquation<double>;

}
