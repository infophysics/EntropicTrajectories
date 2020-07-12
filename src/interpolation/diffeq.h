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

//------------------------------------------------------------------------------
//  Forward declaration of Approximator, Integrator and ScalarField
//------------------------------------------------------------------------------
namespace ET
{
  template<typename T> class Approximator;
  template<typename T> class ScalarField;
  template<typename T> class Integrator;
}
#include "approximator.h"
#include "scalarfield.h"
#include "integrator.h"

namespace ET
{
  //----------------------------------------------------------------------------
  //  Base class for differential equations
  //----------------------------------------------------------------------------
  template<typename T>
  class diffEQ
  {
  public:
    diffEQ();
    ~diffEQ();
  private:
  };
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Base class for first order differential equations
  //----------------------------------------------------------------------------
  template<typename T>
  class firstOrderODE : public diffEQ<T>
  {
  public:
    firstOrderODE();
    ~firstOrderODE();

    //--------------------------------------------------------------------------
    //  Setters and getters
    //--------------------------------------------------------------------------
    void setScalarField(std::shared_ptr<ScalarField<T>> field);
    void setVectorField(std::shared_ptr<VectorField<T>> field);

    std::shared_ptr<ScalarField<T>> getScalarField();
    std::shared_ptr<VectorField<T>> getVectorField();

    //--------------------------------------------------------------------------
    //  evaluation function
    //--------------------------------------------------------------------------
    std::vector<T>  dt(std::vector<T> point, double t, T k);
  private:
    std::shared_ptr<ScalarField<T>> _sfield;
    std::shared_ptr<VectorField<T>> _vfield;
  };
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Base class for second order differential equations
  //----------------------------------------------------------------------------
  template<typename T>
  class secondOrderODE : public diffEQ<T>
  {
  public:
    secondOrderODE();
    ~secondOrderODE();

    //--------------------------------------------------------------------------
    //  Setters and getters
    //--------------------------------------------------------------------------
    void setScalarField(std::shared_ptr<ScalarField<T>> field);
    void setVectorField(std::shared_ptr<VectorField<T>> field);

    std::shared_ptr<ScalarField<T>> getScalarField();
    std::shared_ptr<VectorField<T>> getVectorField();

    //--------------------------------------------------------------------------
    //  evaluation functions
    //--------------------------------------------------------------------------
    std::vector<T>  dt(std::vector<T> point, double t, T k);
    std::vector<T>  d2t(std::vector<T> point, double t, T k);
  private:
    std::shared_ptr<ScalarField<T>> _sfield;
    std::shared_ptr<VectorField<T>> _vfield;
  };
  //----------------------------------------------------------------------------


  //----------------------------------------------------------------------------
  //  Various first order differential equations
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Heat equation
  //----------------------------------------------------------------------------
  template<typename T>
  class heatEquation : public firstOrderODE<T>
  {
  public:
    heatEquation();
    ~heatEquation();
    //--------------------------------------------------------------------------
    //  Setters and getters
    //--------------------------------------------------------------------------
    void setAlpha(T alpha);
    T getAlpha();

    //--------------------------------------------------------------------------
    //  evaluation function
    //--------------------------------------------------------------------------
    std::vector<T>  dt(std::vector<T> point, double t, T k);
  private:
    T _alpha;
  };
  //----------------------------------------------------------------------------

  template class diffEQ<double>;
  template class firstOrderODE<double>;
  template class secondOrderODE<double>;
  template class heatEquation<double>;

}
