//------------------------------------------------------------------------------
//  integrator.h
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

#include <vector>
#include <iostream>
#include <map>
#include <memory>

#include "rksuite.h"
//#include <boost/numeric/odeint.hpp>

#include "ugrid.h"
#include "params.h"
#include "utils.h"
#include "matrix.h"

//------------------------------------------------------------------------------
//  Forward declaration of ScalarField
//------------------------------------------------------------------------------
namespace ET
{
  template<typename T> class ScalarField;
  template<typename T> class VectorField;
}
#include "scalarfield.h"
#include "vectorfield.h"
namespace ET
{
  //----------------------------------------------------------------------------
  //  Integrator class.
  //----------------------------------------------------------------------------
  template<typename T>
  class Integrator
  {
  public:
    Integrator();
    ~Integrator();
    Integrator(int type);
    Integrator(std::string type);

    Integrator(std::shared_ptr<Log> log);

    void setF(Vector<T>(ScalarField<T>::*f)(const Vector<T>&,double,Vector<T>));

    //--------------------------------------------------------------------------
    //  Integration methods
    //--------------------------------------------------------------------------
    void scalarIntegrate(const UGrid<T>& ugrid, ScalarField<T>& field);
    //--------------------------------------------------------------------------
    void scalarRK4Step(const UGrid<T>& ugrid, ScalarField<T>& field, double dt);

  private:
    IntegratorType _type;
    struct IntegratorParams _params;
    Vector<T>(ScalarField<T>::*_f)(const Vector<T>&,double,Vector<T>);
    std::shared_ptr<Log> _log;
  };
  //----------------------------------------------------------------------------

  template class Integrator<double>;
}
