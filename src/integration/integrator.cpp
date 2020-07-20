//------------------------------------------------------------------------------
//  integrator.cpp
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
#include "integrator.h"

namespace ET
{
  template<typename T>
  Integrator<T>::Integrator()
  {
  }
  template<typename T>
  Integrator<T>::~Integrator()
  {
  }
  template<typename T>
  Integrator<T>::Integrator(int type)
  {
  }
  template<typename T>
  Integrator<T>::Integrator(std::string type)
  {
  }

  template<typename T>
  Integrator<T>::Integrator(std::shared_ptr<Log> log)
  {
    _log = log;
  }
  template<typename T>
  Integrator<T>::Integrator(std::shared_ptr<Grid<T>>, std::shared_ptr<Log> log)
  {
    _log = log;
  }
  // template<typename T>
  // void Integrator<T>::setF(Vector<T>(ScalarField<T>::*f)(const Vector<T>&,double,Vector<T>))
  // {
  //   _f = f;
  // }

  // //----------------------------------------------------------------------------
  // //  Integration methods
  // //----------------------------------------------------------------------------
  // template<typename T>
  // void Integrator<T>::scalarIntegrate(const Grid<T>& ugrid,
  //                                     ScalarField<T>& field)
  // {
  //   std::vector<T> result(field.getN());
  //   const double dt = 0.01;
  //   //integrate(field.DiffEQ, field.data(), 0.0, 1.0, dt);
  // }
  // //----------------------------------------------------------------------------
  //
  // //----------------------------------------------------------------------------
  // //  Runge-Kutta 4th-order single time step
  // //----------------------------------------------------------------------------
  // template<typename T>
  // void Integrator<T>::scalarRK4Step(const Grid<T>& ugrid,
  //                    ScalarField<T>& field, double dt)
  // {
  //   // //  grab the initial field values
  //   // Vector<T> init(field.getField());
  //   // Vector<T> k0(field.getN(),0.0);
  //   // //  generate k1-k4
  //   // Vector<T> k1 = field.DiffEQ(init,0.0,k0);
  //   // Vector<T> k2 = field.DiffEQ(init,dt/2,k1);
  //   // Vector<T> k3 = field.DiffEQ(init,dt/2,k2);
  //   // Vector<T> k4 = field.DiffEQ(init,dt,k3);
  //   // Vector<T> final = init + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
  //   // //  set new field
  //   // field.setField(final.getVec());
  // }
  // //----------------------------------------------------------------------------
}
