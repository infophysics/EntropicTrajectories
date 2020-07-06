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

  //----------------------------------------------------------------------------
  //  Integration methods
  //----------------------------------------------------------------------------
  template<typename T>
  void Integrator<T>::integrate(UGrid<T>& ugrid, ScalarField<T>& field)
  {
    std::vector<T> result(field.getN());
    const double dt = 0.01;
    //integrate(field.diffEQ, field.data(), 0.0, 1.0, dt);
  }
  //----------------------------------------------------------------------------
}
