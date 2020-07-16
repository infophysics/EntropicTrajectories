//------------------------------------------------------------------------------
//  dynamicalsystem.h
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
#include <vector>
#include <string>

#include "ugrid.h"
#include "params.h"
#include "utils.h"
#include "matrix.h"
#include "scalarfield.h"
#include "interpolator.h"
#include "params.h"
#include "log.h"

namespace ET
{
  template<typename T>
  class DynamicalSystem
  {
  public:
    DynamicalSystem();
    ~DynamicalSystem();
    DynamicalSystem(std::string name);

  private:
    std::string _name;
    std::shared_ptr<UGrid<T>> _ugrid;
    std::shared_ptr<Interpolator<T>> _approx;
    std::vector<std::shared_ptr<ScalarField<T>>> _scalarFields;
    std::shared_ptr<Log> _log;
  };

  template class DynamicalSystem<double>;

}
