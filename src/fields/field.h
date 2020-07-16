//------------------------------------------------------------------------------
//  field.h
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
#include <memory>

#include "ugrid.h"
#include "log.h"

//------------------------------------------------------------------------------
//  Forward declaration of Interpolator, Integrator and diffEQ
//------------------------------------------------------------------------------
namespace ET
{
  template<typename T> class Interpolator;
  template<typename T> class diffEQ;
  template<typename T> class Integrator;
}
#include "interpolator.h"
#include "diffeq.h"
#include "integrator.h"

namespace ET
{
  //! Field Class
  /*! A Base class for various types of fields, such as ET::ScalarField,
   *  ET::VectorField and ET::FrameField.
   */
  template<typename T>
  class Field
  {
  public:
    Field();
    ~Field();
    Field(std::shared_ptr<Log> t_log);
    Field(std::shared_ptr<UGrid<T>> t_ugrid);
    Field(std::shared_ptr<Interpolator<T>> t_interpolator);
    Field(std::shared_ptr<UGrid<T>> t_ugrid, std::shared_ptr<Log> t_log);

  protected:
    /*! Name.  The name of the Interpolator. */
    std::string m_name {""};
    /*! Dim.  Dimension of the field. */
    size_t m_dim {0};
    /*! N.  Number of samples of the field. */
    size_t m_N {0};
    /*! Logger.  A shared instance of a Logger.*/
    std::shared_ptr<Log> m_log;
    /*! UGrid.  A shared instance of a UGrid.*/
    std::shared_ptr<UGrid<T>> m_ugrid;
    /*! Interpolator.  A shared instance of an Interpolator. */
    std::shared_ptr<Interpolator<T>> m_interpolator;
    /*! diffEQ.  A shared instance of a diffEQ. */
    std::shared_ptr<diffEQ<T>> m_diffeq;
    /*! Integrator.  A shared instance of an Integrator. */
    std::shared_ptr<Integrator<T>> m_integrator;
  };

}
