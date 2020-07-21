//------------------------------------------------------------------------------
//  dynamicalsystem.cpp
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
#include "dynamicalsystem.h"


namespace ET
{
  template<typename T>
  DynamicalSystem<T>::DynamicalSystem() : _name("default")
  {
    _log = std::make_shared<Log>();
    _log->init("ET:DS:" + _name,".logs/system_default.txt");
    _log->TRACE("Dynamical System created at location " + address_to_string(*this));
  }
  template<typename T>
  DynamicalSystem<T>::~DynamicalSystem()
  {
    _log->TRACE("Dynamical System destroyed at location " + address_to_string(*this));
  }

  template<typename T>
  DynamicalSystem<T>::DynamicalSystem(std::string name) : _name(name)
  {
    _log = std::make_shared<Log>();
    _log->init("ET:DS:" + _name,".logs/system_" + _name + ".txt");
    _log->TRACE("Dynamical System created at location " + address_to_string(*this));
  }


}
