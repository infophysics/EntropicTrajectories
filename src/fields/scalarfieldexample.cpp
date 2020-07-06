//------------------------------------------------------------------------------
//  scalarfieldexample.cpp
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
#include "scalarfield.h"


namespace ET
{
  //----------------------------------------------------------------------------
  //  KleinGordon1D Example
  //----------------------------------------------------------------------------
  template<typename T>
  KleinGordon1D<T>::KleinGordon1D() : ScalarField<T>(), _mass(0)
  {
  }
  template<typename T>
  KleinGordon1D<T>::KleinGordon1D(std::shared_ptr<UGrid<T>> ugrid)
  : ScalarField<T>("KleinGordon1D",ugrid), _mass(0)
  {
  }
  template<typename T>
  KleinGordon1D<T>::KleinGordon1D(std::shared_ptr<UGrid<T>> ugrid, T mass)
  : ScalarField<T>("KleinGordon1D",ugrid), _mass(mass)
  {
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Getters and Setters
  //----------------------------------------------------------------------------
  template<typename T>
  T KleinGordon1D<T>::getMass()
  {
    return _mass;
  }
  template<typename T>
  void KleinGordon1D<T>::setMass(T mass)
  {
    _mass = mass;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Diff EQ function
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<std::vector<T>> KleinGordon1D<T>::diffEQ()
  {

  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Diff EQ function for a point
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<T> KleinGordon1D<T>::diffEQ(uint64_t index)
  {

  }
  //----------------------------------------------------------------------------
}
