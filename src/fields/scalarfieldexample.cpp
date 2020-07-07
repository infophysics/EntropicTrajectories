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
  //  Gaussian1D Example
  //----------------------------------------------------------------------------
  template<typename T>
  Gaussian1D<T>::Gaussian1D() : ScalarField<T>(), _mass(0)
  {
  }
  template<typename T>
  Gaussian1D<T>::Gaussian1D(std::shared_ptr<UGrid<T>> ugrid)
  : ScalarField<T>("Gaussian1D",ugrid), _mass(0)
  {
  }
  template<typename T>
  Gaussian1D<T>::Gaussian1D(std::shared_ptr<UGrid<T>> ugrid, T mass)
  : ScalarField<T>("Gaussian1D",ugrid), _mass(mass)
  {
  }
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //  Getters and Setters
  //----------------------------------------------------------------------------
  template<typename T>
  T Gaussian1D<T>::getMass()
  {
    return _mass;
  }
  template<typename T>
  void Gaussian1D<T>::setMass(T mass)
  {
    _mass = mass;
  }
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //  Diff EQ function
  //----------------------------------------------------------------------------
  template<typename T>
  Vector<T> Gaussian1D<T>::diffEQ(const Vector<T>& f,
                                       double dt, Vector<T> k)
  {
    Vector<T> y_i = f;
    Vector<T> d(derivative(0,1));
    Vector<T> y_f = -d;
    return y_f;
  }
  //----------------------------------------------------------------------------

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
  Vector<T> KleinGordon1D<T>::diffEQ(const Vector<T>& f,
                                          double dt, Vector<T> k)
  {
    //  The differential equation must be defined by the user
    std::vector<T> dd = derivative(0,2);
    for (uint64_t i = 0; i < getN(); i++)
    {
      dd[i] -= _mass*(*this)(i);
    }
  }
  //----------------------------------------------------------------------------

}
