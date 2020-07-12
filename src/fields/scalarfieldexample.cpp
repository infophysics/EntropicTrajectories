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
#include "scalarfieldexample.h"


namespace ET
{
  //----------------------------------------------------------------------------
  //  WaveEQ1D Example
  //----------------------------------------------------------------------------
  template<typename T>
  WaveEQ1D<T>::WaveEQ1D() : ScalarField<T>(), _A(1), _k(1), _w(1)
  {
    //getIntegrator()->setF(&diffEQ);
  }
  template<typename T>
  WaveEQ1D<T>::WaveEQ1D(std::shared_ptr<UGrid<T>> ugrid)
  : ScalarField<T>("WaveEQ1D",ugrid), _A(1), _k(1), _w(1)
  {
    //getIntegrator()->setF(&diffEQ);
  }
  template<typename T>
  WaveEQ1D<T>::WaveEQ1D(std::shared_ptr<UGrid<T>> ugrid, T A, T k, T w)
  : ScalarField<T>("WaveEQ1D",ugrid), _A(A), _k(k), _w(w)
  {
    //getIntegrator()->setF(&diffEQ);
  }
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //  Getters and Setters
  //----------------------------------------------------------------------------
  template<typename T>
  T WaveEQ1D<T>::getA()
  {
    return _A;
  }
  template<typename T>
  T WaveEQ1D<T>::getk()
  {
    return _k;
  }
  template<typename T>
  T WaveEQ1D<T>::getw()
  {
    return _w;
  }
  template<typename T>
  T WaveEQ1D<T>::getv()
  {
    return _w/_k;
  }
  template<typename T>
  void WaveEQ1D<T>::setA(T A)
  {
    _A = A;
  }
  template<typename T>
  void WaveEQ1D<T>::setk(T k)
  {
    _k = k;
  }
  template<typename T>
  void WaveEQ1D<T>::setw(T w)
  {
    _w = w;
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

}
