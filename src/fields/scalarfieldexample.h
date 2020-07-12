//------------------------------------------------------------------------------
//  scalarfieldexample.h
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
  //  Example scalar fields
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Wave equation in one dimension
  //  f(x) = Acos(kx - wt)
  //       = Acos(k(x - vt))
  //----------------------------------------------------------------------------
  template<typename T>
  class WaveEQ1D : public ScalarField<T>
  {
  public:
    WaveEQ1D();
    WaveEQ1D(std::shared_ptr<UGrid<T>> ugrid);
    WaveEQ1D(std::shared_ptr<UGrid<T>> ugrid, T A, T k, T w);
    //  Getters
    T getA();
    T getk();
    T getw();
    T getv();
    //  Setters
    void setA(T A);
    void setk(T k);
    void setw(T w);
  private:
    T _A;
    T _k;
    T _w;
  };
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Klein Gordon field in 1D
  //----------------------------------------------------------------------------
  template<typename T>
  class KleinGordon1D : public ScalarField<T>
  {
  public:
    KleinGordon1D();
    KleinGordon1D(std::shared_ptr<UGrid<T>> ugrid);
    KleinGordon1D(std::shared_ptr<UGrid<T>> ugrid, T mass);
    //  Getters
    T getMass();
    //  Setters
    void setMass(T mass);
  private:
    T _mass;
  };
  //----------------------------------------------------------------------------

  template class WaveEQ1D<double>;
  template class KleinGordon1D<double>;
}
