//------------------------------------------------------------------------------
//  scalarfield.cpp
//  The Entropic Trajectories Framework
//  -----------------------------------
//  Copyright (C) [2020] by [N. Carrara, F. Costa, P. Pessoa]
//  [ncarrara@albany.edu,felipecosta.physics@gmail.com,
//    pedroh.pessoa100@gmail.com]
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
  //  Constructors
  //----------------------------------------------------------------------------
  template<typename T>
  ScalarField<T>::ScalarField()
  {
  }
  template<typename T>
  ScalarField<T>::~ScalarField()
  {
  }
  template<typename T>
  ScalarField<T>::ScalarField(UGrid<T>* micro) : _micro(micro)
  {
    _N = _micro->getN();
    _approx = new Approximator<T>();
    _name = " ";
  }
  template<typename T>
  ScalarField<T>::ScalarField(std::string name, UGrid<T>* micro)
  : _name(name), _micro(micro)
  {
    _N = _micro->getN();
    _approx = new Approximator<T>();
  }
  template<typename T>
  ScalarField<T>::ScalarField(UGrid<T>* micro, std::vector<T> field)
  : _micro(micro), _field(field)
  {
    _N = micro->getN();
    _approx = new Approximator<T>();
    _name = " ";
  }
  template<typename T>
  ScalarField<T>::ScalarField(std::string name, UGrid<T>* micro,
                              std::vector<T> field)
  : _name(name), _micro(micro), _field(field)
  {
    _N = micro->getN();
    _approx = new Approximator<T>();
    _name = " ";
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Getters and Setters
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<T> ScalarField<T>::getField() const
  {
    return _field;
  }
  template<typename T>
  std::vector<T>* ScalarField<T>::accessField()
  {
    return &_field;
  }
  template<typename T>
  T* ScalarField<T>::data()
  {
    return _field.data();
  }
  template<typename T>
  std::string ScalarField<T>::getName() const
  {
    return _name;
  }
  template<typename T>
  uint64_t ScalarField<T>::getN() const
  {
    return _N;
  }
  template<typename T>
  Approximator<T>* ScalarField<T>::getApproximator()
  {
    return _approx;
  }
  template<typename T>
  void ScalarField<T>::setUGrid(UGrid<T>* micro)
  {
    _micro = micro;
  }
  template<typename T>
  void ScalarField<T>::setField(std::vector<T> field)
  {
    _field = field;
  }
  template<typename T>
  void ScalarField<T>::setName(std::string name)
  {
    _name = name;
  }
  template<typename T>
  void ScalarField<T>::setApproxType(std::string type)
  {
    _approx->setApproxType(type);
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Operator overloads
  //----------------------------------------------------------------------------
  template<typename T>
  T& ScalarField<T>::operator()(const uint32_t& i)
  {
    return _field[i];
  }
  template<typename T>
  const T& ScalarField<T>::operator()(const uint32_t& i) const
  {
    return _field[i];
  }
  //----------------------------------------------------------------------------

}
