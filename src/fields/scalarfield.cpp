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
    std::cout << "\nScalar Field created at location " << this;
  }
  template<typename T>
  ScalarField<T>::~ScalarField()
  {
    std::cout << "\nScalar Field at location " << this << " destroyed.";
  }
  template<typename T>
  ScalarField<T>::ScalarField(std::shared_ptr<UGrid<T>> ugrid) : _ugrid(ugrid)
  {
    //std::cout << "\nScalar Field created at location " << this;
    _N = _ugrid->getN();
    _approx = std::make_shared<Approximator<T>>();
    _name = " ";
  }
  template<typename T>
  ScalarField<T>::ScalarField(std::string name, std::shared_ptr<UGrid<T>> ugrid)
  : _name(name), _ugrid(ugrid)
  {
    //std::cout << "\nScalar Field created at location " << this;
    _N = _ugrid->getN();
    _approx = std::make_shared<Approximator<T>>();
  }
  template<typename T>
  ScalarField<T>::ScalarField(std::shared_ptr<UGrid<T>> ugrid, std::vector<T> field)
  : _ugrid(ugrid), _field(field)
  {
    //std::cout << "\nScalar Field created at location " << this;
    _N = ugrid->getN();
    _approx = std::make_shared<Approximator<T>>();
    _name = " ";
  }
  template<typename T>
  ScalarField<T>::ScalarField(std::string name, std::shared_ptr<UGrid<T>> ugrid,
                              std::vector<T> field)
  : _name(name), _ugrid(ugrid), _field(field)
  {
    //std::cout << "\nScalar Field created at location " << this;
    _N = ugrid->getN();
    _approx = std::make_shared<Approximator<T>>();
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
  uint32_t ScalarField<T>::getDim() const
  {
    return _dim;
  }
  template<typename T>
  std::shared_ptr<Approximator<T>> ScalarField<T>::getApproximator() const
  {
    return _approx;
  }
  template<typename T>
  int ScalarField<T>::getFlag() const
  {
    return _flag;
  }
  template<typename T>
  std::string ScalarField<T>::getInfo() const
  {
    return _info;
  }
  template<typename T>
  void ScalarField<T>::setUGrid(std::shared_ptr<UGrid<T>> ugrid)
  {
    _ugrid = ugrid;
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

  //----------------------------------------------------------------------------
  //  Methods for calculating derivatives
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<std::vector<T>> ScalarField<T>::gradient()
  {
    return _approx->scalarGradient(_ugrid,std::make_shared<ScalarField<T>>(*this));
  }
  //----------------------------------------------------------------------------
}
