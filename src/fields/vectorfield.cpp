//------------------------------------------------------------------------------
//  vectorfield.cpp
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
#include "vectorfield.h"


namespace ET
{
  //----------------------------------------------------------------------------
  //  VectorField constructors
  //----------------------------------------------------------------------------
  template<typename T>
  VectorField<T>::VectorField()
  : Field<T>()
  {
  }
  //----------------------------------------------------------------------------
  template<typename T>
  VectorField<T>::~VectorField()
  {
  }
  //----------------------------------------------------------------------------
  template<typename T>
  VectorField<T>::VectorField(std::shared_ptr<UGrid<T>> ugrid)
  : Field<T>(ugrid)
  {
  }
  //----------------------------------------------------------------------------
  template<typename T>
  VectorField<T>::VectorField(std::string name,
                              std::shared_ptr<UGrid<T>> ugrid)
  : Field<T>(name,ugrid)
  {
  }
  //----------------------------------------------------------------------------
  template<typename T>
  VectorField<T>::VectorField(std::shared_ptr<UGrid<T>> ugrid,
                              std::vector<std::vector<T>> field)
  : Field<T>(ugrid), m_field(field)
  {
  }
  //----------------------------------------------------------------------------
  template<typename T>
  VectorField<T>::VectorField(std::string name,
                              std::shared_ptr<UGrid<T>> ugrid,
                              std::vector<std::vector<T>> field)
  : Field<T>(name,ugrid), m_field(field)
  {
  }
  //----------------------------------------------------------------------------
  template<typename T>
  VectorField<T>::VectorField(std::shared_ptr<Log> log)
  : Field<T>(log)
  {
  }
  //----------------------------------------------------------------------------
  template<typename T>
  VectorField<T>::VectorField(std::shared_ptr<UGrid<T>> ugrid,
                              std::shared_ptr<Log> log)
  : Field<T>(ugrid,log)
  {
  }
  //----------------------------------------------------------------------------
  template<typename T>
  VectorField<T>::VectorField(std::string name,
                              std::shared_ptr<UGrid<T>> ugrid,
                              std::shared_ptr<Log> log)
  : Field<T>(name,ugrid,log)
  {
  }
  //----------------------------------------------------------------------------
  template<typename T>
  VectorField<T>::VectorField(std::shared_ptr<UGrid<T>> ugrid,
                              std::vector<std::vector<T>> field,
                              std::shared_ptr<Log> log)
  : Field<T>(ugrid,log), m_field(field)
  {
  }
  //----------------------------------------------------------------------------
  template<typename T>
  VectorField<T>::VectorField(std::string name,
                              std::shared_ptr<UGrid<T>> ugrid,
                              std::vector<std::vector<T>> field,
                              std::shared_ptr<Log> log)
  : Field<T>(name,ugrid,log), m_field(field)
  {
  }
  //----------------------------------------------------------------------------
  //  Getters and Setters
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<std::vector<T>> VectorField<T>::getField() const
  {
    return m_field;
  }
  template<typename T>
  std::vector<std::vector<T>>* VectorField<T>::accessField()
  {
    return &m_field;
  }
  template<typename T>
  T* VectorField<T>::data()
  {
    return m_field[0].data();
  }
  template<typename T>
  void VectorField<T>::setField(std::vector<std::vector<T>> field)
  {
    m_field = field;
  }
  //----------------------------------------------------------------------------
  //  Operator overloads
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<T>& VectorField<T>::operator()(const uint32_t& i)
  {
    return m_field[i];
  }
  template<typename T>
  const std::vector<T>& VectorField<T>::operator()(const uint32_t& i) const
  {
    return m_field[i];
  }
  template<typename T>
  T& VectorField<T>::operator()(const uint32_t& i, const uint32_t& j)
  {
    return m_field[i][j];
  }
  template<typename T>
  const T& VectorField<T>::operator()(const uint32_t& i, const uint32_t& j) const
  {
    return m_field[i][j];
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Derivative operations
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Derivative of a vector field along a particular direction
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<std::vector<T>>
  VectorField<T>::derivative(uint32_t dir, uint32_t n)
  {
    //return _approx->vectorDerivative(_ugrid, (*this), dir, n);
  }
  //----------------------------------------------------------------------------
}
