//------------------------------------------------------------------------------
//  interpolant.cpp
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
#include "interpolant.h"

namespace ET
{
  //----------------------------------------------------------------------------
  //  Constructors
  //----------------------------------------------------------------------------
  template<typename T>
  Interpolant<T>::Interpolant()
  : m_name("default"), m_dim(0), m_ranges({{0}})
  {
    m_log = std::make_shared<Log>();
    m_log->init("ET:Interpolant:" + m_name, ".logs/interpolant_" + m_name + ".txt");
		m_log->TRACE(NAME() + " created at location "
		            + address_to_string(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  Interpolant<T>::~Interpolant()
  {
    m_log->TRACE(NAME() + " destroyed at location " + address_to_string(*this));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::string Interpolant<T>::getName() const
  {
    return m_name;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  size_t Interpolant<T>::getDim() const
  {
    return m_dim;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<std::vector<T>> Interpolant<T>::getRanges() const
  {
    return m_ranges;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<T> Interpolant<T>::getRange(const size_t t_i) const
  {
    //  check that index exists
    if (t_i >= m_dim) {
      m_log->ERROR(NAME() + "Attempted to get range at index " + std::to_string(t_i)
                  + " for Inteprolant with dimension " + std::to_string(m_dim));
      m_log->INFO(NAME() + "Returning empty vector.");
      return std::vector<T>();
    }
    else {
      return m_ranges[t_i];
    }
  }
  //----------------------------------------------------------------------------
  template<typename T>
  T Interpolant<T>::getRangeMin(const size_t t_i) const
  {
    //  check that index exists
    if (t_i >= m_dim) {
      m_log->ERROR(NAME() + "Attempted to get range min at indeInterpolatorx " + std::to_string(t_i)
                  + " for Inteprolant with dimension " + std::to_string(m_dim));
      m_log->INFO(NAME() + "Returning zero.");
      return 0;
    }
    else {
      return m_ranges[t_i][0];
    }
  }
  //----------------------------------------------------------------------------
  template<typename T>
  T Interpolant<T>::getRangeMax(const size_t t_i) const
  {
    //  check that index exists
    if (t_i >= m_dim) {
      m_log->ERROR(NAME() + "Attempted to get range max at index " + std::to_string(t_i)
                  + " for Inteprolant with dimension " + std::to_string(m_dim));
      m_log->INFO(NAME() + "Returning zero.");
      return 0;
    }
    else {
      return m_ranges[t_i][0];
    }
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Interpolant<T>::setName(const std::string t_name)
  {
    m_name = t_name;
    m_log->INFO(NAME() + "Set 'name' to " + m_name);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Interpolant<T>::setDim(const size_t t_dim)
  {
    m_dim = t_dim;
    m_log->INFO(NAME() + "Set 'dim' to " + std::to_string(m_dim));
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Interpolant<T>::setRanges(const std::vector<std::vector<T>> t_ranges)
  {
    m_ranges = t_ranges;
    m_log->INFO(NAME() + "Set new ranges.");
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Interpolant<T>::setRange(const size_t t_i, const std::vector<T> t_range)
  {
    if (t_i >= m_dim) {
      m_log->ERROR(NAME() + "Attempted to set ranges at index " + std::to_string(t_i)
                  + " for Inteprolant with dimension " + std::to_string(m_dim));
      return;
    }
    if (t_range.size() > 2) {
      m_log->ERROR(NAME() + "Attempted to set ranges at index " + std::to_string(t_i)
                   + " with array with more than two values.");
      return;
    }
    else {
      m_ranges[t_i] = t_range;
    }
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Interpolant<T>::setRangeMin(const size_t t_i, const T t_min)
  {
    if (t_i >= m_dim) {
      m_log->ERROR(NAME() + "Attempted to set range min at index " + std::to_string(t_i)
                  + " for Inteprolant with dimension " + std::to_string(m_dim));
      return;
    }
    else {
      m_ranges[t_i][0] = t_min;
    }
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void Interpolant<T>::setRangeMax(const size_t t_i, const T t_max)
  {
    if (t_i >= m_dim) {
      m_log->ERROR(NAME() + "Attempted to set range max at index " + std::to_string(t_i)
                  + " for Inteprolant with dimension " + std::to_string(m_dim));
      return;
    }
    else {
      m_ranges[t_i][1] = t_max;
    }
  }
  //----------------------------------------------------------------------------
  template<typename T>
  T Interpolant<T>::operator()(const std::vector<T>& t_point)
  {
    T result;
    return result;
  }
  //----------------------------------------------------------------------------
  template<typename T>
  std::vector<T>& Interpolant<T>::operator[](const size_t& t_i)
  {
    //  check that index exists
    if (t_i >= m_dim) {
      m_log->ERROR(NAME() + "Attempted to get range at index " + std::to_string(t_i)
                  + " for Inteprolant with dimension " + std::to_string(m_dim));
      m_log->INFO(NAME() + "Returning empty vector.");
      std::vector<T> result;
      return result;
    }
    else {
      return m_ranges[t_i];
    }
  }
  //----------------------------------------------------------------------------
  template<typename T>
  const std::vector<T>& Interpolant<T>::operator[](const size_t& t_i) const
  {
    //  check that index exists
    if (t_i >= m_dim) {
      m_log->ERROR(NAME() + "Attempted to get range at index " + std::to_string(t_i)
                  + " for Inteprolant with dimension " + std::to_string(m_dim));
      m_log->INFO(NAME() + "Returning empty vector.");
      return std::vector<T>();
    }
    else {
      return m_ranges[t_i];
    }
  }
  //----------------------------------------------------------------------------
}
