//------------------------------------------------------------------------------
//  utilities.h
//  The Entropic Trajectories Framework
//  -----------------------------------
//  Copyright (C) [2020] by [N. Carrara]
//  [ncarrara@albany.edu]
//
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
#include <type_traits>
#include <typeinfo>
#ifndef _MSC_VER
#include <cxxabi.h>
#endif
#include <memory>
#include <vector>
#include <string>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include <iterator>
#include <math.h>
#include <iostream>
#include <numeric>

#include "monomial.hpp"

namespace ET
{
  //--------------------------------------------------------------------------
  //  Function for getting the 'type' of a variable
  //  This was taken from a comment thread: https://stackoverflow.com/
  //  questions/81870/is-it-possible-to-print-a-variables-type-in-standard-c
  //--------------------------------------------------------------------------
  template <class T>
  std::string
  type_name()
  {
    typedef typename std::remove_reference<T>::type TR;
    std::unique_ptr<char, void(*)(void*)> own
    (
#ifndef _MSC_VER
      abi::__cxa_demangle(typeid(TR).name(), nullptr, nullptr, nullptr),
#else
      nullptr,
#endif
      std::free
    );
    std::string r = own != nullptr ? own.get() : typeid(TR).name();
    if (std::is_const<TR>::value)
        r += " const";
    if (std::is_volatile<TR>::value)
        r += " volatile";
    if (std::is_lvalue_reference<T>::value)
        r += "&";
    else if (std::is_rvalue_reference<T>::value)
        r += "&&";
    return r;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Function for turning a memory address into a string
  //----------------------------------------------------------------------------
  template<class T>
  std::string address_to_string(const T& object)
  {
    std::ostringstream address;
    address << &object;
    return address.str();
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Method for turning a double into scientific
  //  notation with a certain number of decimal places
  //----------------------------------------------------------------------------
  template<typename T>
  std::string scientific_not(T x, size_t dec);

  //----------------------------------------------------------------------------
  //  Cartesian product between two vectors
  //----------------------------------------------------------------------------
  std::vector<std::vector<double>> cartesianProduct(std::vector<double> a,
                                                    std::vector<double> b);

  //----------------------------------------------------------------------------
  //  Method for generating a set of taylor polynomials for
  //  a delta p = (x - p).
  //----------------------------------------------------------------------------
  std::vector<double> taylorPolynomial(double p, double x, size_t n);

  //----------------------------------------------------------------------------
  //  Checks that the number of elements in each rows of an
  //  array are the same.
  //----------------------------------------------------------------------------
  bool checkConsistency(std::vector<std::vector<double>>& array,
                        std::vector<std::pair<size_t,size_t>>& rows);

  //! \class Monomial class.
  /*! Class for generating d-dimensional monomials of degree n
   *  for use in a local Taylor expansion algorithm.
   */
  class Monomial
  {
  public:
    Monomial();
    ~Monomial();
    Monomial(size_t t_dim);
    Monomial(size_t t_dim, size_t t_deg);

    //  Getters
    size_t getDim();
    size_t getDeg();
    std::vector<std::vector<size_t>> getMono();
    std::vector<std::vector<size_t>> getMono(size_t t_deg);
    std::vector<std::vector<std::string>> getMonoFactors();
    size_t getMultisetCoefficient(size_t t_deg);
    size_t getTaylorIndex(std::vector<size_t> t_term);
    //  Setters
    void setDim(size_t t_dim);
    void setDeg(size_t t_deg);


    void generateMonomial();
    void generateMonomial(size_t t_deg);

    //  Various functions
    std::vector<double> taylorMonomialExpansion(const std::vector<double>& x1,
                                                const std::vector<double>& x2);
    std::vector<double> taylorMonomialExpansion(const std::vector<double>& x1,
                                                const std::vector<double>& x2,
                                                size_t deg);
    const std::string summary();

  private:
    size_t m_dim;
    size_t m_deg;
    size_t m_numElements;
    std::vector<std::vector<size_t>> m_mono;
    std::vector<std::vector<std::string>> m_monoFactors;
    std::vector<size_t> m_multisetCoeffs;
    //  conatiner for message status
    int m_flag;
    //  container for messages
    std::string m_info;
  };

}
