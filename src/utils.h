//------------------------------------------------------------------------------
//  utils.h
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
                    abi::__cxa_demangle(typeid(TR).name(), nullptr,
                                               nullptr, nullptr),
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
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Method for turning a double into scientific
    //  notation with a certain number of decimal places
    //--------------------------------------------------------------------------
    template<typename T>
    std::string scientific_not(T x, uint32_t dec);
    template std::string scientific_not<double> (double, uint32_t);

    //--------------------------------------------------------------------------
    //  Cartesian product between two vectors
    //--------------------------------------------------------------------------
    std::vector<std::vector<double> > cartesianProduct(std::vector<double> a,
                                                       std::vector<double> b);

    //--------------------------------------------------------------------------
    //  Multiset coefficient
    //--------------------------------------------------------------------------
    double multisetCoeff(uint32_t dim, uint32_t n);

    //--------------------------------------------------------------------------
    //  Generate an index array for all monomials up to order n of some
    //  d dimensional space.
    //--------------------------------------------------------------------------
    std::vector<std::vector<uint32_t> > monomial_n(uint32_t dim, uint32_t n);

    //--------------------------------------------------------------------------
    //  Return the index of the element of order <n_1,...,n_d>
    //  for some monomial expansion.
    //--------------------------------------------------------------------------
    uint32_t getTaylorIndex(uint32_t dim, uint32_t n,
                            std::vector<uint32_t> term);

    //--------------------------------------------------------------------------
    //  Method for generating a set of taylor polynomials for
    //  a delta p = (x - p).
    //--------------------------------------------------------------------------
    std::vector<double> taylorPolynomial(double p, double x, uint32_t n);

    //--------------------------------------------------------------------------
    //  Generate a monomial expansion about a point p up to order n for some d
    //  dimensiona space according to a point x.
    //--------------------------------------------------------------------------
    std::vector<double> taylorMonomialExpansion(std::vector<double> p,
                                                std::vector<double> x,
                                                uint32_t n);
    //--------------------------------------------------------------------------
    std::vector<std::string> taylorMonomialFactors(uint32_t dim, uint32_t n);

    //--------------------------------------------------------------------------
    //  Checks that the number of elements in each rows of an
    //  array are the same.
    //--------------------------------------------------------------------------
    bool checkConsistency(std::vector<std::vector<double>>& array,
                          std::vector<std::pair<uint32_t,uint32_t>>& rows);

    //--------------------------------------------------------------------------
    //  Error messages
    //  Here we have a set of functions for generating various
    //  error messages
    //--------------------------------------------------------------------------
    std::string MATRIX_INCONSISTENT_ARRAY(std::vector<std::pair<uint32_t,uint32_t>>& rows);
    std::string MATRIX_OUT_OF_BOUNDS(bool axis, const uint32_t& bound,
                                     const uint32_t& attempt,
                                     const std::string& name);
    std::string MATRIX_ADD_INCOMPATIBLE_ROWS(const uint32_t& m1,
                                             const uint32_t& m2,
                                             const std::string& name1,
                                             const std::string& name2);
    std::string MATRIX_ADD_INCOMPATIBLE_COLS(const uint32_t& n1,
                                             const uint32_t& n2,
                                             const std::string& name1,
                                             const std::string& name2);
    std::string MATRIX_SUB_INCOMPATIBLE_ROWS(const uint32_t& m1,
                                             const uint32_t& m2);
    std::string MATRIX_SUB_INCOMPATIBLE_COLS(const uint32_t& n1,
                                             const uint32_t& n2);
    std::string MATRIX_MUL_INCOMPATIBLE(const uint32_t& n1,
                                        const uint32_t& m2);
    std::string MATRIX_ZERO_DIV(const uint32_t& m, const uint32_t& n);

    //--------------------------------------------------------------------------


}
