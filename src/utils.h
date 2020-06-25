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
#   include <cxxabi.h>
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

#include "monomial.hpp"

namespace ET
{
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

    //  method for turning a double into scientific
    //  notation with a certain number of decimal places
    std::string scientific_not(double x, uint64_t dec);
    //  cartesian product between two vectors
    std::vector<std::vector<double> > cartesianProduct(std::vector<double> a,
                                         std::vector<double> b);
    //  generate an index array for all monomials up to order n of some
    //  d dimensional space.
    std::vector<std::vector<uint64_t> > monomial_n(uint64_t dim, uint64_t n);

    //  method for generating a set of taylor polynomials for
    //  a delta p = (x - p).
    std::vector<double> taylorPolynomial(double p, double x, uint64_t n);
    //  generate a monomial expansion about a point p up to order n for some d
    //  dimensiona space according to a point x.
    std::vector<double> taylorMonomialExpansion(std::vector<double> p,
        std::vector<double> x, uint64_t n);
    std::vector<std::string> taylorMonomialFactors(uint64_t dim, uint64_t n);

}
