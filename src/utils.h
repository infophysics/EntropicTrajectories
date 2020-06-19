//  various utilities
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
