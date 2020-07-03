//------------------------------------------------------------------------------
//  utils.cpp
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
#include "utils.h"

namespace ET
{
  //----------------------------------------------------------------------------
  //  Function for turning a number into scientific notation
  //----------------------------------------------------------------------------
  template<typename T>
  std::string scientific_not(T x, uint32_t dec)
  {
    T y;
    int exp = int(floor(log10(abs(x))));
    T factor = pow(10,exp);
    if (factor != 0)
      y = x / abs(factor);
    else
      y = x;
    std::stringstream s;
    s << std::fixed << std::setprecision(dec) << y;
    std::string result = s.str() + "e";
    if (exp < -150)
      exp = 0;
    if (exp > 150)
    {
      if ( x < 0)
        return "  -inf   ";
      else
        return "   inf   ";
    }
    if (exp >= 0)
      if (exp < 10)
        result += "+0" + std::to_string(abs(exp));
      else
        result += "+" + std::to_string(abs(exp));
      else
    if (exp > -10)
      result += "-0" + std::to_string(abs(exp));
    else
      result += "-" + std::to_string(abs(exp));
    return result;
  }
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Function for generating a Cartesian product
    //--------------------------------------------------------------------------
    std::vector<std::vector<double> > cartesianProduct(std::vector<double> a,
                                                       std::vector<double> b)
    {
        std::vector<std::vector<double> > cart(a.size());
        for (uint32_t i = 0; i < a.size(); i++)
        {
            cart[i].resize(b.size());
            for (uint32_t j = 0; j < b.size(); j++)
            {
                cart[i][j] = a[i]*b[j];
            }
        }
        return cart;
    }
    //--------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //    Monomial class
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //    Constructors
  //----------------------------------------------------------------------------
  Monomial::Monomial()
  {
  }
  Monomial::~Monomial()
  {
    //std::cout << "\nMonomial destroyed.";
  }
  Monomial::Monomial(uint32_t dim) : _dim(dim)
  {
  }
  Monomial::Monomial(uint32_t dim, uint32_t deg) : _dim(dim), _deg(deg)
  {
    generateMonomial();
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Getters and Setters
  //----------------------------------------------------------------------------
  uint32_t Monomial::getDim()
  {
    return _dim;
  }
  uint32_t Monomial::getDeg()
  {
    return _deg;
  }
  std::vector<std::vector<uint32_t>> Monomial::getMono()
  {
    return _mono;
  }
  std::vector<std::vector<uint32_t>> Monomial::getMono(uint32_t deg)
  {
    generateMonomial(_deg);
    return _mono;
  }
  std::vector<std::vector<std::string>> Monomial::getMonoFactors()
  {
    return _monoFactors;
  }
  uint32_t Monomial::getMultisetCoefficient(uint32_t deg)
  {
    //  use the i4_choose function from monomial.hpp
    //  the degree of the kth-monomial for a space of
    //  dimension d is the binomial coefficent (k + d - 1 // k)
    return i4_choose(deg + _dim - 1,deg);
  }
  uint32_t Monomial::getTaylorIndex(std::vector<uint32_t> term)
  {
    if (term.size() != _dim)
    {
      std::cout << "Number of elements in exp.";
      std::cout << " array does not equal dimension!" << std::endl;
      return 0;
    }
    std::vector<int> x(_dim);
    std::vector<int> y(_dim);
    for (uint32_t i = 0; i < _dim; i++)
    {
      x[i] = 0;
      y[i] = term[i];
    }
    int index = 0;
    while (x != y)
    {
      mono_between_next_grlex(_dim,0,_deg,x.data());
      index++;
    }
    return index;
  }
  void Monomial::setDim(uint32_t dim)
  {
    _dim = dim;
  }
  void Monomial::setDeg(uint32_t deg)
  {
    _deg = deg;
  }
  void Monomial::generateMonomial()
  {
    int x[_dim];
    for (uint32_t i = 0; i < _dim; i++)
    {
      x[i] = 0;
    }
    std::vector<std::vector<uint32_t>> mono;
    std::vector<uint32_t> temp(x, x + _dim);
    mono.push_back(temp);
    int i = 1;
    for ( ; ; )
    {
      if ( x[0] == _deg )
      {
        break;
      }
      mono_between_next_grlex(_dim, 0, _deg, x);
      std::vector<uint32_t> temp(x, x + _dim);
      mono.push_back(temp);
      i = i + 1;
    }
    _mono = mono;
  }
  void Monomial::generateMonomial(uint32_t deg)
  {
    _deg = deg;
    generateMonomial();
  }
  //----------------------------------------------------------------------------
  //  Generates a set of monomial expansions
  //  x1 is the point being expanded around.
  //----------------------------------------------------------------------------
  std::vector<double> Monomial::taylorMonomialExpansion(const std::vector<double>& x1,
                                                        const std::vector<double>& x2)
  {
    std::vector<double> taylorExp(x1.size(),0.0);
    if (x2.size() != x1.size() || x1.size() != _dim)
    {
      std::cout << "Dimension of the points (x1,x2) are not compatible with";
      std::cout << " the dimension of the monomial!" << std::endl;
      return taylorExp;
    }
    for (uint32_t i = 0; i < x1.size(); i++)
    {
      taylorExp[i] = pow((x2[i]-x1[i]),1);
    }
    std::vector<double> taylor(_mono.size(),1.0);
    for (uint32_t i = 0; i < _mono.size(); i++)
    {
      for (uint32_t j = 0; j < x1.size(); j++)
      {
        if (_mono[i][j] > 0)
        {
          for (uint32_t k = 0; k < _mono[i][j]; k++)
          {
            taylor[i] *= taylorExp[j];
          }
        }
      }
    }
    return taylor;
  }
  //----------------------------------------------------------------------------
  //  Generates a set of monomial expansions by also
  //  specifying a new degree
  //----------------------------------------------------------------------------
  std::vector<double> Monomial::taylorMonomialExpansion(const std::vector<double>& x1,
                                                        const std::vector<double>& x2,
                                                        uint32_t deg)
  {
    generateMonomial(deg);
    return taylorMonomialExpansion(x1,x2);
  }

  const std::string Monomial::summary()
  {
    std::string s;
    s += "\n----";
    for (uint32_t i = 0; i < _dim; i++)
    {
      s += "-----";
    }
    s += "-----";
    s += "\n   Table of monomials";
    s += "\n     (n,d) = (" + std::to_string(_deg) + ",";
    s += std::to_string(_dim) + ")";
    s += "\n----";
    for (uint32_t i = 0; i < _dim; i++)
    {
      s += "-----";
    }
    s += "-----";
    s += "\nid.  ";
    if(_dim < 10)
    {
      for (uint32_t i = 0; i < _dim; i++)
      {
        s += "x_" + std::to_string(i) + "  ";
      }
      s += "deg";
      s += "\n   +";
      for (uint32_t i = 0; i < _dim; i++)
      {
        s += "-----";
      }
      s += "----+";
      uint32_t currentOrder = 0;
      for (uint32_t i = 0; i < _mono.size(); i++)
      {
        uint32_t order = std::accumulate(_mono[i].begin(),_mono[i].end(),0);
        if (currentOrder < order)
        {
          s += "\n   +";
          for (uint32_t i = 0; i < _dim; i++)
          {
            s += "-----";
          }
          s += "----+";
          currentOrder = order;
        }
        if (i < 10)
        {
          s += "\n" + std::to_string(i) + "  |";
        }
        else
        {
          s += "\n" + std::to_string(i) + " |";
        }
        for (uint32_t j = 0; j < _dim; j++)
        {
          s += "  " + std::to_string(_mono[i][j]) + "  ";
        }
        s += "  " + std::to_string(order) + " |";
      }
      s += "\n----";
      for (uint32_t i = 0; i < _dim; i++)
      {
        s += "-----";
      }
      s += "-----";
    }
    return s;
  }
  //----------------------------------------------------------------------------

  //--------------------------------------------------------------------------
  //  Generates a Taylor polynomial
  //--------------------------------------------------------------------------
  std::vector<double> taylorPolynomial(double p, double x, uint32_t n)
  {
      std::vector<double> taylor(n);
      for (uint32_t i = 1; i < n+1; i++)
      {
          taylor[i-1] = pow((x-p),i);
      }
      return taylor;
  }
  //--------------------------------------------------------------------------

  //--------------------------------------------------------------------------
  //  Checks if an array is consistently defined.
  //--------------------------------------------------------------------------
  bool checkConsistency(std::vector<std::vector<double>>& array,
                        std::vector<std::pair<uint32_t,uint32_t>>& rows)
  {
    //  check that the array is well defined
    bool well_defined = true;
    auto p = std::make_pair(array[0].size(),1);
    rows.push_back(p);
    for (uint32_t i = 0; i < array.size(); i++)
    {
      if (array[0].size() != array[i].size())
      {
        well_defined = false;
        bool unique = false;
        for (uint32_t j = 0; j < rows.size(); j++)
        {
          if (std::get<0>(rows[j]) == array[i].size())
          {
            unique = true;
            std::get<1>(rows[j]) += 1;
          }
        }
        if (unique == false)
        {
          auto p = std::make_pair(array[i].size(),1);
          rows.push_back(p);
        }
      }
    }
    return well_defined;
  }
  //--------------------------------------------------------------------------

  //--------------------------------------------------------------------------
  //  Error messages
  //  Here we have a set of functions for generating various
  //  error messages
  //--------------------------------------------------------------------------
  std::string MATRIX_INCONSISTENT_ARRAY(std::vector<std::pair<uint32_t,uint32_t>>& rows)
  {
    std::string error = "ERROR: Attempted to construct matrix with";
    error += " inconsistent numbers of columns.";
    for (uint32_t i = 0; i < rows.size(); i++)
    {
      error += "  There are " + std::to_string(std::get<1>(rows[i]))
            +  " rows of size " + std::to_string(std::get<0>(rows[i])) + ".";
    }
    return error;
  }
  std::string MATRIX_OUT_OF_BOUNDS(bool axis, const uint32_t& bound,
                                   const uint32_t& attempt,
                                   const std::string& name)
  {
    std::string error = "ERROR: Attempted to access index "
                      + std::to_string(attempt);
    if (name != " ")
    {
      error += " of matrix '" + name + "'";
    }
    error += "; out of bounds for axis " + std::to_string(int(axis))
           + " with size " + std::to_string(bound);
    return error;
  }
  std::string MATRIX_ADD_INCOMPATIBLE_ROWS(const uint32_t& m1,
                                           const uint32_t& m2,
                                           const std::string& name1,
                                           const std::string& name2)
  {
    std::string error = "ERROR: Attempted to add incompatible matrices";
    if (name1 != " " && name2 != " ")
    {
      error += "'" + name1 + "' and '" + name2 + "'";
    }
    error += "; row sizes are m1 = " + std::to_string(m1)
           + " and m2 = " + std::to_string(m2);
    return error;
  }
  std::string MATRIX_ADD_INCOMPATIBLE_COLS(const uint32_t& n1,
                                           const uint32_t& n2,
                                           const std::string& name1,
                                           const std::string& name2)
  {
    std::string error = "ERROR: Attempted to add incompatible matrices";
    if (name1 != " " && name2 != " ")
    {
      error += "'" + name1 + "' and '" + name2 + "'";
    }
    error += "; column sizes are n1 = " + std::to_string(n1)
           + " and n2 = " + std::to_string(n2);
    return error;
  }
  std::string MATRIX_SUB_INCOMPATIBLE_ROWS(const uint32_t& m1,
                                           const uint32_t& m2,
                                           const std::string& name1,
                                           const std::string& name2)
  {
    std::string error = "ERROR: Attempted to subtract incompatible matrices";
    if (name1 != " " && name2 != " ")
    {
      error += "'" + name1 + "' and '" + name2 + "'";
    }
    error += "; row sizes are m1 = " + std::to_string(m1)
           + " and m2 = " + std::to_string(m2);
    return error;
  }
  std::string MATRIX_SUB_INCOMPATIBLE_COLS(const uint32_t& n1,
                                           const uint32_t& n2,
                                           const std::string& name1,
                                           const std::string& name2)
  {
    std::string error = "ERROR: Attempted to subtract incompatible matrices";
    if (name1 != " " && name2 != " ")
    {
      error += "'" + name1 + "' and '" + name2 + "'";
    }
    error += "; column sizes are n1 = " + std::to_string(n1)
           + " and n2 = " + std::to_string(n2);
    return error;
  }
  std::string MATRIX_MUL_INCOMPATIBLE(const uint32_t& n1,
                                      const uint32_t& m2)
  {
    return " ";
  }
  std::string MATRIX_ZERO_DIV(const uint32_t& m, const uint32_t& n)
  {
    return " ";
  }
  //--------------------------------------------------------------------------
}
