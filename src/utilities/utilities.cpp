//------------------------------------------------------------------------------
//  utilities.cpp
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
#include "utilities.h"

namespace ET
{

  size_t UniqueID::nextID = 0;

  //----------------------------------------------------------------------------
  //  Function for turning a number into scientific notation
  //----------------------------------------------------------------------------
  template<typename T>
  std::string scientific_not(T x, size_t dec)
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
    if (exp > 150) {
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
        for (size_t i = 0; i < a.size(); i++) {
            cart[i].resize(b.size());
            for (size_t j = 0; j < b.size(); j++) {
                cart[i][j] = a[i]*b[j];
            }
        }
        return cart;
    }
    //--------------------------------------------------------------------------


  //--------------------------------------------------------------------------
  //  Checks if an array is consistently defined.
  //--------------------------------------------------------------------------
  bool checkConsistency(std::vector<std::vector<double>>& array,
                        std::vector<std::pair<size_t,size_t>>& rows)
  {
    //  check that the array is well defined
    bool well_defined = true;
    auto p = std::make_pair(array[0].size(),1);
    rows.push_back(p);
    for (size_t i = 0; i < array.size(); i++) {
      if (array[0].size() != array[i].size()) {
        well_defined = false;
        bool unique = false;
        for (size_t j = 0; j < rows.size(); j++) {
          if (std::get<0>(rows[j]) == array[i].size()) {
            unique = true;
            std::get<1>(rows[j]) += 1;
          }
        }
        if (unique == false) {
          auto p = std::make_pair(array[i].size(),1);
          rows.push_back(p);
        }
      }
    }
    return well_defined;
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
  }
  Monomial::Monomial(size_t t_dim)
  : m_dim(t_dim)
  {
  }
  Monomial::Monomial(size_t t_dim, size_t t_deg)
  : m_dim(t_dim), m_deg(t_deg)
  {
    generateMonomial();
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  Getters and Setters
  //----------------------------------------------------------------------------
  size_t Monomial::getDim()
  {
    return m_dim;
  }
  size_t Monomial::getDeg()
  {
    return m_deg;
  }
  std::vector<std::vector<size_t>> Monomial::getMono()
  {
    return m_mono;
  }
  std::vector<std::vector<size_t>> Monomial::getMono(size_t t_deg)
  {
    generateMonomial(t_deg);
    return m_mono;
  }
  std::vector<std::vector<std::string>> Monomial::getMonoFactors()
  {
    return m_monoFactors;
  }
  size_t Monomial::getMultisetCoefficient(size_t t_deg)
  {
    //  use the i4_choose function from monomial.hpp
    //  the degree of the kth-monomial for a space of
    //  dimension d is the binomial coefficent (k + d - 1 // k)
    return i4_choose(t_deg + m_dim - 1,t_deg);
  }
  size_t Monomial::getTaylorIndex(std::vector<size_t> t_term)
  {
    if (t_term.size() != m_dim) {
      std::cout << "Number of elements in exp.";
      std::cout << " array does not equal dimension!" << std::endl;
      return 0;
    }
    std::vector<int> x(m_dim);
    std::vector<int> y(m_dim);
    for (size_t i = 0; i < m_dim; i++) {
      x[i] = 0;
      y[i] = t_term[i];
    }
    int index = 0;
    while (x != y) {
      mono_between_next_grlex(m_dim,0,m_deg,x.data());
      index++;
    }
    return index;
  }
  void Monomial::setDim(size_t t_dim)
  {
    m_dim = t_dim;
  }
  void Monomial::setDeg(size_t t_deg)
  {
    m_deg = t_deg;
  }
  void Monomial::generateMonomial()
  {
    int x[m_dim];
    for (size_t i = 0; i < m_dim; i++) {
      x[i] = 0;
    }
    std::vector<std::vector<size_t>> mono;
    std::vector<size_t> temp(x, x + m_dim);
    mono.push_back(temp);
    int i = 1;
    for ( ; ; ) {
      if ( x[0] == m_deg ) {
        break;
      }
      mono_between_next_grlex(m_dim, 0, m_deg, x);
      std::vector<size_t> temp(x, x + m_dim);
      mono.push_back(temp);
      i = i + 1;
    }
    m_mono = mono;
  }
  void Monomial::generateMonomial(size_t t_deg)
  {
    m_deg = t_deg;
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
    if (x2.size() != x1.size() || x1.size() != m_dim) {
      std::cout << "Dimension of the points (x1,x2) are not compatible with";
      std::cout << " the dimension of the monomial!" << std::endl;
      return taylorExp;
    }
    for (size_t i = 0; i < x1.size(); i++) {
      taylorExp[i] = pow((x2[i]-x1[i]),1);
    }
    std::vector<double> taylor(m_mono.size(),1.0);
    for (size_t i = 0; i < m_mono.size(); i++) {
      for (size_t j = 0; j < x1.size(); j++) {
        if (m_mono[i][j] > 0) {
          for (size_t k = 0; k < m_mono[i][j]; k++) {
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
                                                        size_t deg)
  {
    generateMonomial(deg);
    return taylorMonomialExpansion(x1,x2);
  }

  const std::string Monomial::summary()
  {
    std::string s;
    s += "\n----";
    for (size_t i = 0; i < m_dim; i++) {
      s += "-----";
    }
    s += "-----";
    s += "\n   Table of monomials";
    s += "\n     (n,d) = (" + std::to_string(m_deg) + ",";
    s += std::to_string(m_dim) + ")";
    s += "\n----";
    for (size_t i = 0; i < m_dim; i++) {
      s += "-----";
    }
    s += "-----";
    s += "\nid.  ";
    if(m_dim < 10) {
      for (size_t i = 0; i < m_dim; i++) {
        s += "x_" + std::to_string(i) + "  ";
      }
      s += "deg";
      s += "\n   +";
      for (size_t i = 0; i < m_dim; i++) {
        s += "-----";
      }
      s += "----+";
      size_t currentOrder = 0;
      for (size_t i = 0; i < m_mono.size(); i++) {
        size_t order = std::accumulate(m_mono[i].begin(),m_mono[i].end(),0);
        if (currentOrder < order) {
          s += "\n   +";
          for (size_t i = 0; i < m_dim; i++) {
            s += "-----";
          }
          s += "----+";
          currentOrder = order;
        }
        if (i < 10) {
          s += "\n" + std::to_string(i) + "  |";
        }
        else {
          s += "\n" + std::to_string(i) + " |";
        }
        for (size_t j = 0; j < m_dim; j++) {
          s += "  " + std::to_string(m_mono[i][j]) + "  ";
        }
        s += "  " + std::to_string(order) + " |";
      }
      s += "\n----";
      for (size_t i = 0; i < m_dim; i++) {
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
  std::vector<double> taylorPolynomial(double p, double x, size_t n)
  {
      std::vector<double> taylor(n);
      for (size_t i = 1; i < n+1; i++) {
          taylor[i-1] = pow((x-p),i);
      }
      return taylor;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //  UniqueID Class
  //----------------------------------------------------------------------------
  UniqueID::UniqueID()
  {
    id = ++nextID;
  }
  //----------------------------------------------------------------------------
  UniqueID::UniqueID(const UniqueID& orig)
  {
    id = orig.id;
  }
  //----------------------------------------------------------------------------
  UniqueID& UniqueID::operator=(const UniqueID& orig)
  {
    id = orig.id;
    return(*this);
  }
  //----------------------------------------------------------------------------
  template std::string scientific_not<double> (double, size_t);

}
