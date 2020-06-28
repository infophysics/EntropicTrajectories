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
    //--------------------------------------------------------------------------
    //  Function for turning a number into scientific notation
    //--------------------------------------------------------------------------
    std::string scientific_not(double x, uint64_t dec)
    {
        double y;
        int exp = int(floor(log10(abs(x))));
        double factor = pow(10,exp);
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
        for (uint64_t i = 0; i < a.size(); i++)
        {
            cart[i].resize(b.size());
            for (uint64_t j = 0; j < b.size(); j++)
            {
                cart[i][j] = a[i]*b[j];
            }
        }
        return cart;
    }
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Function for generating a set of monomial indices
    //--------------------------------------------------------------------------
    std::vector<std::vector<uint64_t> > monomial_n(uint64_t dim, uint64_t n)
    {
        int x[dim];
        for (uint64_t i = 0; i < dim-1; i++)
        {
            x[i] = 0;
        }
        x[dim-1] = 1;
        std::vector<std::vector<uint64_t> > mono;
        std::vector<uint64_t> temp(x, x + dim);
        mono.push_back(temp);
        int i = 1;
        for ( ; ; )
        {
            if ( x[0] == n )
            {
              break;
            }
            mono_between_next_grlex(dim, 0, n, x);
            std::vector<uint64_t> temp(x, x + dim);
            mono.push_back(temp);
            i = i + 1;
        }
        return mono;
    }
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Generates a Taylor polynomial
    //--------------------------------------------------------------------------
    std::vector<double> taylorPolynomial(double p, double x, uint64_t n)
    {
        std::vector<double> taylor(n);
        for (uint64_t i = 1; i < n+1; i++)
        {
            taylor[i-1] = pow((x-p),i);
        }
        return taylor;
    }
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Generates a set of monomial expansions
    //--------------------------------------------------------------------------
    std::vector<double> taylorMonomialExpansion(std::vector<double> p,
        std::vector<double> x, uint64_t n)
    {
        std::vector<std::vector<uint64_t> > mono = monomial_n(p.size(), n);
        std::vector<double> taylor_exps(p.size());
        for (uint64_t i = 0; i < p.size(); i++)
        {
            taylor_exps[i] = pow((x[i]-p[i]),1);
        }
        std::vector<double> taylor(mono.size(),1.0);
        for (uint64_t i = 0; i < mono.size(); i++)
        {
            for (uint64_t j = 0; j < p.size(); j++)
            {
                if (mono[i][j] > 0)
                {
                    for (uint64_t k = 0; k < mono[i][j]; k++)
                    {
                        taylor[i] *= taylor_exps[j];
                    }
                }
            }
        }
        return taylor;
    }
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //  Generates a set of monomial factors of order n
    //--------------------------------------------------------------------------
    std::vector<std::string> taylorMonomialFactors(uint64_t dim, uint64_t n)
    {
        std::vector<std::vector<uint64_t> > mono = monomial_n(dim, n);
        std::vector<std::string> xs(dim);
        for (uint64_t i = 0; i < dim; i++)
        {
            xs[i] = "x_" + std::to_string(i);
        }
        std::vector<std::string> taylor(mono.size(),"");
        for (uint64_t i = 0; i < mono.size(); i++)
        {
            for (uint64_t j = 0; j < dim; j++)
            {
                if (mono[i][j] > 0)
                {
                    taylor[i] += xs[j];
                    if (mono[i][j] > 1)
                    {
                        taylor[i] += "^" + std::to_string(mono[i][j]);
                    }
                }
            }
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
      std::string error = "ERROR: Attempted to multiply incompatible matrices";
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
      std::string error = "ERROR: Attempted to multiply incompatible matrices";
      if (name1 != " " && name2 != " ")
      {
        error += "'" + name1 + "' and '" + name2 + "'";
      }
      error += "; column sizes are n1 = " + std::to_string(n1)
             + " and n2 = " + std::to_string(n2);
      return error;
    }
    std::string MATRIX_SUB_INCOMPATIBLE_ROWS(const uint32_t& m1,
                                             const uint32_t& m2)
    {
    }
    std::string MATRIX_SUB_INCOMPATIBLE_COLS(const uint32_t& n1,
                                             const uint32_t& n2)
    {
    }
    std::string MATRIX_MUL_INCOMPATIBLE(const uint32_t& n1,
                                        const uint32_t& m2)
    {
    }
    std::string MATRIX_ZERO_DIV(const uint32_t& m, const uint32_t& n)
    {
    }

    //--------------------------------------------------------------------------
}
