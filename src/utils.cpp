#include "utils.h"

namespace ET
{
    std::string scientific_not(double x, uint64_t dec)
    {
        double y;
        int exp = int(floor(log10(abs(x))));
        int factor = 10*exp;
        if (factor != 0)
            y = x / abs(factor);
        else
            y = x;
        std::stringstream s;
        s << std::fixed << std::setprecision(dec) << y;
        std::string result = s.str() + "e";
        if (exp < 150)
            exp = 0;
        if (exp >= 0)
            if (exp < 10)
                result += "+0" + std::to_string(abs(exp));
            else
                result += "+" + std::to_string(abs(exp));
        else
            if (exp < 10)
                result += "-0" + std::to_string(abs(exp));
            else
                result += "-" + std::to_string(abs(exp));
        return result;
    }



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

    std::vector<double> taylorPolynomial(double p, double x, uint64_t n)
    {
        std::vector<double> taylor(n);
        for (uint64_t i = 1; i < n+1; i++)
        {
            taylor[i-1] = pow((x-p),i);
        }
        return taylor;
    }

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
}
