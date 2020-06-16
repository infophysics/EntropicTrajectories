#include "utils.h"

namespace ET
{
    std::string scientific_not(double x, unsigned int dec)
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
}
