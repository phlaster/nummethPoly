#include "functions.hpp"
#include <cmath>
#include <cstddef>

double sum(const Vec& v)
{
    double S = 0;
    for (auto i : v)
        S += i;
    return S;
}

double mean(const Vec& v)
{
    return sum(v)/v.size();
}

double amplitude(Vec vector)
{
    double min = INFINITY;
    double max = -INFINITY;
    for (auto n: vector)
    {
        if (n > max) max = n;
        if (n < min) min = n;
    }
    return (max - min)/2.0;
}

// Профиль ошибки: абсолютные ошибки каждого компонента
Vec errorProfile(const Vec& v1, const Vec& v2)
{
    assert(v1.size()==v2.size() && "Для вычисления ошибки векторы должны иметь одинаковую длину!\n");
    Vec errs(v1.size());
    for (size_t i=0; i<v1.size(); i++)
        errs[i] = fabs(v1[i]-v2[i]);
    return errs;
}


double f1(double x, bool derivative){
    return !derivative ? x*x + cos(x)/sin(x) :
                         x+x - 1/pow(sin(x),2);
}

double f2(double x, bool derivative){
    return !derivative ? pow(x,5) - 3.2*pow(x,3) + 2.5*x*x - 7*x + 1.5 :
                       5*pow(x,4) - 9.6*pow(x,2) +   5*x   - 7;
}
