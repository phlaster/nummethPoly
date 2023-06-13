/*
    вариант 8 Б
    
    x^4 - x^3 - 2x^2 + 3x - 3 = 0

    3exp(x) = 5x + 3

*/

#include <iostream>
#include <cmath>
using namespace std;

double p4(double x,
          double a4 = 1.,
          double a3 = -1.,
          double a2 = -2.,
          double a1 = 3.,
          double a0 = -3.)
{
    return a4 * pow(x, 4) +
           a3 * pow(x, 3) +
           a2 * pow(x, 2) +
           a1 * x +
           a0;
}


int main()
{
    // Применение теоремы о верхней границы
    return 0;
}