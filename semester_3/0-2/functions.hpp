#ifndef FUNCS_HPP
#define FUNCS_HPP

#include <iostream>
#include <utility>
// #include <string>
#include <cfloat>
#include <cmath>
#include <math.h>
#include <vector>

using namespace std;
using Vec = vector<double>;
using Coords = pair<Vec, Vec>;

struct Graphic
{
    Vec xVals;
    Vec yVals;
    double dx;
    int N;
};


double cot(double x);
double f1(double x, bool derivative=false);
double f2(double x, bool derivative=false);

double F1(double a, double b);
double F2(double a, double b);

Graphic calculateGraphic(double (*f)(double, bool), double left, double right, int N);
Graphic calculateGraphic(double (*f)(double, bool), double left, double right, double dx);

Graphic calculateDerivativeAnalitical(double (*f)(double, bool), const Vec& x_grid);
Graphic calculateDerivativeNumerical(double (*f)(double, bool), const Graphic& main, int LagrangeDegree=3);

double lagrange(double x, const Vec& x_i, const Vec& y_i);
double hermitePolynomial(double a, double b, double fa, double fb, double dfa, double dfb, double t);
Graphic hermiteSpline(const Graphic& function, const Graphic& derivative, double dt);


#endif