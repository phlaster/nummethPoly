#ifndef FUNCS_HPP
#define FUNCS_HPP

#include <iostream>
#include <utility>
#include <string>
#include <cfloat>
#include <cmath>
#include <math.h>
#include <vector>

using namespace std;
using Vec = vector<double>;
using Coords = pair<Vec, Vec>;
using Str = string;

struct Graphic
{
    Vec xVals;
    Vec yVals;
    double dx;
    int N;
};


double ctg(double x);
double f1(double x, bool derivative=false);
double f2(double x, bool derivative=false);

double F1(double a, double b);
double F2(double a, double b);

Graphic tabulateFunction(double (*f)(double, bool), double left, double right, int N);
Graphic tabulateFunction(double (*f)(double, bool), double left, double right, double dx);

Graphic tabulateDerivative(double (*f)(double, bool), const Graphic& nodes);
Graphic tabulateDerivativeNum(double (*f)(double, bool), const Graphic& main, int LagrangePoints=3);

double lagrangeTerm(double x, const Vec& x_i, const Vec& y_i);
Graphic lagrangeInterpol(const Graphic& func, double left, double right, int N);
Graphic lagrangeInterpol(const Graphic& func, double left, double right, double dx);


double hermiteTerm(double a, double b, double fa, double fb, double dfa, double dfb, double t);
Graphic hermiteSpline(const Graphic& function, const Graphic& derivative, int N);
Graphic hermiteSpline(const Graphic& function, const Graphic& derivative, double dt);

#endif