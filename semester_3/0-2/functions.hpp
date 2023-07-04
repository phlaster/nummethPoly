#ifndef FUNCS_HPP
#define FUNCS_HPP

#include <iostream>
#include <cassert>
#include <utility>
#include <algorithm>
#include <string>
#include <cfloat>
#include <cmath>
#include <math.h>
#include <vector>
#include <ios>
#include <fstream>

using namespace std;
using Vec = vector<double>;
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

Vec errorProfile(const Vec& v1, const Vec& v2);
double mean(const Vec& v);

Graphic tabulateFunction(double (*f)(double, bool), double left, double right, int N);
Graphic tabulateFunction(double (*f)(double, bool), double left, double right, double dx);
Graphic tabulateFunction(double (*f)(double, bool), const Vec& grid);
Vec chebyshevGrid(double a, double b, int n);

Graphic tabulateDerivative(double (*f)(double, bool), const Graphic& nodes);
Graphic tabulateDerivativeNum(double (*f)(double, bool), const Graphic& main, int LagrangePoints=3);

double lagrangeTerm(double x, const Vec& x_i, const Vec& y_i);
Graphic lagrangeInterpol(const Graphic& exactNodes, int N);
// Graphic lagrangeInterpol(const Graphic& func, double left, double right, double dx);


double hermiteTerm(double a, double b, double fa, double fb, double dfa, double dfb, double t);
Graphic hermiteSpline(const Graphic& function, const Graphic& derivative, int N);
Graphic hermiteSpline(const Graphic& function, const Graphic& derivative, double dt);

void task2Print(Graphic exactNodes, Graphic exactFunction, Graphic exactDerivative, Graphic derivativeNum, Graphic interpolL, Graphic interpolH, Str filename);
void task3Print(Graphic evenlySpacedNodes, Graphic chebyshevNodes, Graphic exactForEvenly, Graphic exactForChebyshev, Graphic evenlyL, Graphic chebyshevL, int nNodes, int nInterpol, Str filename);
void print(
    const Str& filename,
    const Graphic& uniformNodes,
    const Graphic& uniformNodes_DerNum,
    const Graphic& chebNodes,
    const Graphic& uniformTab,
    const Graphic& uniformTab_Der,
    const Graphic& chebTab,
    const Graphic& uniformLagr,
    const Graphic& uniformHerm,
    const Graphic& chebLagr,
    const Vec& err_uniformLagr,
    const Vec& err_uniformHerm,
    const Vec& err_chebLagr
);
#endif