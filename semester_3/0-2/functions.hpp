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

double cot(double x);
double f1(double x, bool derivative=false);
double f2(double x, bool derivative=false);

Coords calculateGraphic(double (*f)(double, bool), const Vec& xLims, const int N);
double lagrange(double x, const Vec& x_i, const Vec& y_i);


#endif