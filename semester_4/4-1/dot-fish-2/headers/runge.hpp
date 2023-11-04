#ifndef RUNGE_HPP
#define RUNGE_HPP
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>
using namespace std;

//struct to avoid having 2 dimentional arrays/vectors
typedef struct point {
    double x;
    double y;
} point;

vector<point> RK2(
    const double y0,
    const double a,
    const double b,
    const double h,
    double (*func)(double, double)
);

tuple<double, double, int> RK2_adaptive(
    const double y0,
    const double a,
    const double b,
    const double h,
    const double eps,
    double (*func)(double, double),
    double (*etalon)(double)
);

tuple<double, double, int> RK2_maxerror(
    const double y0,
    const double a,
    const double b,
    const double h,
    double (*func)(double, double),
    double (*etalon)(double)
);

pair<double, int> RK2_recursive(
    double (*f)(double, double),
    double (*exact)(double),
    const double y0,
    const double x0,
    const double b,
    const double eps
);

pair<double, int> RK2_adpt(
    double (*f)(double, double),
    double (*exact)(double),
    const double y0,
    const double x0,
    const double b,
    const double h,
    const double eps
);

pair<double, int> RK2___(
    const double y0,
    const double x0,
    const double b,
    const double h0,
    const double eps,
    double (*f)(double, double),
    double (*exact)(double)
);

#endif