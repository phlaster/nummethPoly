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

vector<point> RK2_adaptive(
    const double y0,
    const double a,
    const double b,
    const double h,
    const double p,
    double (*func)(double, double)
);

double RK2_maxerror(
    const double y0,
    const double a,
    const double b,
    const double h,
    double (*func)(double, double),
    double (*etalon)(double)
);
#endif