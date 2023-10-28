#ifndef RUNGE_HPP
#define RUNGE_HPP
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>
using namespace std;

typedef struct point
//struct to avoid having 2 dimentional arrays/vectors
{
    double x;
    double y;
} point;

vector<point> RK2(
    double y0, double a, double b, double h, double (*func)(double, double)
);

vector<point> RK2_adaptive(
    double y0, double a, double b, double h, double p, double (*func)(double, double)
);

double RK2_maxerror(
    double y0, double a, double b, double h, double (*func)(double, double), double (*etalon)(double)
);
#endif