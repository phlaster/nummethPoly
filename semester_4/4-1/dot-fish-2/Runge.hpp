#pragma once
#include <iostream>
#include <cmath>
#include <vector>

typedef struct point
//struct to avoid having 2 dimentional arrays/vectors
{
    double x;
    double y;
} point;

std::vector<point> RK2(double y0, double a, double b, double h, double (*func)(double, double));
std::vector<point> RK2_adaptive(double y0, double a, double b, double h, double p, double (*func)(double, double));
double RK2_maxerror(double y0, double a, double b, double h, double (*func)(double, double), double (*etalon)(double));