#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP
#include "runge.hpp"

double exact_y(const double x);
double dydx(const double x, double y);
double error_max(const vector<point> &vec);
int max_error_index(const vector<point> &vec);
#endif