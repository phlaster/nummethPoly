#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP
#include "Runge.hpp"

double exact_y(double x);
double dydx(double x, double y);
double error_max(vector<point> &vec);
int max_error_index(vector<point> &vec);
#endif