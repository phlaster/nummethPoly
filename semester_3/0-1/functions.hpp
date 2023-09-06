#ifndef FUNCS_HPP
#define FUNCS_HPP

#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>

using namespace std;
using Vec = vector<double>;

struct ans {
    double root;
    double absErr;
    int nSteps;
};

double f1(double x);
double f1d1(double x);
double f1d2(double x);

double f2(double x);
double f2d1(double x);
double f2d2(double x);

int sign(double x);
ans bisection(double (*f)(double), double exact, Vec lims, double eps);
ans newton(double (*f)(double), double exact, Vec lims, double eps);
ans newton_broken(double exact, Vec lims, double x0, double eps);

void errorConverg(double exact, Vec interval,
                  ans    (*method)(double (*f)(double), double, Vec, double),
                  double (*f)     (double),
                  string filename);
void brokenConverg(double exact, double x0, Vec lims, string filename);
#endif