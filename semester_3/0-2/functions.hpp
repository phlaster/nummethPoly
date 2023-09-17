#ifndef FUNCS_HPP
#define FUNCS_HPP

#include <iostream>
#include <iomanip>
#include <cassert>
#include <string>
#include <cfloat>
#include <vector>
#include <fstream>
#include <random>

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

struct Buffer {
    Str data;

    void append(const Str& str);
    void append(int value);
    void append(double value, bool last=false);
private:
    static Str doubleToString(double value, int digits=5);
};


double f1(double x, bool derivative=false);
double f2(double x, bool derivative=false);
double sum(const Vec& v);
double mean(const Vec& v);
double amplitude(Vec vector);

Vec chebyshevGrid(double a, double b, int n);
Vec uniformGrid(double a, double b, int n);
Vec errorProfile(const Vec& v1, const Vec& v2);

Graphic tabulateFunction(double (*f)(double, bool), const Vec& grid);
Graphic tabulateDerivative(double (*f)(double, bool), const Graphic& nodes);
Graphic tabulateDerivativeNum(const Graphic& main, int LagrangePoints=3);

double lagrangeTerm(double x, const Vec& x_i, const Vec& y_i);
Graphic lagrangeInterpol(const Graphic& nodes, const Vec& grid);

double hermiteTerm(double a, double b, double fa, double fb, double dfa, double dfb, double t);
Graphic hermiteSpline(const Graphic& function, const Graphic& derivative, int N);

Graphic deviate(const Graphic& g, double sigma);


void task2_3(double (*f)(double, bool), int minNodes, int maxNodes, int coeffMult, Vec lims, Str threadname, Buffer& buffer);
void task4(double (*f)(double, bool), int minNodes, int maxNodes, int coeffMult, Vec lims, Str threadname, Buffer& buffer);
void task5(double (*f)(double, bool), int minNodes, int maxNodes, int coeffMult, double deviation, Vec lims, Str threadname, Buffer& buffer);
void task5progression(double (*f)(double, bool), int nNodes, int coeffMult, Vec lims, Str threadname, Buffer& buffer);

void write2_3(
    const Str& filename,
    const Graphic& uniformNodes,
    const Graphic& uniformNodes_DerNum,

    const Graphic& uniformTab,
    const Graphic& uniformTab_Der,

    const Graphic& uniformLagr,
    const Graphic& uniformHerm,

    const Vec& err_uniformLagr,
    const Vec& err_uniformHerm
);

void write4(
    const Str& filename,
    const Graphic& uniformNodes,
    const Graphic& chebNodes,

    const Graphic& uniformTab,

    const Graphic& uniformLagr,
    const Graphic& chebLagr,

    const Vec& err_uniformLagr,
    const Vec& err_chebLagr
);

void write5(
    const Str& filename,
    const Graphic& uniformNodes_Noize,
    const Graphic& chebNodes_Noize,
    const Graphic& uniformNodes_Noize_DerNum,

    const Graphic& uniformTab,
    const Graphic& uniformTab_Der,

    const Graphic& uniformLagr_Noize,
    const Graphic& uniformHerm_Noize,
    const Graphic& chebLagr_Noize,

    const Vec& err_uniformLagr_Noize,
    const Vec& err_uniformHerm_Noize,
    const Vec& err_chebLagr_Noize
);

#endif