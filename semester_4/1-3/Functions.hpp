#ifndef VEC_MTR_FUNCS_HPP
#define VEC_MTR_FUNCS_HPP

#include <stdexcept>
#include <vector>
#include <random>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>


using namespace std;
using vInt = vector<int>;
using Vec = vector<double>;
using Mtr = vector<Vec>;

Mtr upperTriang(const Mtr& A);
Mtr transpose(const Mtr& matrix);

Vec scalarMultiply(double scalar, const Vec& vec);
Vec vecSum(double c1, const Vec& vec1, double c2, const Vec& vec2);

double dotProduct(const Vec& vec1, const Vec& vec2);
double euclideanNorm(const Vec& vec);
double mean(vInt V);

Vec multiplyMatrixVector(const Mtr& A, const Vec& x);
Mtr matMul(const Mtr& A, const Mtr& B);

Vec generateRandomVector(int n, double lower=0, double upper=1);
Mtr generateRandomMatrix(int rows, int cols=-1, double lower=0, double upper=10);
Mtr generateRndSymPos(int n, double cond);
Mtr generateRndSymPos(int n);

void print(const Vec& v);
void print(const Mtr& A);
void print(const Mtr& A, const Vec& x, const Vec& b);

pair<Vec, int> conjugateGradientMethod(const Mtr& A, const Vec& b, double eps, int maxIter=1000, bool verbose=false);
void PCG(const Mtr& A, const Vec& b, int lowestDeg=-15, int maxIter=1000, bool verbose=false);

double PCG_multi(int n_repeats, int mtr_size, double eps, int maxIter);
void PCG_multi_convergence(int n_repeats, int maxIter, double cond=0);

#endif