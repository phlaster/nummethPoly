#ifndef VEC_MTR_FUNCS_HPP
#define VEC_MTR_FUNCS_HPP

#include <array>
#include <stdexcept>
#include <vector>
#include <random>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>

using namespace std;

using Vec = vector<double>;
using Mtr = vector<Vec>;
using LU_result = tuple <Mtr, Mtr, vector<int>>;


bool isapprox(const Vec& V1, const Vec& V2, double tol=1e-6);
bool isapprox(const Mtr& M1, const Mtr& M2, double tol=1e-6);

Vec fill(double val, size_t length);
Mtr fill(double val, size_t rows, size_t cols);
Mtr E(size_t size);
Mtr diag(const Vec& V);

Mtr upperTrSymmetric(const Mtr& A);
Mtr T(const Mtr& M);

double max(const Vec& v);
double min(const Vec& v);
Vec mul(double scalar, const Vec& V);
Vec sum(const Vec& V1, const Vec& V2, double c1=1, double c2=1);
Mtr sum(const Mtr& M1, const Mtr& M2, double c1=1, double c2=1);
double sum(const Vec& v);
double mean(const Vec& v);
Vec div(const Vec& v1, const Vec& v2, double eps=1e-16);

double dot(const Vec& V1, const Vec& V2);
double euclideanNorm(const Vec& V);
// double euclideanNorm(const Mtr& M);
double opnorm_1(const Mtr& M);
double det(const Mtr& M);
double cond(const Mtr& M);

double maxdiff(const Mtr& m1, const Mtr& m2);
double maxdiff(const Vec& v1, const Vec& v2);

Mtr inv(const Mtr& M);
Vec normalize(const Vec& V);

Mtr mul(double scalar, const Mtr& M);
Vec mul(const Mtr& M, const Vec& x);
Mtr mul(const Mtr& M1, const Mtr& M2);

Vec randVec(int n, double lower=0, double upper=1);
Mtr randMtr(int rows, int cols=-1, double lower=0, double upper=1);
Mtr randSymmPositive(int n);
Mtr generateRndSymPos(int n, double cond);

Mtr toCol(const Vec& v);
Mtr toRow(const Vec& v);

Mtr householder(const Vec& V);

LU_result LU_decomposition(Mtr M);
Mtr apply_row_permutation(const Mtr& M, const vector<int>& perm);
Vec solveLinearEquation(const Mtr& L,
                        const Mtr& U,
                        const std::vector<int>& permutation,
                        const Vec& b);
Vec residual(const Mtr& A, const Vec& x, const Vec& b);

pair<int, double> naive_PM(const Mtr& A, double eps=1e-7);
pair<int, double> normed_PM(const Mtr& A, double eps=1e-7);





void print(const Vec& v);
void print(const vector<int>& v);
void print(const Mtr& A);
void print(const Mtr& A, const Vec& x, const Vec& b);
void print(const double value);
void print(const bool b);
void print(const char* s);
void print(const LU_result& res);



#endif