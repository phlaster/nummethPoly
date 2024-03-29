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
#include <math.h>
#include <tuple>
#include <utility>


using namespace std;
using Vec = vector<double>;
using vInt = vector<int>;
using Mtr = vector<Vec>;
using LU_result = tuple <Mtr, Mtr, vInt>;


// Boolean checks
bool issquare(const Mtr& M);
bool isapprox(const Vec& V1, const Vec& V2, double tol=1e-6);
bool isapprox(const Mtr& M1, const Mtr& M2, double tol=1e-6);
bool issingleval(const Mtr& M);


// Generators
Vec fill(double val, size_t length);
Mtr fill(double val, size_t rows, size_t cols);
Mtr E(size_t size);
Mtr upperTrSymmetric(const Mtr& A);
Mtr T(const Mtr& M);
Vec normalize(const Vec& V);


// Special
Mtr inv(const Mtr& M);
Mtr toCol(const Vec& v);
Mtr toRow(const Vec& v);
Mtr householder(const Vec& V);


// Random generators
Vec randVec(int n, double lower=0, double upper=1);
Mtr randMtr(int rows, int cols=-1, double lower=0, double upper=1);
Mtr randSymmPositive(int n);


// Matrix algebra
Vec sum(const Vec& V1, const Vec& V2, double c1=1, double c2=1);
Mtr sum(const Mtr& M1, const Mtr& M2, double c1=1, double c2=1);
Vec mul(double scalar, const Vec& V);
Mtr mul(double scalar, const Mtr& M);
Vec mul(const Mtr& M, const Vec& x);
Mtr mul(const Mtr& M1, const Mtr& M2);


// Measures
double dot(const Vec& V1, const Vec& V2);
double euclideanNorm(const Vec& V);
// double euclideanNorm(const Mtr& M);
double det(const Mtr& M);
double maxdiff(const Mtr& m1, const Mtr& m2);
double cond(const Mtr& M, double p=2);


// LUP
LU_result LUP(Mtr M);
void find_max_and_swap(Mtr& M, vInt& perm, int i);
void subtract_product(Mtr& M, int i);
vInt inversePermutation(const vInt& permutation);
Mtr permutationMatrix(const vInt& permutation);

// SLAE stuff
Vec solveLinearEquation(const Mtr& L,
                        const Mtr& U,
                        const vInt& permutation,
                        const Vec& b);
Vec residual(const Mtr& A, const Vec& x, const Vec& b);


// Printers
void print(const Vec& v);
void print(const vInt& v);
void print(const Mtr& A);
void print(const Mtr& A, const Vec& x, const Vec& b);
void print(const Mtr& A, const vector<string>& S, const Vec& b);
void print(const double value);
void print(const bool b);
void print(const char* s);
void print(const LU_result& res);

#endif