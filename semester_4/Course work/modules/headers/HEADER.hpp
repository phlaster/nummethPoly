#ifndef HEADER_HPP
#define HEADER_HPP

#include "Types.hpp"
#include "Printers.hpp"
#include "operators.hpp"
#include "Sparse.hpp"
#include <cmath>

Vec randVec(size_t n, double lower=0, double upper=1);

double euclideanNorm(const Vec& V);
double euclideanNorm(const spMtr& M);
Vec normalize(const Vec& V);
spMtr normalize(const spMtr& M);
Mtr toCol(const Vec& v);
Mtr toRow(const Vec& v);
Mtr householder(const Vec& V);
Mtr generateRndSymPos(size_t n, double cond);

Mtr E(size_t size);
Mtr T(const Mtr& M);

Mtr to_dense(const spMtr& SP);

spMtr T(const spMtr& M);
spMtr E(size_t size, bool sparse);
spMtr householder(const spMtr& M);
spMtr generateRndSymPos(int n, double cond, double density);
Mtr generateRndSymPos(const Vec& V, const Vec& conds);
pair<double, double> mean_std(const Mtr& A);

double maximum(const spMtr& A);
double maxabs(const spMtr& A);
pair<double, double> mean_std(const spMtr& A);
spMtr erase_above_diag(spMtr A, bool below=false);
spMtr sparsen(spMtr M, double prob);
spMtr chol(const spMtr& A);
spMtr ichol(const spMtr& A, double theta=0.0);


Vec solve_L(const spMtr& L, const Vec& b);
Vec solve_U(const spMtr& U, const Vec& b);
Vec solve_L_U(const spMtr& L, const spMtr& U, const Vec& b);

void spy(const spMtr& sA);

pair<Vec, int> cg(
    const spMtr& A,
    const Vec& b,
    double eps=1e-16,
    int maxIter=500
);
pair<Vec, int> pcg(
    const spMtr& A,
    const spMtr& L, // Cholesky reconditioner
    const Vec& b,
    double eps=1e-16,
    int maxIter=500
);

spMtr block5diag(size_t m, size_t n);
spMtr block5diag(size_t n);
Vec estimate_b(size_t n);
Vec exact_x(double n);
#endif