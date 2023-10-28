#ifndef HEADER_HPP
#define HEADER_HPP

#include "Types.hpp"
#include "Printers.hpp"
#include "operators.hpp"
#include "Sparse.hpp"

Vec randVec(int n, double lower=0, double upper=1);

double euclideanNorm(const Vec& V);
double euclideanNorm(const spMtr& M);
Vec normalize(const Vec& V);
spMtr normalize(const spMtr& M);
Mtr toCol(const Vec& v);
Mtr toRow(const Vec& v);
Mtr householder(const Vec& V);
Mtr generateRndSymPos(int n, double cond);







Mtr E(size_t size);
Mtr T(const Mtr& M);

Mtr to_dense(const spMtr& SP);
spMtr T(const spMtr& M);
spMtr E(size_t size, bool sparse);
spMtr householder(const spMtr& M);
spMtr generateRndSymPos(int n, double cond, double sparsity);

void erase_above_diag(spMtr& A, bool below=false);

// void choleskyDecomposition(Mtr& A);
void chol(spMtr& A);
void incomp_chol_zero_tol(spMtr& A, const double threshold=1e-16);
void ___incomp_chol_tol(spMtr& A, const double threshold);
#endif