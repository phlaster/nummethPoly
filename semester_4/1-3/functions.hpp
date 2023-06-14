#ifndef VEC_MTR_FUNCS_HPP
#define VEC_MTR_FUNCS_HPP

#include "Types.hpp"

Matrix upperTriang(const Matrix& A);
Matrix transpose(const Matrix& matrix);

Vector scalarMultiply(double scalar, const Vector& vec);
Vector vecSum(double c1, const Vector& vec1, double c2, const Vector& vec2);

double dotProduct(const Vector& vec1, const Vector& vec2);
double euclideanNorm(const Vector& vec);

Vector multiplyMatrixVector(const Matrix& A, const Vector& x);
Matrix matMul(const Matrix& A, const Matrix& B);

Vector generateRandomVector(int n, double lower=0, double upper=10);
Matrix generateRandomMatrix(int rows, int cols=-1, double lower=0, double upper=10);
Matrix generateRndSymPos(int n);

int conjugateGradientMethod(const Matrix& A, const Vector& b, double eps, int maxIter=1000);

void printMatrix(const Matrix& A);
void printSLAE(const Matrix& A, const Vector& x, const Vector& b);
#endif