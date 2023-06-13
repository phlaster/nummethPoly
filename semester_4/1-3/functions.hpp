#include "Types.hpp"

#ifndef GRM_HPP
#define GRM_HPP
Matrix generateRandomMatrix(Vector& x, int n);
#endif

#ifndef GRV_HPP
#define GRV_HPP
Vector generateRandomVector(int n);
#endif

#ifndef SM__HPP
#define SM__HPP
Vector scalarMultiply(double scalar, const Vector& vec);
#endif

#ifndef VS__HPP
#define VS__HPP
Vector vecSum(double c1, const Vector& vec1, double c2, const Vector& vec2);
#endif

#ifndef DP__HPP
#define DP__HPP
double dotProduct(const Vector& vec1, const Vector& vec2);
#endif

#ifndef EN__HPP
#define EN__HPP
double euclideanNorm(const Vector& vec);
#endif

#ifndef MM_HPP
#define MM_HPP
Matrix matmul(const Matrix& A, const Matrix& B);
#endif

#ifndef MMV_HPP
#define MMV_HPP
Vector multiplyMatrixVector(const Matrix& A, const Vector& x);
#endif

#ifndef CSN_HPP
#define CSN_HPP
bool checkSolutionNorm(const Matrix& A, const Vector& b, const Vector& x, double tolerance);
#endif

#ifndef CGM_HPP
#define CGM_HPP
int conjugateGradientMethod(const Matrix& A, const Vector& b, double tolerance);
#endif

#ifndef PMX_HPP
#define PMX_HPP
void printMatrix(const Matrix& A, const Vector& x, const Vector& b);
#endif
