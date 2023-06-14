#include "Types.hpp"
#include <stdexcept>
#include <vector>
#include <random>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

Matrix upperTriang(const Matrix& A)
{
    if (A.size() != A[0].size()) {
        throw invalid_argument("Только квадратные матрцы!");
    }

    int n = A.size();
    Matrix UT(n, Vector(n));

    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            UT[i][j] = A[i][j];
            UT[j][i] = A[i][j];
        }
    }
    return UT;
}
Matrix transpose(const Matrix& matrix)
{
    int rows = matrix.size();
    int cols = matrix[0].size();

    Matrix result(cols, Vector(rows));

    for (int i = 0; i < rows; ++i) {
        for(int j = 0; j < cols; ++j) {
            result[j][i] = matrix[i][j];
        }
    }

    return result;
}


Vector scalarMultiply(double scalar, const Vector& vec)
{
    Vector result;
    for (double value : vec) {
        result.push_back(scalar * value);
    }
    return result;
}
Vector vecSum(double c1, const Vector& vec1, double c2, const Vector& vec2)
{
   if (vec1.size() != vec2.size()) {
      throw std::invalid_argument("Можно складывать только векторы одинаковой длины!");
   }
   
   Vector v1 = scalarMultiply(c1, vec1);
   Vector v2 = scalarMultiply(c2, vec2);

   int len = v1.size();
   Vector result(len);

   for (int i = 0; i < len; ++i) {
      result[i] = v1[i] + v2[i];
   }
   return result;
}


double dotProduct(const Vector& vec1, const Vector& vec2)
{
    if (vec1.size() != vec2.size()) {
        throw invalid_argument("Скалярное произведение только для векторов одинаковой длины!");
    }
    double result = 0.0;
    for (int i = 0; i < vec1.size(); ++i) {
        result += vec1[i] * vec2[i];
    }
    return result;
}
double euclideanNorm(const Vector& vec)
{
    return sqrt(dotProduct(vec, vec));
}


Vector multiplyMatrixVector(const Matrix& A, const Vector& x)
{
    if (A[0].size() != x.size()) {
        throw std::invalid_argument("Для умножения матрицы на вектор размеры должны согласовываться!");
    }
    Vector result;
    for (const auto& row : A) {       
        result.push_back(dotProduct(row, x));
    }
    return result;
}
Matrix matMul(const Matrix& A, const Matrix& B)
{
    int n = A.size();
    int m = B.size();
    int k = B[0].size();
    
    // Проверка возможности умножения матриц
    if (A[0].size() != m) {
        throw invalid_argument("Размеры матриц несовместимы для умножения!");
    }
    Matrix At = transpose(A);

    // Создание новой матрицы с нулевыми значениями
    Matrix res(n, Vector(k));
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < k; ++j) {
            res[i][j] = dotProduct(At[i], B[j]);
        }         
    }
     
     return res;
}


Vector generateRandomVector(int n, double lower=0, double upper=10)
{
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> distrib(lower, upper);
    
    Vector vec;
    for (int i = 0; i < n; ++i) {
        vec.push_back(distrib(gen));
    }

    return vec;
}
Matrix generateRandomMatrix(int rows, int cols=-1, double lower=0, double upper=10)
{
    cols = cols <= 0 ? rows : cols;
    Matrix A;
    for (int i = 0; i < rows; ++i) {
        A.push_back(generateRandomVector(cols, lower, upper));
    }
    return A;
}
Matrix generateRndSymPos(int n)
{
    Matrix M = generateRandomMatrix(n);
    return upperTriang(matMul(M, transpose(M)));
}


int conjugateGradientMethod(const Matrix& A, const Vector& b, double eps, int maxIter=1000)
{
    int n = A.size();
    int k = 0;
    Vector q;
    Vector x = Vector(n);
    Vector r = vecSum(1, b, -1, multiplyMatrixVector(A, x));
    Vector p = r;
    double alpha, beta, pq_denom;
    while (euclideanNorm(r) > eps && k <= maxIter)
    {
        q = multiplyMatrixVector(A, p);
        pq_denom = dotProduct(p, q);
        alpha = dotProduct(r, p) / pq_denom; 
        r = vecSum(1, r, -alpha, q);        
        x = vecSum(1, x, alpha, p);
        beta = dotProduct(r, q) / pq_denom;
        p = vecSum(1, r, -beta, p);
        k++;
    }
    return k <= maxIter ? k : -1;
}


void printMatrix(const Matrix& A)
{
    int n = A.size();
	for(int i=0; i<n; ++i){
		for(int j=0; j<n; ++j)
			cout << setw(4) << A[i][j] << " ";
		cout << "\n";
	}
    cout << "\n";
}
void printSLAE(const Matrix& A, const Vector& x, const Vector& b)
{
    int n = A.size();
    
    // Выводим заголовки столбцов
    cout << "A:" << setw(9*n - 1) << "| x:" << setw(10) << "| b:" << '\n';

	// Выводим значения матрицы, вектора x и вектора b
	for(int i=0; i<n; ++i){
		for(int j=0; j<n; ++j)
			cout << setw(4) << A[i][j] << " ";
		cout<<"| "<< setw(4) << x[i] << " |"<< setw(4) << b[i] << "\n";
	}
    cout << "\n";
}
