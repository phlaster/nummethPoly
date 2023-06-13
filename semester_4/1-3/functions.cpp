#include "Types.hpp"
#include <algorithm>
#include <stdexcept>
#include <vector>
#include <random>
#include <iostream>
#include <iomanip>
#include <cmath>


using namespace std;

Matrix generateRandomMatrix(Vector& x, int n) {
    // Инициализация генератора случайных чисел
    random_device rd;
    mt19937 gen(rd());
    
    // Создание распределения равномерного распределения на интервале [0, 10]
    uniform_real_distribution<> distrib(0.0, 10.0);
    
    Matrix A;
    A.resize(n);
    for (int i = 0; i < n; ++i) {
        A[i].resize(n);
        x.push_back(distrib(gen)); // Модификация вектора решений
    }
     
    // Заполнение матрицы случайными значениями из заданного распределения 
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A[i][j] = distrib(gen);
        }
    }

    return A;
}

Vector generateRandomVector(int n) {
    random_device rd;
    mt19937 gen(rd());
    
    // Создание распределения равномерного распределения на интервале [0, 10]
    uniform_real_distribution<> distrib(0, 10.0);
    
    Vector x0;

    for (int i = 0; i < n; ++i) {
        x0.push_back(distrib(gen));
    }

    return x0;
}

void printMatrix(const Matrix& A, const Vector& x,
                 const Vector& b) {
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

Matrix matmul(const Matrix& A, const Matrix& B) {
    int n = A.size();
    int m = B.size();
    int k = B[0].size();
    
    // Проверка возможности умножения матриц
    if (A[0].size() != m) {
        throw invalid_argument("Invalid matrix sizes");
    }
    
     // Создание новой матрицы с нулевыми значениями
     Matrix res(n, Vector(k));
     
     for (int i = 0; i < n; ++i) {
         for (int j = 0; j < k; ++j) {
             double sum = 0.0;
             
             for (int l = 0; l < m; ++l) {                
                 sum += A[i][l] * B[l][j];
             }
             
             res[i][j] = sum;
         }         
     }
     
     return res;
}

Vector scalarMultiply(double scalar, const Vector& vec) {
    Vector result;
    for (double value : vec) {
        result.push_back(scalar * value);
    }
    return result;
}

Vector vecSum(double c1, const Vector& vec1, double c2, const Vector& vec2) {
   if (vec1.size() != vec2.size()) {
      throw std::invalid_argument("Vectors must be same length.");
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

double dotProduct(const Vector& vec1, const Vector& vec2){
    if (vec1.size() != vec2.size()) {
        throw invalid_argument("Dot product is only defined for two vectors of same length.");
    }
    double result = 0.0;
    for (int i = 0; i < vec1.size(); ++i) {
        result += vec1[i] * vec2[i];
    }
    return result;
}

double euclideanNorm(const Vector& vec) {
    return sqrt(dotProduct(vec, vec));
}

Vector multiplyMatrixVector(const Matrix& A, const Vector& x) {
    if (A[0].size() != x.size()) {
        throw std::invalid_argument("Matrix and vector must be compatible.");
    }
    Vector result;
    for (const auto& row : A) {       
        result.push_back(dotProduct(row, x));
    }
    return result;
}

bool goodSolution(const Vector& r, double tolerance) {       
    return euclideanNorm(r) <= tolerance;
}

// Вычисление вектора невязки
Vector computeResiduals(const Matrix& A, const Vector& b, const Vector& x_i)
{
    Vector r = vecSum(1, b, -1, multiplyMatrixVector(A, x_i));
    return r;
}

// Вычисления направления поиска p на каждой итерации
Vector computeDirection(const Matrix& A, const Vector& p) {
    Vector Ap = multiplyMatrixVector(A, p);
    double alpha = dotProduct(p, p) / dotProduct(p, Ap);
    return scalarMultiply(alpha, p);
}

// Проверка деления на ноль при расчете коэффициента alpha_k или beta
void checkDivideByZero(double value, double tolerance)
{
    if(fabs(value)<tolerance){
        throw runtime_error("Error: division by zero");
    }
}

int conjugateGradientMethod(const Matrix& A, const Vector& b, double tolerance)
{
    int n = A.size();
    Vector x_i = generateRandomVector(n);
    Vector r = computeResiduals(A, b, x_i);
    Vector p = r;
    int k = 0;

    while (!goodSolution(r, tolerance) && k < 10) {
        Vector q = multiplyMatrixVector(A, p);
        double alpha_num = dotProduct(r, p);
        double alpha_denom = dotProduct(q, p);
        
        checkDivideByZero(alpha_denom, tolerance);
        double alpha = alpha_num/alpha_denom; 
        // cout << "alpha: " << alpha << "\n";
        x_i = vecSum(1, x_i, alpha, p);

        Vector r_next = vecSum(1, r, -alpha, q);

        if(goodSolution(r_next, tolerance)) break;

        double beta_num = dotProduct(r_next,q);        
        double beta_denom = dotProduct(p,q);  
        
        checkDivideByZero(beta_denom, tolerance);
        
        double beta=(beta_num/beta_denom);

        cout <<"k=" << k <<
            ", |r|=" << euclideanNorm(r) <<
            ", a=" << alpha <<
            ", b=" << beta << "\n";

        for (auto c : x_i)
        {
            cout << c << " ";
        }
        cout << endl;
        
        // cout << "beta: " << beta << "\n";
        p = vecSum(1, r_next, -beta, p);
        swap(r,r_next);
        k++;
    }
    return k;
}


