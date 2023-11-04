#include "headers/operators.hpp"


Vec operator+(Vec a, const Vec& b) {
    size_t n = a.size();
    if (n != b.size()) {
        throw invalid_argument("Vector sizes are not aligned for addition.");
    }
    
    for (size_t i = 0; i < n; ++i) {
        a[i] += b[i];
    }
    return a;
}

Vec operator-(Vec a){
    size_t n = size(a);
    for (size_t i=0; i<n; i++){
        a[i] *= -1.0;
    }
    return a;
}

Vec operator-(Vec a, const Vec& b){
    size_t n = a.size();
    if (n != b.size()) {
        throw invalid_argument("Vector sizes are not aligned for subtraction.");
    }
    for (size_t i=0; i<n; i++){
        a[i] -= b[i];
    }
    return a;
}

Vec operator*(const double c, Vec a){
    size_t n = a.size();
    
    for (size_t i = 0; i < n; ++i) {
        a[i] *= c;
    }
    return a;
}

// Dot product
double operator*(const Vec& a, const Vec& b){
    if (a.size() != b.size())
        throw invalid_argument("Vector sizes are not aligned for dot product.");

    double res = 0.0;
    for (size_t i = 0; i < a.size(); ++i)
        res += a[i] * b[i];
    return res;
}

Vec operator/(Vec a, const double c){
    if (fabs(c)<1e-15) {
        throw invalid_argument("Can not devide vector by zero.");
    }
    size_t n = size(a);
    for (size_t i=0; i<n; i++){
        a[i] /= c;
    }
    return a;
}



Mtr operator+(Mtr A, const Mtr& B) {
    size_t n = A.size(), m = A[0].size();
    if (n != B.size() || m != B[0].size()) {
        throw invalid_argument("Matrix dimensions are not aligned for addition.");
    }

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            A[i][j] += B[i][j];
        }
    }
    return A;
}

Mtr operator-(Mtr A) {
    size_t n = A.size(), m = A[0].size();
    
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            A[i][j] *= -1.0;
        }
    }
    return A;
}

Mtr operator-(Mtr A, const Mtr& B) {
    size_t n = A.size(), m = A[0].size();
    if (n != B.size() || m != B[0].size()) {
        throw invalid_argument("Matrix dimensions are not aligned for subtraction.");
    }

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            A[i][j] -= B[i][j];
        }
    }
    return A;
}

Mtr operator*(const double c, Mtr A) {
    size_t n = A.size(), m = A[0].size();
    
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            A[i][j] *= c;
        }
    }
    return A;
}

Vec operator*(const Mtr& A, const Vec& x) {
    if (A[0].size() != x.size())
        throw invalid_argument("Matrix and vector dimentions are not aligned for multiplication.");
    size_t res_len = A.size();

    Vec res(res_len);
    for (size_t i=0; i<res_len; i++){
        res[i] = A[i] * x;
    }
    return res;
}
// Mtr mul
Mtr operator*(const Mtr& A, const Mtr& B) {
    size_t n = A.size();
    size_t m = B.size();
    size_t k = B[0].size();

    if (A[0].size() != m)
        throw invalid_argument("Matrix dimentions are not aligned for matrix multiplication.");

    Mtr res(n, Vec(k));

    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < k; ++j) {
            double S = 0.0;
            for (size_t l = 0; l < m; ++l)
                S += A[i][l] * B[l][j];
            res[i][j] = S;
        }
    return res;
}

Mtr operator/(Mtr A, const double c){
    if (fabs(c)<1e-15) {
        throw invalid_argument("Can not devide matrix by zero.");
    }
    size_t n = A.size(), m = A[0].size();
    
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            A[i][j] /= c;
        }
    }
    return A;
}




spMtr operator+(spMtr A, const spMtr& B) {
    size_t n = A.rows, m = A.cols;
    if (n != B.rows || m != B.cols) {
        throw invalid_argument("Matrix dimensions are not aligned for addition.");
    }
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            double res = A.get(i,j) + B.get(i,j);
            A.set(res, i,j);
        }
    }
    return A;
}

spMtr operator-(spMtr A) {
    size_t n = A.rows, m = A.cols;
    
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            double val = -1.0 * A.get(i,j);
            A.set(val,i,j);
        }
    }
    return A;
}

spMtr operator-(const spMtr& A, const spMtr& B) {
    size_t n = A.rows, m = A.cols;
    if (n != B.rows || m != B.cols) {
        throw invalid_argument("Matrix dimensions are not aligned for subtraction.");
    }
    return A + (-B);
}

spMtr operator*(const double c, spMtr A) {
    size_t n = A.rows, m = A.cols;
    
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            double val = A.get(i,j) * c;
            A.set(val, i, j);
        }
    }
    return A;
}

Vec operator*(const spMtr& A, const Vec& x) {
    if (A.cols != x.size())
        throw invalid_argument("Matrix and vector dimentions are not aligned for multiplication.");

    Vec res(A.rows);
    for (size_t i=0; i<A.rows; i++){
        res[i] = A.get(i) * x;
    }
    return res;
}
// Mtr mul
spMtr operator*(const spMtr& A, const spMtr& B) {
    size_t n = A.rows;
    size_t m = B.rows;
    size_t k = B.cols;

    if (A.cols != m)
        throw invalid_argument("Matrix dimentions are not aligned for matrix multiplication.");

    spMtr res(n, k);

    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < k; ++j) {
            double S = 0.0;
            for (size_t l = 0; l < m; ++l){
                S += A.get(i,l) *B.get(l, j);
            }
            res.set(S, i, j);
        }
    return res;
}

spMtr operator/(spMtr A, const double c){
    if (fabs(c)<1e-15) {
        throw invalid_argument("Can not devide matrix by zero.");
    }
    size_t n = A.rows, m = A.cols;
    
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            double val = A.get(i,j) / c;
            A.set(val, i, j);
        }
    }
    return A;
}


void operator+=(spMtr& A, const spMtr& B) {
    size_t n = A.rows, m = B.cols;
    if (n != B.rows || m != B.cols) {
        throw std::invalid_argument("Matrix dimensions are not aligned.");
    }
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            double val = A.get(i,j) + B.get(i,j);
            A.set(val, i, j);
        }
    }
}