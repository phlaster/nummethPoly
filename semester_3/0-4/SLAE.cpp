#include "LinearAlgebra.hpp"

Vec solveLinearEquation(const Mtr& L,
                        const Mtr& U,
                        const vInt& permutation,
                        const Vec& b) {
    int n = L.size();
    
    // Прямая подстановка (Ly = Pb)
    Vec y(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j <= i - 1; ++j)
            sum += L[i][j] * y[j];
        y[i] = b[permutation[i]] - sum;
    }

    // Обратная подстановка (Ux = y)
    Vec x(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j= i + 1 ; j < n ; ++j)
            sum += U[i][j] * x[j];
        x[i] =(y[i]-sum)/U[i][i];
    }
    return x;
}

Vec residual(const Mtr& A, const Vec& x, const Vec& b) {
    int n = A.size();
    Vec residual(n, 0.0);
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            residual[i] += A[i][j] * x[j];
        residual[i] -= b[i];
    }
    return residual;
}