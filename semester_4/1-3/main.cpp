#include "functions.hpp"
#include "Types.hpp"
#include <iostream>
#include <cmath>

using namespace std;

void PCG(const Matrix& A, const Vector& b, int lowestDeg=-15, int maxIter=1000, bool verbose=false)
{
    cout << "eps,iters,n="<<b.size()<<"\n";
    for (int p = -1; p >= lowestDeg; p--)
    {
        double eps = pow(10, p);
        int n_iters = conjugateGradientMethod(A, b, eps, maxIter, verbose);
        if (n_iters == -1)
        {
            cout << eps << " и далее " <<">"<< maxIter <<" итераций.\n";
            break;
        } else {
            cout << eps << "," << n_iters<<"\n"; 
        }
    }
}

int main() {
    Matrix A;
    Vector x, b;
    int n = 100; // Размер матрицы
    A = generateRndSymPos(n);
    x = generateRandomVector(n);
    b = multiplyMatrixVector(A, x);
    PCG(A, b);

    /* Можно запустить матрицу из примера в отчёте:
    A = {
        { 5,-1,-1, 0},
        {-1, 6, 1,-1},
        {-1, 1, 6, 1},
        { 0,-1, 1, 7}
    };
    x = {-2, 0, 2, 1};
    b = {-12, 3, 15, 9};
    cout << "x = "; printVector(x);
    PCG(A, b, -1, 5, true);
    */
    

    return 0;
}