#include "functions.hpp"
#include "Types.hpp"
#include <iostream>
#include <cmath>

using namespace std;

void PCG(const Matrix& A, const Vector& b, int lowestDeg, int maxIter = 1000)
{
    cout << "eps,iters,n="<<b.size()<<"\n";
    for (int p = -1; p >= lowestDeg; p--)
    {
        double eps = pow(10, p);
        int n_iters = conjugateGradientMethod(A, b, eps, maxIter);
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
    int n = 100;
    A = generateRndSymPos(n);
    x = generateRandomVector(n);
    b = multiplyMatrixVector(A, x);
    PCG(A, b, -16, 500);

    return 0;
}