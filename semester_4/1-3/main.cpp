#include "functions.hpp"

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;


int main() {
    Matrix A;
    Vector x, b;
    int n = 4;

    A = generateRandomMatrix(x, n);
    b = multiplyMatrixVector(A, x);
    printMatrix(A, x, b);

    for (int p = -1; p > -2; p--)
    {
        double tolerance = pow(10, p);
        int numberOfIterations = conjugateGradientMethod(A, b, tolerance);
        cout << tolerance << " | " << numberOfIterations<<"\n\n"; 
    }

   return 0;
}