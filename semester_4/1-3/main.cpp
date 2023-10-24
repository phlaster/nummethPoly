#include "Functions.hpp"


int main() {
    int size = 100; // Размер матрицы
    double cond = 3;
    const int MAXITERS = 1000;

    Mtr A = generateRndSymPos(size, 1);
    Vec x = generateRandomVector(size);
    Vec b = multiplyMatrixVector(A, x);

    print(A, x, b);
    PCG(A, b, -15, MAXITERS);
    return 0;
}
