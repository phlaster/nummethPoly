#include <vector>

#ifndef MATRIX_HPP
#define MATRIX_HPP
class Matrix {
private:
    int rows;
    int cols;
    int **data;
public:
    Matrix(int rows, int cols);
    ~Matrix();
    void fill();
    void print();
    Matrix matrixMultiply(Matrix other);
};
#endif

// #ifndef VECTOR_HPP
// #define VECTOR_HPP
// typedef std::vector<double> Vector;
// #endif

// #ifndef MATRIX_HPP
// #define MATRIX_HPP
// typedef std::vector<Vector> Matrix;
// #endif


// #ifndef FILL_MATRIX_HPP
// #define FILL_MATRIX_HPP
// void fillMatrix(Matrix & matrix);
// #endif

// #ifndef PRINT_MATRIX_HPP
// #define PRINT_MATRIX_HPP
// void printMatrix(const Matrix& matrix, int d);
// #endif


