#include <random>
#include <vector>
#include <stdexcept>
#include <cmath>

#include <iostream>
#include <iomanip>

#include "functions.hpp"


class Matrix {
public:
    Matrix(int rows, int cols) : data(rows, std::vector<double>(cols)) {}

    void fill(double low, double high) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dist(low, high);

        for (auto& row : data) {
            for (auto& element : row) {
                element = dist(gen);
            }
        }
    }

    void print(int d) const {
        for (const auto& row : data) {
            for (const auto& element : row) {
                std::cout << std::fixed << std::setprecision(d) << element << " ";
            }
            std::cout << std::endl;
        }
    }

    Matrix matrixMultiply(const Matrix& matrix2) const {
        // Check if dimensions are correct
        if (data[0].size() != matrix2.data.size()) {
            throw std::invalid_argument("Matrix dimensions are not compatible for multiplication");
        }

        // Create result matrix
        Matrix result(data.size(), std::vector<double>(matrix2.data[0].size()));

        // Perform multiplication
        for (int i = 0; i < data.size(); i++) {
            for (int j = 0; j < matrix2.data[0].size(); j++) {
                double sum = 0;
                for (int k = 0; k < matrix2.data.size(); k++) {
                    sum += data[i][k] * matrix2.data[k][j];
                }
                result.data[i][j] = sum;
            }
        }

        return result;
    }


private:
    std::vector<std::vector<double>> data;
};




// Vector vectorMultiply(const Vector& vec, double scalar) {
//     Vector result(vec.size());
//     for (int i = 0; i < vec.size(); i++) {
//         result[i] = vec[i] * scalar;
//     }
//     return result;
// }

// Vector vectorAdd(const Vector& vec1, const Vector& vec2) {
//     if (vec1.size() != vec2.size()) {
//         throw std::invalid_argument("Vectors are not the same length");
//     }
//     Vector result(vec1.size());
//     for (int i = 0; i < vec1.size(); i++) {
//         result[i] = vec1[i] + vec2[i];
//     }
//     return result;
// }

// double dotProduct(const Vector& vec1, const Vector& vec2) {
//     if (vec1.size() != vec2.size()) {
//         throw std::invalid_argument("Vectors are not the same length");
//     }
//     double result = 0;
//     for (int i = 0; i < vec1.size(); i++) {
//         result += vec1[i] * vec2[i];
//     }
//     return result;
// }

// double euclideanNorm(const Vector& vec) {
//     double result = 0;
//     for (int i = 0; i < vec.size(); i++) {
//         result += pow(vec[i], 2);
//     }
//     return sqrt(result);
// }