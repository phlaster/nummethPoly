#include <random>

void fillMatrix(std::vector<std::vector<double>>& matrix) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    for (auto& row : matrix) {
        for (auto& element : row) {
            element = dist(gen);
        }
    }
}


#include <vector>
#include <stdexcept>

std::vector<std::vector<double>> matrixMultiply(const std::vector<std::vector<double>>& matrix1, const std::vector<std::vector<double>>& matrix2) {
    // Check if dimensions are correct
    if (matrix1[0].size() != matrix2.size()) {
        throw std::invalid_argument("Matrix dimensions are not compatible for multiplication");
    }

    // Create result matrix
    std::vector<std::vector<double>> result(matrix1.size(), std::vector<double>(matrix2[0].size()));

    // Perform multiplication
    for (int i = 0; i < matrix1.size(); i++) {
        for (int j = 0; j < matrix2[0].size(); j++) {
            double sum = 0;
            for (int k = 0; k < matrix2.size(); k++) {
                sum += matrix1[i][k] * matrix2[k][j];
            }
            result[i][j] = sum;
        }
    }

    return result;
}

#include <vector>
#include <stdexcept>

std::vector<double> vectorMultiply(const std::vector<double>& vec, double scalar) {
    std::vector<double> result(vec.size());
    for (int i = 0; i < vec.size(); i++) {
        result[i] = vec[i] * scalar;
    }
    return result;
}

std::vector<double> vectorAdd(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    if (vec1.size() != vec2.size()) {
        throw std::invalid_argument("Vectors are not the same length");
    }
    std::vector<double> result(vec1.size());
    for (int i = 0; i < vec1.size(); i++) {
        result[i] = vec1[i] + vec2[i];
    }
    return result;
}

double dotProduct(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    if (vec1.size() != vec2.size()) {
        throw std::invalid_argument("Vectors are not the same length");
    }
    double result = 0;
    for (int i = 0; i < vec1.size(); i++) {
        result += vec1[i] * vec2[i];
    }
    return result;
}

#include <cmath>

double euclideanNorm(const std::vector<double>& vec) {
    double result = 0;
    for (int i = 0; i < vec.size(); i++) {
        result += pow(vec[i], 2);
    }
    return sqrt(result);
}