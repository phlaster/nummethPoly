#include "Functions.hpp"

Mtr upperTriang(const Mtr& A){
    if (A.size() != A[0].size())
        throw invalid_argument("Только квадратные матрцы!");

    int n = A.size();
    Mtr UT(n, Vec(n));

    for (int i = 0; i < n; ++i)
        for (int j = i; j < n; ++j) {
            UT[i][j] = A[i][j];
            UT[j][i] = A[i][j];
        }
    return UT;
}

Mtr transpose(const Mtr& matrix){
    int rows = matrix.size();
    int cols = matrix[0].size();

    Mtr result(cols, Vec(rows));

    for (int i = 0; i < rows; ++i)
        for(int j = 0; j < cols; ++j)
            result[j][i] = matrix[i][j];
    return result;
}

Vec scalarMultiply(double scalar, const Vec& vec){
    Vec result;
    for (double value : vec)
        result.push_back(scalar * value);
    return result;
}

Vec vecSum(double c1, const Vec& vec1, double c2, const Vec& vec2){
    if (vec1.size() != vec2.size())
        throw std::invalid_argument("Можно складывать только векторы одинаковой длины!");
   
    Vec v1 = scalarMultiply(c1, vec1), v2 = scalarMultiply(c2, vec2);

    int len = v1.size();
    Vec result(len);

    for (int i = 0; i < len; ++i)
        result[i] = v1[i] + v2[i];
    return result;
}

double dotProduct(const Vec& vec1, const Vec& vec2){
    if (vec1.size() != vec2.size())
        throw invalid_argument("Скалярное произведение только для векторов одинаковой длины!");
    double result = 0.0;
    for (int i = 0; i < vec1.size(); ++i)
        result += vec1[i] * vec2[i];
    return result;
}

double euclideanNorm(const Vec& vec){
    return sqrt(dotProduct(vec, vec));
}

Vec multiplyMatrixVector(const Mtr& A, const Vec& x){
    if (A[0].size() != x.size())
        throw std::invalid_argument("Для умножения матрицы на вектор размеры должны согласовываться!");
    Vec result;
    for (const auto& row : A)   
        result.push_back(dotProduct(row, x));
    return result;
}

Mtr matMul(const Mtr& A, const Mtr& B) {
    if (A.empty() || B.empty() || A[0].size() != B.size())
        throw invalid_argument("Размеры матриц несовместимы для умножения!");

    int n = A.size();
    int m = B[0].size();
    int p = B.size();

    Mtr res(n, vector<double>(m, 0));

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            for (int k = 0; k < p; ++k)
                res[i][j] += A[i][k] * B[k][j];

    return res;
}


Vec generateRandomVector(int n, double lower, double upper) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> distrib(lower, upper);

    Vec vec;
    for (int i = 0; i < n; ++i) {
        double r = distrib(gen);
        vec.push_back(r);
    }
    return vec;
}

Mtr generateRandomMatrix(int rows, int cols, double lower, double upper){
    cols = cols <= 0 ? rows : cols;
    Mtr A;
    for (int i = 0; i < rows; ++i)
        A.push_back(generateRandomVector(cols, lower, upper));
    return A;
}

double mean(vInt V){
    double S = 0.0;
    for (auto x : V)
        S += x;
    return S / V.size();
}

Mtr E(size_t size){
    Mtr m = Mtr(size, Vec(size, 0.0));
    for (size_t i=0; i<size; i++)
        m[i][i] = 1.0;
    return m;
}

Vec mul(double scalar, const Vec& V){
    Vec result;
    for (double value : V)
        result.push_back(scalar * value);
    return result;
}

Mtr mul(double scalar, const Mtr& mtr){
    Mtr result;
    for (Vec row : mtr)
        result.push_back(mul(scalar, row));
    return result;
}

Vec normalize(const Vec& V){
    double norm = euclideanNorm(V);
    if (norm==0.0)
        throw invalid_argument("Нельзя нормировать нулевой вектор!");
    return mul(1/norm, V);
}


Mtr sum(const Mtr& m1, const Mtr& m2, double c1, double c2) {
    if (m1.size() != m2.size() ||
        m1.empty() ||
        m2.empty() ||
        m1[0].size() != m2[0].size())
        throw std::invalid_argument("Можно складывать только матрицы одинакового размера!");

    Mtr mat1 = mul(c1, m1);
    Mtr mat2 = mul(c2, m2);

    size_t rows = m1.size();
    size_t cols = m1[0].size();
    Mtr result(rows, Vec(cols));

    for (size_t i = 0; i < rows; ++i)
        for (size_t j = 0; j < cols; ++j)
            result[i][j] = mat1[i][j] + mat2[i][j];
    return result;
}


Mtr householder(const Vec& w0){
    
    auto fill = [](double val, size_t rows, size_t cols) -> Mtr{
        return Mtr(rows, Vec(cols, val));
    };

    auto toCol = [&](const Vec& v) -> Mtr{
        Mtr alloc = fill(0, size(v), 1);
        for (size_t i=0; i< size(v); i++)
            alloc[i][0] = v[i];
        return alloc;
    };

    auto toRow = [](const Vec& v) -> Mtr{
        return Mtr(1, v);
    };

    size_t n = size(w0);
    Vec w = normalize(w0);
    Mtr Q = sum(E(n), matMul(toCol(w), toRow(w)), 1, -2);
    return Q;
}

Mtr generateRndSymPos(int n, double cond){
    Vec w0 = generateRandomVector(n);
    Mtr H = householder(w0);
    Mtr D = E(n);

    Vec diag_elems = generateRandomVector(n-1, 1, cond);
    for (int i=0; i<n-1; i++)
        D[i][i] = diag_elems[i];
    D[0][0] = cond;

    Mtr A = matMul(matMul(transpose(H), D), H);
    return A;
}

Mtr generateRndSymPos(int n){
    Mtr M = generateRandomMatrix(n);
    Mtr A = matMul(M, transpose(M));
    return A;
}