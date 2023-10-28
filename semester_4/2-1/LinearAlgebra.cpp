#include "LinearAlgebra.hpp"
#include <cmath>
#include <math.h>
#include <stdexcept>
#include <tuple>
#include <utility>

bool issquare(const Mtr& M){
    return M.size() == M[0].size();
}
bool isapprox(const Vec& v1, const Vec& v2, double tol) {
    if (v1.size() != v2.size())
        return false;
    for (size_t i = 0; i < v1.size(); ++i)
        if (fabs(v1[i] - v2[i]) > tol)
            return false;
    return true;
}
bool isapprox(const Mtr& m1, const Mtr& m2, double tol) {
    if (m1.size() != m2.size() || m1[0].size() != m2[0].size())
        return false;
    for (size_t i = 0; i < m1.size(); ++i)
        if (!isapprox(m1[i], m2[i], tol))
            return false;
    return true; 
}
bool issingleval(const Mtr& M){
    return issquare(M) && M.size() == 1;
}


Mtr upperTrSymmetric(const Mtr& A){
    if (!issquare(A))
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


Mtr T(const Mtr& matrix){
    int rows = matrix.size();
    int cols = matrix[0].size();

    Mtr result(cols, Vec(rows));

    for (int i = 0; i < rows; ++i)
        for(int j = 0; j < cols; ++j)
            result[j][i] = matrix[i][j];
    return result;
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

Vec sum(const Vec& v1, const Vec& v2, double c1, double c2){
    if (v1.size() != v2.size())
        throw invalid_argument("Можно складывать только векторы одинаковой длины!");
   
    Vec V1 = mul(c1, v1);
    Vec V2 = mul(c2, v2);

    int len = v1.size();
    Vec result(len);

    for (int i = 0; i < len; ++i)
        result[i] = V1[i] + V2[i];
    return result;
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
double sum(const Vec& v){
    double S = 0;
    for (auto x : v) S+= x;
    return S;
}

double dot(const Vec& V1, const Vec& V2){
    if (V1.size() != V2.size())
        throw invalid_argument("Скалярное произведение только для векторов одинаковой длины!");

    double result = 0.0;
    for (size_t i = 0; i < V1.size(); ++i)
        result += V1[i] * V2[i];
    return result;
}

double euclideanNorm(const Vec& V){
    return sqrt(dot(V, V));
}

// double euclideanNorm(const Mtr& M) {
//     double sum_squares = 0.0;
    
//     for (const auto& row : M)
//         for (const auto& elem : row)
//             sum_squares += elem * elem;
//     return sqrt(sum_squares);
// }

double opnorm_1(const Mtr& M) {
    size_t numRows = M.size();
    size_t numCols = M[0].size();
    double maxSum = -std::numeric_limits<double>::infinity();

    for (size_t j = 0; j < numCols; j++) {
        double sum = 0.0;
        for (size_t i = 0; i < numRows; i++) {
            sum += std::abs(M[i][j]);
        }

        if (sum > maxSum) {
            maxSum = sum;
        }
    }
    return maxSum;
}


Vec normalize(const Vec& V){
    double norm = euclideanNorm(V);
    if (norm==0.0)
        throw invalid_argument("Нельзя нормировать нулевой вектор!");
    return mul(1/norm, V);
}

Vec mul(const Mtr& A, const Vec& x){
    if (A[0].size() != x.size())
        throw invalid_argument("Для умножения матрицы на вектор размеры должны согласовываться!");

    Vec result;
    for (const auto& row : A)  
        result.push_back(dot(row, x));
    return result;
}


Mtr mul(const Mtr& A, const Mtr& B) {
    if (issingleval(A)) return mul(A[0][0], B);
    if (issingleval(B)) return mul(B[0][0], A);

    size_t n = A.size();
    size_t m = B.size();
    size_t k = B[0].size();

    if (A[0].size() != m)
        throw invalid_argument("Размеры матриц несовместимы для умножения!");

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


Vec randVec(int n, double lower, double upper){
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> distrib(lower, upper);
    
    Vec V;
    for (int i = 0; i < n; ++i)
        V.push_back(distrib(gen));
    return V;
}

double rand_double(){
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> distrib(0, 1);
    return distrib(gen);
}

Mtr randMtr(int rows, int cols, double lower, double upper){
    cols = cols <= 0 ? rows : cols;
    Mtr A;
    for (int i = 0; i < rows; ++i)
        A.push_back(randVec(cols, lower, upper));

    return A;
}

// Mtr randSymmPositive(int n){
//     Mtr M = randMtr(n);
//     return upperTrSymmetric(mul(M, T(M)));
// }

Vec fill(double val, size_t length){
    return Vec(length, val);
}

Mtr fill(double val, size_t rows, size_t cols){
    return Mtr(rows, Vec(cols, val));
}

Mtr E(size_t size){
    Mtr m = fill(0.0, size, size);
    for (size_t i=0; i<size; i++)
        m[i][i] = 1.0;
    return m;
}

Mtr _getSubmatrix(const Mtr& matrix, size_t i, size_t j) {
    Mtr submatrix;
    for (size_t row = 0; row < matrix.size(); ++row)
        if (row != i) {
            Vec temp_row;
            for (size_t col = 0; col < matrix[row].size(); ++col)
                if (col != j)
                    temp_row.push_back(matrix[row][col]);
            submatrix.push_back(temp_row);
        }
    return submatrix;
}

double det(const Mtr& M) {
    if (!issquare(M))
        throw invalid_argument("Определитель считают только для квадратных матриц!");

    size_t n = M.size();
    if (n == 1) return M[0][0];
  
  	double determinant = 0.0;
  	int sign = 1;

    for(size_t i=0; i<n; ++i){
        determinant += sign * M[0][i] * det(_getSubmatrix(M, 0, i));
        sign *= -1;
    }

    return determinant;
}

// double cond(const Mtr& M){
//     return euclideanNorm(M) * euclideanNorm(inv(M));
// }

Mtr inv(const Mtr& A) {
    if (!issquare(A))
        throw invalid_argument("Обратная матрица существует только для квадратных!");

    double determinant = det(A);
    if (determinant == 0)
        throw invalid_argument("Обратной матрицы нет, det=0!");


    int n = A.size();
    Mtr inv = fill(0.0, n, n);

    if (n == 1){
        inv[0][0] = 1/A[0][0];
        return inv;
    }

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            Mtr submatrix = _getSubmatrix(A, i, j);
            double submatrixDeterminant = det(submatrix);
            double cofactor = pow(-1, i + j) * submatrixDeterminant;
            inv[j][i] = cofactor / determinant;
        }
    return inv;
}

Mtr toCol(const Vec& v){
    Mtr alloc = fill(0, size(v), 1);
    for (size_t i=0; i< size(v); i++)
        alloc[i][0] = v[i];
    return alloc;
}

Mtr toRow(const Vec& v){
    return Mtr(1, v);
}

double maxdiff(const Mtr& m1, const Mtr& m2) {
    double diff = -INFINITY;
    if (m1.size() != m2.size() || m1[0].size() != m2[0].size())
        return diff;
    for (size_t i = 0; i < m1.size(); ++i)
        for (size_t j = 0; j < m1[0].size(); ++j){
            double d = fabs(m1[i][j] - m2[i][j]);
            diff = diff < d ? d : diff;
        }
    return diff; 
}

Mtr householder(const Vec& V){
    size_t n = size(V);
    Vec W = normalize(V);
    Mtr H = sum(E(n), mul(toCol(W), toRow(W)), 1, -2);
    return H;
}

pair<Mtr, Mtr> _separateLU(const Mtr& M){
    if (!issquare(M))
        throw invalid_argument("Только квадратные матрцы!");
    size_t n = M.size();
    Mtr L = fill(0.0, n,n);
    Mtr U = fill(0.0, n,n);

    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j){
            if (j > i)
                U[i][j] = M[i][j];
            else if (j < i)
                L[i][j] = M[i][j];
            else{
                L[i][j] = 1.0;
                U[i][j] = M[i][j];
            }
        }
    return make_pair(L, U);
}

Mtr apply_row_permutation(const Mtr& M, const vector<int>& perm){
    int n = M.size();
    Mtr permuted(n, Vec(n));
    for (int i = 0; i < n; i++) {
        permuted[i] = M[perm[i]];
    }
    return permuted;
}

Vec solveLinearEquation(const Mtr& L,
                        const Mtr& U,
                        const std::vector<int>& permutation,
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

// Функция для поиска максимального элемента в столбце и его строки
int _findMaxElementRow(const Mtr& M, int startRow, int mtr_size) {
    double maxVal = abs(M[startRow][startRow]);
    int maxRow = startRow;
    
    for (int i = startRow + 1; i < mtr_size; i++) {
        if (abs(M[i][startRow]) > maxVal) {
            maxVal = abs(M[i][startRow]);
            maxRow = i;
        }
    }
    return maxRow;
}
Mtr _eliminationStepForColumn(const Mtr& M, int col_num) {
    int mtr_size = M.size();
    Mtr result = M;

    for (int row_num = col_num + 1; row_num < mtr_size; ++row_num) {
        result[row_num][col_num] /= result[col_num][col_num];
        for (int j = col_num + 1; j < mtr_size; ++j) {
            result[row_num][j] -= result[row_num][col_num] * result[col_num][j];
        }
    }
    return result;
}

// Основная функция для разложения матрицы на L и U с перестановкой элементов
LU_result LU_decomposition(Mtr M) {
    int mtr_size = M.size();

    // Вектор для хранения текущей перестановки элементов
    vector<int> perm;

    // Инициализация начальной перестановки [0,1,...,(n-1)]
    for (int i = 0; i < mtr_size; i++)
        perm.push_back(i);


    for (int col_num=0; col_num<mtr_size; ++col_num){
        int maxRow = _findMaxElementRow(M,col_num,mtr_size);
        if (M[col_num][col_num] == 0){
            cerr << "Невозможно диагонализовать матрицу! Нулевой элемент на диагонали. Пропуск шага.\n";
            continue;
        }

        if(maxRow!=col_num){
            for(int j=col_num;j<mtr_size;++j)
                swap(M[col_num][j], M[maxRow][j]);
            swap(perm[col_num], perm[maxRow]);
        }
        M = _eliminationStepForColumn(M, col_num);
    }
    auto [L, U] = _separateLU(M);
    return {L, U, perm};
}

Mtr diag(const Vec& V){

    Mtr res = fill(0, size(V), size(V));
    for (size_t i=0; i<size(V); i++)
        res[i][i] = V[i];
    return res;
}

double maxdiff(const Vec& v1, const Vec& v2) {
    double diff = -INFINITY;
    if (v1.size() != v2.size())
        return diff;
    for (size_t i = 0; i < v1.size(); ++i){
        double d = fabs(v1[i] - v2[i]);
        diff = diff < d ? d : diff;
    }
    return diff; 
}

Vec div(const Vec& v1, const Vec& v2, double eps){
    if (v1.size() != v2.size())
        throw invalid_argument("Операция применима только к векторам одинаковой длины");
    Vec res;
    for (size_t i = 0; i < v1.size(); i++){
        if (fabs(v2[i]) >= eps)
            res.push_back(v1[i]/v2[i]);
        else
            cerr << "Одно из значений в знаменателе оказалось меньше порогового, сокращаем вектор результата\n";
    }        
    return res;
}

double max(const Vec& v){
    double res = -INFINITY;
    for (auto x : v)
        res = (x > res) ? x : res;
    return res;
}

double min(const Vec& v){
    double res = INFINITY;
    for (auto x : v)
        res = (x < res) ? x : res;
    return res;
}

double mean(const Vec& v){
    return sum(v)/size(v);
}



// Mtr householder(const Vec& w0){
    
//     // auto fill = [](double val, size_t rows, size_t cols) -> Mtr{
//     //     return Mtr(rows, Vec(cols, val));
//     // };

//     // auto toCol = [&](const Vec& v) -> Mtr{
//     //     Mtr alloc = fill(0, size(v), 1);
//     //     for (size_t i=0; i< size(v); i++)
//     //         alloc[i][0] = v[i];
//     //     return alloc;
//     // };

//     // auto toRow = [](const Vec& v) -> Mtr{
//     //     return Mtr(1, v);
//     // };

//     size_t n = size(w0);
//     Vec w = normalize(w0);
//     Mtr Q = sum(E(n), mul(toCol(w), toRow(w)), 1, -2);
//     return Q;
// }

Mtr generateRndSymPos(int n, double cond){
    Vec w0 = randVec(n);
    Mtr H = householder(w0);
    Vec eigens = {cond};
    for (int i=0; i<n-2; i++) eigens.push_back(rand_double() + 1.0);
    eigens.push_back(1);
    
    Mtr D = diag(eigens);
    Mtr A = mul(mul(T(H), D), H);
    return A;
}