#include "LinearAlgebra.hpp"
#include <stdexcept>

void find_max_and_swap(Mtr& M, vInt& perm, int i) {
    size_t N = M.size();
    double maxVal = 0.0;
    int maxRow = i;

    for (int j = i; j < N; j++) {
        double absVal = fabs(M[j][i]);
        if (absVal > maxVal) {
            maxVal = absVal;
            maxRow = j;
        }
    }
    if (fabs(maxVal) < 1e-16){
        cerr << "i = " << i << ", pivot = " << maxVal << endl;
        throw runtime_error("Невозможно найти максимальнный элемент в нулевом столбце!");
    }
    if (maxRow != i){
        swap(M[i], M[maxRow]);
        swap(perm[i], perm[maxRow]);
    }
}

void div_under_main_diag(Mtr& M, int i) {
    size_t N = M.size();
    double pivot = M[i][i];
    if (fabs(pivot) < 1e-16)
        throw runtime_error("Опорный элемент близок к нулю! Невозможно прозвести шаг деления!");
    for (size_t j = i + 1; j < N; j++)
        M[j][i] /= pivot;
}

void subtract_product(Mtr& M, int i) {
    size_t N = M.size();
    for (size_t row = i + 1; row < N; row++)
        for (size_t col = i + 1; col < N; col++)
            M[row][col] -= M[row][i] * M[i][col];
}

pair<Mtr, Mtr> split_LUP(Mtr& U) {
    size_t N = U.size();
    Mtr L(N, Vec(N, 0.0));
    for (size_t i = 0; i < N; i++) {
        L[i][i] = 1.0;
        for (size_t j = 0; j < i; j++) {
            L[i][j] = U[i][j];
            U[i][j] = 0.0;
        }
    }

    return make_pair(L, U);
}

LU_result LUP(Mtr M){
    if (!issquare(M))
        throw invalid_argument("Только квадратные матрицы!");
    
    size_t N = M.size();
    vInt perm(N); for (int i = 0; i < N; i++) perm[i] = i;
    cout << "Исходная:\n";
    print(M);
    for (size_t i=0; i<N; i++){ // По каждому ряду
        cout << "Столбец " << i << ":\n";
        find_max_and_swap(M, perm, i);
        cout << "После перестановки:\n";
        print(M);
        div_under_main_diag(M, i);
        cout << "После деления:\n";
        print(M);
        subtract_product(M, i);
        cout << "После вычитания:\n";
        print(M);
    }
    auto [L, U] = split_LUP(M);
    return {L, U, perm};
}

Mtr apply_row_permutation(const Mtr& M, const vInt& perm){
    int n = M.size();
    Mtr permuted(n, Vec(n));
    for (int i = 0; i < n; i++)
        permuted[i] = M[perm[i]];
    return permuted;
}

vInt inversePermutation(const vInt& permutation) {
    int n = permutation.size();
    vInt inverse(n);
    for (int i = 0; i < n; i++)
        inverse[permutation[i]] = i;
    return inverse;
}

Mtr permutationMatrix(const vInt& permutation) {
    int n = permutation.size();
    Mtr matrix(n, Vec(n, 0));
    for (int i = 0; i < n; i++)
        matrix[i][permutation[i]] = 1.0;
    return matrix;
}
