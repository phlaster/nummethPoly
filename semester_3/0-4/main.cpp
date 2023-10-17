// Решение СЛАУ прямыми методами
#include "LinearAlgebra.hpp"

int main() {
    size_t n;
    cout << "Введите размер системы уравнений: ";
    cin >> n;
    n = n > 9 ? 9 : n;
    
    Vec r = randVec(n);
    print("\nСгенерирован случайный вектор:");
    print(r);

    print("Перед применением алгорима Хаусхолдера ветор будет нормирован:");
    print(normalize(r));

    Mtr A = householder(r);
    print("Матрица A из преобразования Хаусхолдера:");
    print(A);

    print("Зададим вектор неизвестных x в уравнении: Ax=b");
    Vec x = randVec(n);
    print(x);

    print("Вычислим вектор b и запишем систему целиком:");
    Vec b = mul(A, x);
    print(A, x, b);

    print("Для решения системы нужно произвести LU-разложение матрицы A:");
    auto [L, U, perm] = LUP(A); //LU_decomposition(A);
    print({L, U, perm});

    print("Подтвердим верность расчётов, переможнив L*U и подставив строки в правильном порядке:");
    // Mtr reverse_LU = apply_row_permutation(mul(L, U), perm);
    vInt inv_perm = inversePermutation(perm);
    Mtr reverse_LU = mul(permutationMatrix(inv_perm), mul(L, U));
    print(reverse_LU);

    print("Элементы реконструированной матрицы отличаются от A не более, чем на:");
    print(maxdiff(A, reverse_LU));

    print("Вычислим вектор неизвестных, основываясь на LU-разложении A:");
    Vec x_LU = solveLinearEquation(L, U, perm, b);
    print(x_LU);

    print("Вектор невязки найденного решения:");
    Vec res = residual(A, x_LU, b);
    print(res);

    print("2-норма вектора невязки:");
    print(euclideanNorm(res));

    return 0;
}