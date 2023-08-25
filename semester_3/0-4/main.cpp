// Решение СЛАУ прямыми методами
#include "LinearAlgebra.hpp"
#include <algorithm>

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

    Mtr Q = householder(r);
    print("Матрица Q из преобразования Хаусхолдера:");
    print(Q);
    
    print("Число обусловленности Q:");
    print(cond(Q));

    print("Зададим вектор неизвестных x в уравнении: Qx=b");
    Vec x = randVec(n);
    print(x);

    print("Вычислим вектор b и запишем систему целиком:");
    Vec b = mul(Q, x);
    print(Q, x, b);

    print("Для решения системы нужно произвести LU-разложение матрицы Q:");
    auto [L, U, perm] = LU_decomposition(Q);
    print({L, U, perm});

    print("Подтвердим верность расчётов, переможнив L*U и подставив строки в правильном порядке:");
    Mtr reverse_LU = apply_row_permutation(mul(L, U), perm);
    print(reverse_LU);

    print("Элементы реконструированной матрицы отличаются от Q не более, чем на:");
    print(maxdiff(Q, reverse_LU));

    print("Вычислим вектор неизвестных, основываясь на LU-разложении Q:");
    Vec x_LU = solveLinearEquation(L, U, perm, b);
    print(x_LU);

    print("Вектор невязки найденного решения:");
    Vec res = residual(Q, x_LU, b);
    print(res);

    print("2-норма вектора невязки:");
    print(euclideanNorm(res));

    return 0;
}