// Повторение решения из примера
#include "LinearAlgebra.hpp"

int main() {
    Mtr A = {
        {-0.210, -0.155,  -2.15,   -0.129},
        {-0.155,  2.001,   0.058,   0.16},
        {-2.150,  0.058,   2.637,   3.525},
        {-0.129,  0.160,   3.525,   3.572}
    };

    print("Матрица A использованная в отчёте в качестве примера:");
    print(A);

    print("Вектор неизвестных x, использованный в отчёте: Ax=b");
    Vec x = {1, 2, 3, 4};
    print(x);

    print("Вычислим вектор b и запишем систему целиком:");
    Vec b = mul(A, x);
    print(A, x, b);

    print("Для решения системы нужно произвести LU-разложение матрицы Q:");
    auto [L, U, perm] = LUP(A);
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
