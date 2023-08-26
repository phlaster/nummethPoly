// Повторение решения из примера
#include "LinearAlgebra.hpp"

int main() {
    Mtr Q = {
        {0.21429,  -0.0370467,   -0.675735,   -0.704336},
        {-0.0370467,    0.998253,  -0.0318613,  -0.0332099},
        {-0.675735,  -0.0318613,    0.418847,   -0.605751},
        {-0.704336,  -0.0332099,   -0.605751,     0.36861}
    };

    print("Матрица Q использованная в отчёте в качестве примера:");
    print(Q);
    
    print("Число обусловленности Q:");
    print(cond(Q));

    print("Вектор неизвестных x, использованный в отчёте: Qx=b");
    Vec x = {0.131611, 0.442564,0.487704,0.956933};
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
