// Повторение решения из примера
#include "LinearAlgebra.hpp"

void routine(const Mtr& A, const Vec& b){

    print(A, {"x1", "x2", "x3"}, b);
    print("Определитель:");
    print(det(A));
    
    

    print("Для решения системы нужно произвести LU-разложение матрицы:");
    auto [L, U, perm] = LUP(A);
    print({L, U, perm});

    print("Соберём матрицу назад:");
    vInt inv_perm = inversePermutation(perm);
    Mtr reverse_LU = mul(permutationMatrix(inv_perm), mul(L, U));
    print(reverse_LU);

    print("Вычислим вектор неизвестных, основываясь на LU-разложении:");
    Vec x_LU = solveLinearEquation(L, U, perm, b);
    print(x_LU);

    print("Вектор невязки найденного решения:");
    Vec res = residual(A, x_LU, b);
    print(res);

    print("2-норма вектора невязки:");
    print(euclideanNorm(res));
}

int main() {
    Mtr A = {
        {1,2,3},
        {4,5,6},
        {7,8,9}
    };
    Vec b = randVec(3, -10, 10);


    print("Проверка матрицы A:");
    try {
        routine(A, b);
    } catch (const runtime_error& e) {
        cerr << "Ошибка: " << e.what() << endl;
    }


    print("Проверка матрицы B:");
    try {
        routine(mul(1e8, A), b);
    } catch (const runtime_error& e) {
        cerr << "Ошибка: " << e.what() << endl;
    }

    return 0;
}
