// Повторение решения из примера
#include "LinearAlgebra.hpp"
#include <stdexcept>

void routine(const Mtr& A, const Vec& b){

    print(A);
    print("Определитель:");
    print(det(A));
    
    try{
        double c = cond(A);
        print("Число обусловленности:");
        print(c);
    }
    catch (invalid_argument){
        print("Расчитать число обусловленности не удалось!\n\n");
    }
    

    print("Для решения системы нужно произвести LU-разложение матрицы:");
    auto [L, U, perm] = LU_decomposition(A);
    print({L, U, perm});

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
    routine(A, b);

    print("Проверка матрицы B:");
    routine(mul(1e8, A), b);

    


    return 0;
}
