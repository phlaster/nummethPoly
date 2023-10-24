#include "Functions.hpp"


int main() {
    //Можно запустить матрицу из примера в отчёте:
    Mtr A = {
        { 5,-1,-1, 0},
        {-1, 6, 1,-1},
        {-1, 1, 6, 1},
        { 0,-1, 1, 7}
    };
    Vec x = {-2, 0, 2, 1};
    Vec b = {-12, 3, 15, 9};

    cout << "x = "; print(x);
    
    PCG(A, b, -2, 100, true);
    return 0;
}