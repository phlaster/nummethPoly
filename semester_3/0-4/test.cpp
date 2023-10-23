#include "LinearAlgebra.hpp"


int main(){
    Mtr A = {
        {2,2,2},
        {2,3,3},
        {2,3,3}
    };
    vInt perm = {0,1,2};
    print(A);
    find_max_and_swap(A, perm, 0);
    print(A);
    subtract_product(A, 0);
    print(A);
    print(perm);

    return 0;
}