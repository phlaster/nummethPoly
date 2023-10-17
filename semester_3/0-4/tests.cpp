#include "LinearAlgebra.hpp"


int main(){
    Mtr A = {
        {1,2,3},
        {4,5,6},
        {7,8,9}
    };
    vInt perm = {0,1,2};
    print(A);
    find_max_and_swap(A, perm, 0);
    print(A);
    print(perm);

    return 0;
}