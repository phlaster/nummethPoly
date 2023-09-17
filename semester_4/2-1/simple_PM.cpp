#include "LinearAlgebra.hpp"

pair<int, double> simple_power_method(const Mtr& A, double eps=1e-7){
    size_t N = size(A);
    Vec y_i, y = randVec(N);
    int nsteps = 0;
    double lambda_i, lambda = INFINITY;

    for(;;) {
        y_i = mul(A, y);
        lambda_i = y_i[0] / y[0];

        if (fabs(lambda - lambda_i) <= eps)
            return make_pair(nsteps, lambda);
        
        y = y_i;
        lambda = lambda_i;
        nsteps++;
    }
}


int main(){
    size_t N = 4;
    double eps = 1e-3;

    non_singular_needed:
        Mtr B = randMtr(N);
    if (fabs(det(B)) < 1e-7) goto non_singular_needed;
    Vec eigens = randVec(N);
    Mtr D = diag(eigens);
    Mtr A = mul(inv(B), mul(D, B));
    auto [nsteps, eig] = simple_power_method(A, eps);
    double err = fabs(eig - max(eigens));


    print("Random matrix B generated:");
    print(B);
    print("det(B):");
    print(det(B));
    print("Diagonal D:");
    print(D);
    print("A = B^-1 * D * B:");
    print(A);
    cout
        << "Max eigenval: " << eig << "\n"
        << "Number of steps: " << nsteps << "\n"
        << "eps: " << eps << "\n"
        << "err: " << err << "\n"
        << "\n";

    return 0;
}