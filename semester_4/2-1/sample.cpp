#include "LinearAlgebra.hpp"

int main(){
    double eps = 1e-1;
    int N = 4;
    
    Vec eigens = {5, -2, 1, 3};
    Mtr B = {
            { 6,  2, -1,  5},
            { 6,  2,  0, -2},
            {-3, -2,  7, -4},
            {-4,  7, -1,  1},
    };
    Mtr D = diag(eigens);
    Mtr A = mul(inv(B), mul(D, B));

    auto [nsteps_naive, eig_naive] = naive_PM(A, eps);
    double err_naive = fabs(eig_naive - max(eigens));

    auto [nsteps_normed, eig_normed] = normed_PM(A, eps);
    double err_normed = fabs(eig_normed - max(eigens));

    double mu = euclideanNorm(A);
    Mtr C = sum(A, E(N), 1, -mu);
    auto [nsteps_normed_C, eig_normed_C] = normed_PM(C, eps);
    double eig_normed_min = eig_normed_C + mu;
    double err_normed_min = fabs(eig_normed_min - min(eigens));



    print(A);
    cout << "det(A): " << det(A) << "\n\n";

    cout << "Naive PM:" << '\n'
        << "Max eigenval: " << eig_naive << "\n"
        << "Number of steps: " << nsteps_naive << "\n"
        << "eps: " << eps << "\n"
        << "err: " << err_naive << "\n"
        << "\n";

    cout << "Normed PM:" << '\n'
        << "eps: " << eps << "\n"
        << "Max eigenval: " << eig_normed << "\n"
        << "Number of steps: " << nsteps_normed << "\n"
        << "err_max: " << err_normed << "\n"
        << "Min eigenval: " << eig_normed_min << "\n"
        << "Number of steps min: " << nsteps_normed_C << "\n"
        << "err_min: " << err_normed_min << "\n"
        << "\n";

    return 0;
}