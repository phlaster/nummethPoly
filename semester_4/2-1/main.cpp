#include "LinearAlgebra.hpp"

int main(){
    size_t N = 5;
    double delta = 1e-15;
    double max_eig = 12.5, min_eig = 1, cond = max_eig/min_eig;
    Mtr A = generateRndSymPos(N, cond);
    // print(A);
    cout << "det(A): " << det(A) << "\n\n";

    cout << "Normed PM:" << '\n';
    auto [nsteps_normed, eig_normed] = normed_PM(A, delta);
    double err_normed = fabs(eig_normed - max_eig);

    cout << "Max eigenval: " << eig_normed << "\n"
        << "Number of steps: " << nsteps_normed << "\n"
        << "delta: " << delta << "\n"
        << "err: " << err_normed << "\n\n";


    double mu = opnorm_1(A);
    Mtr C = sum(A, E(N), 1, -mu);

    cout << "Normed PM(min):" << '\n';
    auto [nsteps_normed_C, eig_normed_C] = normed_PM(C, delta);
    double eig_normed_min = eig_normed_C + mu;
    double err_normed_min = fabs(eig_normed_min - min_eig);

    cout << "Min eigenval: " << eig_normed_min << "\n"
        << "Number of steps min: " << nsteps_normed_C << "\n"
        << "delta: " << delta << "\n"
        << "err: " << err_normed_min << "\n\n";

    return 0;
}