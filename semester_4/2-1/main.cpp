#include "LinearAlgebra.hpp"

int main(){
    size_t N = 4;
    double eps = 1e-2;

    non_singular_needed:
        Mtr B = randMtr(N);
    if (det(B) < 1e-7) goto non_singular_needed; // && det(B) > 0
    Vec eigens = randVec(N);
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



    print("Random matrix B generated:");
    print(B);
    cout << "det(B): " << det(B) << "\n\n";
    print("Diagonal D:");
    print(D);
    print("A = B^-1 * D * B:");
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