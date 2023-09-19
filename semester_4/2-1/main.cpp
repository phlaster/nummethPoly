#include "LinearAlgebra.hpp"

int main(){
    size_t N = 8;
    double delta = 1e-15;

    non_singular_needed:
        Mtr B = randMtr(N);
    if (det(B) < 1e-7) goto non_singular_needed; // && det(B) > 0
    Vec eigens = randVec(N);
    Mtr D = diag(eigens);
    Mtr A = mul(inv(B), mul(D, B));

    print("Random matrix B generated:");
    print(B);
    cout << "det(B): " << det(B) << "\n\n";
    print("Diagonal D:");
    print(D);
    print("A = B^-1 * D * B:");
    print(A);
    cout << "det(A): " << det(A) << "\n\n";


    cout << "Naive PM:" << '\n';
    auto [nsteps_naive, eig_naive] = naive_PM(A, delta);
    double err_naive = fabs(eig_naive - max(eigens));

    cout << "Max eigenval: " << eig_naive << "\n"
        << "Number of steps: " << nsteps_naive << "\n"
        << "delta: " << delta << "\n"
        << "err: " << err_naive << "\n\n";


    cout << "Normed PM:" << '\n';
    auto [nsteps_normed, eig_normed] = normed_PM(A, delta);
    double err_normed = fabs(eig_normed - max(eigens));

    cout << "Max eigenval: " << eig_normed << "\n"
        << "Number of steps: " << nsteps_normed << "\n"
        << "delta: " << delta << "\n"
        << "err: " << err_normed << "\n\n";


    double mu = euclideanNorm(A);
    Mtr C = sum(A, E(N), 1, -mu);

    cout << "Normed PM(min):" << '\n';
    auto [nsteps_normed_C, eig_normed_C] = normed_PM(C, delta);
    double eig_normed_min = eig_normed_C + mu;
    double err_normed_min = fabs(eig_normed_min - min(eigens));

    cout << "Min eigenval: " << eig_normed_min << "\n"
        << "Number of steps min: " << nsteps_normed_C << "\n"
        << "delta: " << delta << "\n"
        << "err: " << err_normed_min << "\n\n";

    return 0;
}