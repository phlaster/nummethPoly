#include "headers/HEADER.hpp"

spMtr chol(const spMtr& M) {
    int n = M.cols;
    spMtr L(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            if (i == j) {
                double sum = 0.0;
                for (int k = 0; k < j; ++k) {
                    sum += pow(L.get(j,k), 2);
                }
                L.set(sqrt(M.get(j,j) - sum), j, j);
            } else {
                double sum = 0.0;
                for (int k = 0; k < j; ++k) {
                    sum += L.get(i,k) * L.get(j,k);
                }
                L.set((M.get(i,j) - sum) / L.get(j,j), i,j);
            }
        }
    }
    return L;
}
spMtr ichol(const spMtr& A, double theta) {
    int N = A.cols;
    spMtr L(N, N);

    for (int I = 0; I < N; I++) {
        double S = A.get(I, I);

        for (int IP = 0; IP < I; IP++) {
            S -= L.get(I, IP) * L.get(I, IP);
        }

        L.set(sqrt(S), I, I);

        for (int J = I + 1; J < N; J++) {
            S = A.get(J, I);

            for (int IP = 0; IP < I; IP++) {
                S -= L.get(J, IP) * L.get(I, IP);
            }
            if (fabs(S) > theta){
                L.set(S / L.get(I, I), J, I);
            }
        }
    }

    return L;
}