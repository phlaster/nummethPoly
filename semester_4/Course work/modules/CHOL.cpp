#include "headers/HEADER.hpp"

spMtr chol(const spMtr& A) {
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

            L.set(S / L.get(I, I), J, I);
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
            if (S==0.0) continue;

            for (int IP = 0; IP < I; IP++) {
                S -= L.get(J, IP) * L.get(I, IP);
            }

            double toset = S / L.get(I, I);
            if (fabs(toset) <= theta) continue;

            L.set(toset, J, I);
        }
    }
    return L;
}