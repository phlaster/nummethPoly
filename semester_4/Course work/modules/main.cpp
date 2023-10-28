#include "HEADER.hpp"

int main(int argv, char **argc){
    srand((unsigned)time(NULL));

    int msize = stoi(argc[1]);
    float cond = stof(argc[2]);
    float sparsity = stof(argc[3]);
    float threshold = stof(argc[4]);

    spMtr A = generateRndSymPos(msize, cond, sparsity);
    spMtr B = A;
    cout << "A:" << endl;
    print(A);

    ___incomp_chol_tol(A, 0);
    cout << "chol(A):" << endl;
    print(A);
    cout << A.valCounter << endl;
    
    ___incomp_chol_tol(B, threshold);
    cout << "chol(A, delta="<<threshold<<"):" << endl;
    print(B);
    cout << B.valCounter << endl;
    
    // erase_above_diag(A, true);
    // chol(A);
    // print(A);
    // spMtr Ch = chol(A);
    // print(to_dense(Ch));
    // print(A - Ch* T(Ch));

    // Mtr A = generateRndSymPos(4,10);
    // print(A);
    // choleskyDecomposition(A);
    // print(A);
    return 0;
}