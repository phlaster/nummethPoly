#include "headers/HEADER.hpp"

int main(int argv, char **argc){
    srand((unsigned)time(NULL));
    int msize;
    float cond, sparsity, threshold;
    if (argv < 5){
        msize = 10;
        cond = 5;
        sparsity = 1;
        threshold = 0;
    } else if (argv == 5){
        msize = stoi(argc[1]);
        cond = stof(argc[2]);
        sparsity = stof(argc[3]);
        threshold = stof(argc[4]);
    } else {
        print("Слишком много аргументов, завершение!");
        exit(1);
    }
    cout << "msize="<<msize<<"\n"
        << "cond="<<cond<<"\n"
        << "sparsity="<<sparsity<<"\n"
        << "threshold="<<threshold<<endl;

    spMtr A = generateRndSymPos(msize, cond, sparsity);
    Vec x = randVec(msize);
    Vec b = A*x;
    print(to_dense(A), x, b);

    
    auto [_x, nsteps] = pcg(A, b);
    cout << "Шагов: " << nsteps << "\nНорма ошибки: ";
    print(euclideanNorm(_x - x));

    
    spMtr Chol = chol(A, threshold);
    print(Chol);

    print(Chol*T(Chol) - A);

    
    auto [_x2, nsteps2] = pcg(A, b, Chol);
    cout << "Шагов: " << nsteps2 << "\nНорма ошибки: ";
    print(euclideanNorm(_x2 - x));

    return 0;
}