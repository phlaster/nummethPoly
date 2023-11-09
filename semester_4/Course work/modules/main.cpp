#include "headers/HEADER.hpp"
#include "headers/Sparse.hpp"

void clargs_input(int& argv, int& msize, float& cond, float& sparsity, float& threshold, char **argc){
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

}

int main(int argv, char **argc){
    srand((unsigned)time(NULL));
    int msize;
    float cond, sparsity, threshold;
    clargs_input(argv, msize, cond, sparsity, threshold, argc);

    // spMtr L({
    //     {-9,   0,   0,   0, 0},
    //     { 2,  -2,   0,   0, 0},
    //     {-10,   0,  -9,   0, 0},
    //     { 4,   3,   9,   2, 0},
    //     { 2,   9,  -3,  -9, 6}
    // });
    // spMtr U({
    //     {-6,   2,  -10,  -1,  -2},
    //     {0 , -9 ,  -5 , -8 , -6},
    //     {0 ,  0 ,  -4 ,  0 ,  8},
    //     {0 ,  0 ,   0 , -2 , -1},
    //     {0 ,  0 ,   0 ,  0 ,  8}
    // });
    // Vec x = {5,5,1,5,3};
    // Vec bL = L*x;
    // Vec bU = U*x;


    // I
    spMtr A = generateRndSymPos(msize, cond, sparsity);
    Vec x = randVec(msize);
    Vec b = A*x;
    // print(to_dense(A), x, b);
    
    auto [_x, nsteps] = pcg(A, b);
    cout << "Шагов: " << nsteps << "\nНорма ошибки: ";
    print(euclideanNorm(_x - x));
    
    spMtr Chol = chol(A, threshold);
    // print(Chol);

    cout << "Max = " << maximum(Chol*T(Chol) - A) << endl;
    
    auto [_x2, nsteps2] = pcg(A, Chol, b);
    cout << "Шагов: " << nsteps2 <<"\nНорма ошибки: ";
    print(euclideanNorm(_x2 - x));

    return 0;
}