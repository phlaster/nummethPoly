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
    
    // auto [_x, nsteps] = pcg(A, b);
    // cout << "Шагов: " << nsteps << "\nНорма ошибки: ";
    // print(euclideanNorm(_x - x));
    auto [mu, sigma] = mean_std(A);
    double maxval = maximum(A);
    auto [_x1, nsteps1] = pcg(A, b);
    double err1 = euclideanNorm(_x1-x);
    cout << "theta,nsteps2,err2,msize="<<msize<<",cond="<<cond<<",nsteps1="<<nsteps1<<",err1="<<err1<<",mu="<<mu<<",sigma="<<sigma<<",maxval="<<maxval<<endl;
    for (double theta=1e-9; theta <= maxval; theta *= 1.7){
        spMtr Chol = chol(A, theta);
        auto [_x2, nsteps2] = pcg(A, Chol, b);
        double err2 = euclideanNorm(_x2-x);
        cout<<theta<<","<<nsteps2<<","<<err2<<endl;
    }


    // cout << "Max = " << maximum(Chol*T(Chol) - A) << endl;
    
    // cout << "Шагов: " << nsteps2 <<"\nНорма ошибки: ";
    // print(euclideanNorm(_x2 - x));

    return 0;
}