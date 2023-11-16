#include "headers/HEADER.hpp"
#include "headers/Sparse.hpp"
#include <cmath>

void clargs_input(int& argv, int& msize, float& cond, float& density, float& threshold, char **argc){
    if (argv < 5){
        msize = 10;
        cond = 5;
        density = 1;
        threshold = 0;
    } else if (argv == 5){
        msize = stoi(argc[1]);
        cond = stof(argc[2]);
        density = stof(argc[3]);
        threshold = stof(argc[4]);
    } else {
        print("Слишком много аргументов, завершение!");
        exit(1);
    }
    cout << "msize="<<msize<<"\n"
        << "cond="<<cond<<"\n"
        << "density="<<density<<"\n"
        << "threshold="<<threshold<<endl;

}

void ad_hoc_check(){
    print("Ad hoc check!\n\n");
    spMtr AA({
        {3.90235,   0.776848,  -0.653871,  0.383767,  0.352175},
        {0.776848,  1.23975,    0.0,       0.0,       0.0},
        {-0.653871,  0.0,        2.66638,   0.0,       0.0},
        {0.383767,  0.0,        0.0,       2.8985,    0.0},
        {0.352175,  0.0,        0.0,       0.0,       4.29302}}
    );

    print(chol(AA));
    print(ichol(AA, 0));
    print(ichol(AA, 0.5));
    print(ichol(AA, 50));
}
void cholesky1(){
    print("Full Cholesky check\n\n");
    spMtr A = generateRndSymPos(10, 10, 1);
    spMtr C1 = chol(A);
    print(C1);
    print(maxabs(A - C1*T(C1)));
}
void chol_pentadiag(size_t n){
    print("Cholesky 5-diag check:");
    spMtr F = block5diag(n);
    print(F);
    print(chol(F));
    print(ichol(F));
    print(ichol(F, 0.5));
}
void cg_check(int argv, char **argc){
    print("Conjugated grad-s:");

    int msize; float cond, sparsity, threshold;
    
    clargs_input(argv, msize, cond, sparsity, threshold, argc);

    spMtr A = generateRndSymPos(msize, cond, sparsity);
    Vec x = randVec(msize);
    Vec b = A*x;
    
    auto [_x, nsteps] = cg(A, b);
    cout << "Шагов: " << nsteps << "\nНорма ошибки: ";
    print(euclideanNorm(_x - x));
    if (nsteps>msize){
        cerr << "Метод не сошёлся за n шагов!" << endl;
    }
}
void pcg_check(int argv, char **argc){
    print("PRECONDITIONED Conjugated grad-s:");

    int msize; float cond, density, threshold;
    clargs_input(argv, msize, cond, density, threshold, argc);

    spMtr A = generateRndSymPos(msize, cond, density);
    spy(A);
    Vec x = randVec(msize);
    Vec b = A*x;
    spMtr C = chol(A);
    spMtr C0 = ichol(A, 0);
    spMtr C1 = ichol(A, threshold);
    spMtr C2 = ichol(A, INFINITY);

    auto [_, n] = cg(A, b);                     double err_ = euclideanNorm(_ - x);
    auto [_x, nsteps] = pcg(A, C, b);           double err  = euclideanNorm(_x - x);
    auto [_x0, nsteps0] = pcg(A, C0, b);        double err0 = euclideanNorm(_x0 - x);
    auto [_x1, nsteps1] = pcg(A, C1, b);        double err1 = euclideanNorm(_x1 - x);
    auto [_x2, nsteps2] = pcg(A, C2, b);        double err2 = euclideanNorm(_x2 - x);
    auto [_x3, nsteps3] = pcg(A, E(msize), b);  double err3 = euclideanNorm(_x3 - x);

    cout << "Без предобуславливателя: " << n       << "   err="<<err_<<endl;
    cout << "Полное разложение      : " << nsteps  << "   err="<<err <<endl;
    cout << "Разложение с theta=0   : " << nsteps0 << "   err="<<err0<<endl;
    cout << "Разложение с заданным  : " << nsteps1 << "   err="<<err1<<endl;
    cout << "Разложение с theta=inf : " << nsteps2 << "   err="<<err2<<endl;
    cout << "Единичный предобусл.   : " << nsteps3 << "   err="<<err3<<endl;
}
void comp_bs(int n){
    print("Comparing estimated and exact rhs vectors:");
    spMtr F = block5diag(n, n);
    Vec x = exact_x(n);
    Vec b_ = estimate_b(n);
    Vec b = F*x;
    print(euclideanNorm(b-b_));
}
void pentadiag(size_t n){
    Vec x = exact_x(n);

    spMtr F = block5diag(n);
    Vec b = estimate_b(n);
    auto [_x, n_steps] = cg(F, b);
    print(euclideanNorm(x - _x));
}
void pentadiag_check(int n){
    print("Model task:");
    spMtr F = block5diag(n);
    Vec x = exact_x(n);
    Vec b = estimate_b(n);
    // Vec b = F*x;

    spMtr C = chol(F);
    spMtr C0 = ichol(F);
    spMtr C1 = ichol(F, 0.5);
    spMtr C2 = ichol(F, INFINITY);
    // print(to_dense(C));
    spy(C);spy(C0);spy(C1);spy(C2);

    auto [_, nsteps__] = cg(F, b);           double err_ = euclideanNorm(_ - x);
    auto [_x, nsteps] = pcg(F, C, b);        double err  = euclideanNorm(_x - x);
    auto [_x0, nsteps0] = pcg(F, C0, b);     double err0 = euclideanNorm(_x0 - x);
    auto [_x1, nsteps1] = pcg(F, C1, b);     double err1 = euclideanNorm(_x1 - x);
    auto [_x2, nsteps2] = pcg(F, C2, b);     double err2 = euclideanNorm(_x2 - x);
    auto [_x3, nsteps3] = pcg(F, E(n*n), b); double err3 = euclideanNorm(_x3 - x);

    cout << "Без предобуславливателя: " << nsteps__<< "   err="<<err_<<endl;
    cout << "Полное разложение      : " << nsteps  << "   err="<<err <<endl;
    cout << "Разложение с theta=0   : " << nsteps0 << "   err="<<err0<<endl;
    cout << "Разложение с theta=0.5 : " << nsteps1 << "   err="<<err1<<endl;
    cout << "Разложение с theta=inf : " << nsteps2 << "   err="<<err2<<endl;
    cout << "Единичный предобусл.   : " << nsteps3 << "   err="<<err3<<endl;
}
void ichol_test0(int n){
    spMtr A = generateRndSymPos(n, 10, 0.4);
    spMtr C = chol(A);
    spMtr C0 = ichol(A);
    spMtr C1 = ichol(A, 0.5);
    spy(A); spy(C); spy(C0); spy(C1);
}
void ichol_test5diag(int n){
    spMtr A = block5diag(n);
    spMtr C = chol(A);
    spMtr C0 = ichol(A);
    spMtr C1 = ichol(A, 0.5);
    spy(A); spy(C); spy(C0); spy(C1);
}


void task_theta1(size_t n){
    print("Model task:");
    spMtr F = block5diag(n);
    Vec x = exact_x(n);
    Vec b = estimate_b(n);
    // spMtr C = chol(F);
    auto [_x, nsteps] = cg(F, b, 1e-5);
    double err = euclideanNorm(_x-x);
    cout << "CG: steps=" << nsteps << endl;
    cout << "CG: err=" << err << endl;
    cout << "theta,n,err_i"<<endl;
    for(double theta=1e-5; theta<10; theta*=2){
        spMtr C_theta = ichol(F, theta);
        // spy(C_theta);
        auto [x_i, nsteps_i] = pcg(F, C_theta, b, 1e-5);
        double err_i = euclideanNorm(x_i-x);
        cout << theta<< "," <<nsteps_i << "," << err_i << endl;
    }
}

int main(int argv, char **argc){
    srand((unsigned)time(NULL));

    // ad_hoc_check();
    // cholesky1();
    // chol_pentadiag(stoi(argc[1]));
    // cg_check(argv, argc);
    // pcg_check(argv, argc);
    // comp_bs(stoi(argc[1]));
    // pentadiag(stoi(argc[1]));
    // pentadiag_check(stoi(argc[1]));
    // ichol_test0(stoi(argc[1]));
    // ichol_test5diag(stoi(argc[1]));
    task_theta1(stoi(argc[1]));   

    return 0;
}