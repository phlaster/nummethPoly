#include "headers/HEADER.hpp"
#include "headers/Sparse.hpp"
#include <cmath>
#include <string>

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
        {7.7 ,  4.36,  7.79 , 5.77,  0.4},
        {4.36,  7.26,  8.17 , 3.29,  0},
        {7.79,  8.17,  12.87, 7.23,  1.34},
        {5.77,  3.29,  7.23 , 7.87,  1.36},
        {0.4 ,     0,  1.34 , 1.36,  1.2}
    });

    print(chol(AA));
    print(ichol(AA));
    print(ichol(AA, 0.2));
    print(ichol(AA, 50));
}
void cholesky1(int n){
    print("Full Cholesky check");
    for (int i=0; i<=10; i++){
        double density = i/10.;
        spMtr A = generateRndSymPos(n, 10, density);
        spMtr C = chol(A);
        spMtr iC = ichol(A);
        // print(C);
        cout << "density="<<density<<", err_full=" << maxabs(A - C*T(C))<<", err_zero=" << maxabs(A - iC*T(iC)) << endl;
    }
}
void cholesky2(int n){
    print("Full Cholesky check");
    for (int i=0; i<=10; i++){
        double theta = i/10.;
        spMtr A = block5diag(n);
        spMtr C = chol(A);
        spMtr iC = ichol(A, theta);
        // print(C);
        cout << "theta="<<theta<<", err_full=" << maxabs(A - C*T(C))<<", err_zero=" << maxabs(A - iC*T(iC)) << endl;
    }
}
void chol_pentadiag(size_t n){
    print("Cholesky 5-diag check:");
    spMtr F = block5diag(n);
    print(F);
    print(chol(F));
    print(ichol(F));
    print(ichol(F, 1));
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
        cerr << "Метод не сошёлся за 2n шагов!" << endl;
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

    auto [_, nsteps__] = cg(F, b, 1e-5);           double err_ = euclideanNorm(_ - x);
    auto [_x, nsteps] = pcg(F, C, b, 1e-5);        double err  = euclideanNorm(_x - x);
    auto [_x0, nsteps0] = pcg(F, C0, b, 1e-5);     double err0 = euclideanNorm(_x0 - x);
    auto [_x1, nsteps1] = pcg(F, C1, b, 1e-5);     double err1 = euclideanNorm(_x1 - x);
    auto [_x2, nsteps2] = pcg(F, C2, b, 1e-5);     double err2 = euclideanNorm(_x2 - x);
    auto [_x3, nsteps3] = pcg(F, E(n*n), b, 1e-5); double err3 = euclideanNorm(_x3 - x);

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
    spy(A); spy(C); spy(C0);
    for (double t=1e-2;t<1;t+=t){
        spMtr C1 = ichol(A, t);
        cout << "theta=" << t << endl;
        spy(C1);
    }
}

void task_0_theta(size_t n){
    print("Rand sym:");

    spMtr A = generateRndSymPos(n*n, 103, 0.033690476);
    Vec x = randVec(n*n);
    Vec b = A*x;

    auto [x_cg, steps_cg] = cg(A, b, 1e-5);
    double err_cg = euclideanNorm(x_cg - x);
    cout << "CG: steps=" << steps_cg << endl;
    cout << "CG: err=" << err_cg << endl;
    cout << "PCG:" << endl;
    cout << "theta,steps,err" << endl;
    for(double theta=1e-7; theta<10; theta *= 1.4){
        spMtr iC = ichol(A, theta);
        auto [x_i, steps_i] = pcg(A, iC, b, 1e-5);
        double err_i = euclideanNorm(x_i - x);
        cout << theta << "," << steps_i << "," << err_i << endl;
    }
}

void task_1_theta(size_t n){
    print("Model task:");
    spMtr F = block5diag(n);
    Vec x = exact_x(n);
    Vec b = estimate_b(n);

    auto [x_cg, steps_cg] = cg(F, b, 1e-5);
    double err_cg = euclideanNorm(x_cg - x);
    cout << "CG: steps=" << steps_cg << endl;
    cout << "CG: err=" << err_cg << endl;
    cout << "PCG:" << endl;
    cout << "theta,steps,err" << endl;
    for(double theta=1e-7; theta<10; theta *= 1.4){
        spMtr iC = ichol(F, theta);
        auto [x_i, steps_i] = pcg(F, iC, b, 1e-5);
        double err_i = euclideanNorm(x_i - x);
        cout << theta << "," << steps_i << "," << err_i << endl;
    }
}

void getxx(size_t n){
    spMtr F = block5diag(n);
    Vec x = exact_x(n);
    Vec b = estimate_b(n);

    spMtr iC = ichol(F);
    auto [x_i, steps_i] = pcg(F, iC, b, 1e-5);

    print(x);
    print(x_i);
}

int main(int argv, char **argc){
    srand((unsigned)time(NULL));

    // ad_hoc_check();
    // cholesky1(stoi(argc[1]));
    // cholesky2(stoi(argc[1]));

    // chol_pentadiag(stoi(argc[1]));
    // cg_check(argv, argc);
    // pcg_check(argv, argc);
    // comp_bs(stoi(argc[1]));
    // pentadiag(stoi(argc[1]));
    // pentadiag_check(stoi(argc[1]));
    // ichol_test0(stoi(argc[1]));
    ichol_test5diag(stoi(argc[1]));
    
    
    // task_0_theta(stoi(argc[1]));
    // task_1_theta(stoi(argc[1]));

    // getxx(stoi(argc[1]));

    return 0;
}