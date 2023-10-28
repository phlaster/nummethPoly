#include "LinearAlgebra.hpp"


int main(int argv, char **argc){
    int N = atoi(argc[1]);
    double cond = atof(argc[2]);

    double mu = cond * 1.1;
    Mtr A = generateRndSymPos(N, cond);
    Mtr _A_ = sum(A, E(N), 1, -mu);

    // print(A);
    // print(_A_);

    // cout << "N="<<N << "  cond="<<cond<<endl;
    cout << "eps,steps_big,err_big,steps_small,err_small\n";
    for (int power=-2; power>=-15; power--){
        double eps = pow(10, power);

        auto [steps_big, eig_big] = normed_PM(A, eps);
        auto [steps_small, Lambda] = normed_PM(_A_, eps);
        double eig_small = Lambda + mu;
        
        double err_big = fabs(eig_big-cond);
        double err_small = fabs(eig_small-1);

        cout <<
            eps << "," <<
            steps_big << "," <<
            err_big << "," <<
            steps_small << "," <<
            err_small << "\n";
    }

}