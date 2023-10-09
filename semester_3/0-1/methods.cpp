#include "functions.hpp"

ans bisection(double (*f)(double), double exact, Vec lims, double eps){
    double a = lims[0], b = lims[1], c;
    int nSteps = 1; 
    for(;;){
        c = (a+b)/2;
        if (fabs(a-c) <= eps)
            return ans({c, fabs(exact-c), nSteps});
        if (f(a) * f(c) > 0)
            b = c;
        else
            a = c;
        nSteps++;
    }
}

ans newton(double (*f)(double), double exact, Vec lims, double eps){
    auto d1 = (f==f1) ? f1d1 : f2d1;
    auto d2 = (f==f1) ? f1d2 : f2d2;

    double a = lims[0], b = lims[1];
    double x0 = (f(a)*d2(a) > 0) ? a : b;
    
    double alpha = min(fabs(d1(a)), fabs(d1(b)));
    double beta  = max(fabs(d2(a)), fabs(d2(b)));

    cout << "alpha="<<alpha<<", beta="<<beta<<endl;

    int nSteps = 1;
    for (;;){
        double x_n = x0 - f(x0)/d1(x0);
        double delta = fabs(x_n-x0);
        if (beta/alpha/2 * delta * delta <= eps)
            return ans({x_n, fabs(exact-x_n), nSteps});
        x0 = x_n;
        nSteps++;
    }
}

ans newton_broken(double exact, Vec lims, double x0, double eps){
    double a = lims[0], b = lims[1];
    
    double alpha = min(fabs(f1d1(a)), fabs(f1d1(b)));
    double beta  = max(fabs(f1d2(a)), fabs(f1d2(b)));

    cout << "alpha="<<alpha<<", beta="<<beta<<endl;

    int nSteps = 1;
    for (;;){
        double x_n = x0 - f1(x0)/f1d1(x0);
        double delta = fabs(x_n-x0);
        if (beta/alpha/2 * delta * delta <= eps)
            return ans({x_n, fabs(exact-x_n), nSteps});
        x0 = x_n;
        nSteps++;
    }
}