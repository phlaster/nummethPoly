#include "functions.hpp"

ans bisection(double (*f)(double), double exact, Vec lims, double eps){
    double a = lims[0], b = lims[1], c;
    int nSteps = 1; 
    for(;;){
        c = (a+b)/2;
        if (fabs(a-c) <= eps)
            return ans({c, fabs(exact-c), nSteps});
        if (sign(f(a)) != sign(f(c)))
            b = c;
        else
            a = c;
        nSteps++;
    }
}

ans newton(double (*f)(double), double exact, Vec lims, double eps){
    // Анонимные функции для подбора подходящих аналитических
    // производных соответствующего порядка на ходу
    auto d1 = [](double (*f)(double)){ return (f==f1) ? f1d1 : f2d1; };
    auto d2 = [](double (*f)(double)){ return (f==f1) ? f1d2 : f2d2; };

    // Выбор начального приближения из условия сходимости метода
    // (равенство значения и 2-й производной в точке)
    double a = lims[0], b = lims[1], x_n, delta, x0 = (f(a)*d2(f)(a) > 0) ? a : b;
    int nSteps = 1;
    for (;;){
        x_n = x0 - f(x0)/d1(f)(x0);
        delta = fabs(x_n-x0);
        if (delta <= eps)
            return ans({x_n, fabs(exact-x_n), nSteps});
        x0 = x_n;
        nSteps++;
    }
}

ans newton_broken(double exact, Vec lims, double x0, double eps){
    double a = lims[0], b = lims[1], x_n, delta;
    int nSteps = 1;
    for (;;){
        x_n = x0 - f1(x0)/f1d1(x0);
        delta = fabs(x_n-x0);
        if (delta <= eps)
            return ans({x_n, fabs(exact-x_n), nSteps});
        x0 = x_n;
        nSteps++;
    }
}