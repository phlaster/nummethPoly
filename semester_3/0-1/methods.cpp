#include "functions.hpp"

ans bisection(double (*f)(double), double exact, Vec lims, double eps){
    double a = lims[0];
    double b = lims[1];
    double c;

    int nSteps = 0; 
    for(;;)
    {
        c = (a+b)/2;
        if (fabs(a-c) <= eps)
            return ans({c, fabs(a-c), nSteps});
        if (sign(f(a)) != sign(f(c)))
            b = c;
        else
            a = c;
        nSteps++;
    }
}

ans newton(double (*f)(double), double exact, Vec lims, double eps){
    double a = lims[0];
    double b = lims[1];

    // Анонимные функции для подбора подходящих аналитических
    // производных соответствующего порядка на ходу
    auto fd1 = [](double (*f)(double)){
        if (f==f1)
            return f1d1;
        else
            return f2d1;
    };

    auto fd2 = [](double (*f)(double)){
        if (f==f1)
            return f1d2;
        else
            return f2d2;
    };

    // Выбор начального приближения из условия сходимости метода
    // (равенство значения и 2-й производной в точке)
    double x0 = sign(f(a))==sign(fd2(f)(a)) ? a : b;
    int nSteps = 0; 
    double x_n;

    for (;;)
    {
        x_n = x0 - f(x0)/fd1(f)(x0);
        if (fabs(x_n-x0) <= eps)
            return ans({x0, fabs(x_n-x0), nSteps});
        x0 = x_n;
        nSteps++;
    }

    return ans({x0, abs(exact-x0),nSteps});
}
