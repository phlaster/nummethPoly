#include "functions.hpp"

double cot(double x){
	return cos(x)/sin(x);
}

double f1(double x, bool derivative){
    return derivative ? 2*x - 1/pow(sin(x), 2) :
                        cot(x) + x*x;
}

double f2(double x, bool derivative){
    return derivative ? 5*pow(x, 4) - 9.6*x*x + 5*x - 7 : 
                        pow(x, 5) - 3.2*pow(x, 3) + 2.5*x*x - 7*x + 1.5;
}

// Точные значения
double F1(double a, double b){
    auto f = [](double x){
        return log(fabs(sin(x))) + pow(x, 3)/3;
    };
    return f(b) - f(a);
}
double F2(double a, double b){
    auto f = [](double x){
        return pow(x,6)/6 - 0.8*pow(x,4) + 5*x*x*x/6 - 3.5*x*x + 1.5*x;
    };
    return f(b) - f(a);
}
