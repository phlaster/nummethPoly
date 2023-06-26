#include "functions.hpp"

// Функция из условия
double f1(double x){
    return pow(x, 4) - pow(x, 3) -2*pow(x, 2) + 3*x - 3;
}

// Первая производная (аналитически)
double f1d1(double x){
    return 4*pow(x, 3) -3*pow(x, 2) - 4*x + 3;
}

// Вторая производная (аналитически)
double f1d2(double x){
    return 12*pow(x, 2) - 6*x - 4;
}


// Функция из условия
double f2(double x){
    return 3 * exp(x) - 5*x - 3;
}

// Первая производная (аналитически)
double f2d1(double x){
    return 3 * exp(x) - 5;
}

// Вторая производная (аналитически)
double f2d2(double x){
    return 3 * exp(x);
}

// Ункция определения знака
int sign(double x){
    return x>0.0 ?  1 :
           x<0.0 ? -1 : 0;
}
