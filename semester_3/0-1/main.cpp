// Вариант 8 Б
#include "functions.hpp"

const Vec F1_INTERVAL = {-2.732, -1.693},  F2_INTERVAL = {0.6, 1.5};
const double X1 = -1.732050807568877, X2 = 0.947404881615152, X3 = 1.732050807568877;

int main(){
    errorConverg(X1, F1_INTERVAL, bisection, f1, "bisect1.csv");
    errorConverg(X2, F2_INTERVAL, bisection, f2, "bisect2.csv");

    errorConverg(X1, F1_INTERVAL, newton, f1, "newton1.csv");
    errorConverg(X2, F2_INTERVAL, newton, f2, "newton2.csv");

    brokenConverg(X3, 0.5, {0.5,4}, "broken_newton(limits).csv");
    brokenConverg(X3, 1.1, {0.5,4}, "broken_newton(limits+x0).csv");

    return 0;
}
