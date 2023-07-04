#include "functions.hpp"


double lagrangeTerm(double x, const Vec& x_values, const Vec& y_values)
{
    double result = 0;
    for (size_t i = 0; i < x_values.size(); ++i) {
        double term = y_values[i];
        
        if(fabs(x - x_values[i]) < DBL_EPSILON) {
            return term;
        }

        for (size_t j = 0; j < x_values.size(); ++j) {
            if (j != i && fabs(x_values[i] - x_values[j]) > DBL_EPSILON) {
                term *= (x - x_values[j]) / (x_values[i] - x_values[j]);
            }
        }
        result += term;
    }
    return result;
}

// Graphic lagrangeInterpol(const Graphic& g, double left, double right, int N)
// {
//     double dx = (right - left) / (N - 1);

//     Vec xVals(N), yVals(N);
//     for (int i = 0; i < N; ++i) {
//         xVals[i] = left + i * dx;
//         yVals[i] = lagrangeTerm(xVals[i], g.xVals, g.yVals);
//     }
//     return Graphic({xVals, yVals, dx, N});
// }

Graphic lagrangeInterpol(const Graphic& g, int N)
{
    double left = g.xVals[0];
    double right = g.xVals[g.N-1];
    double dx = (right - left) / (N-1);
    Graphic lagrange({Vec(N), Vec(N), dx, N});
    for (int i = 0; i < N; ++i) {
        lagrange.xVals[i] = left + i * dx;
        lagrange.yVals[i] = lagrangeTerm(lagrange.xVals[i], g.xVals, g.yVals);
    }
    return lagrange;
}


Graphic lagrangeInterpol(const Graphic& g, double left, double right, double dx)
{
    int N = ceil((right - left) / dx);
    return lagrangeInterpol(g, left, right, N);
}



double hermiteTerm(double a, double b, double fa, double fb, double dfa, double dfb, double t) {
    // Вычисляем нормализованное значение t на отрезке [0, 1]
    double tNorm = (t - a) / (b - a);
    
    // Вычисляем значения коэффициентов полинома Эрмита
    double h00 = (2 * pow(tNorm , 3)) - (3 * pow(tNorm , 2)) + 1;
    double h10 = pow(tNorm , 3) - (2 * pow(tNorm , 2)) + tNorm ;
    double h01 = (-2 * pow(tNorm , 3)) + (3 * pow(tNorm , 2));
    double h11 = pow(tNorm , 3) - pow(tNorm , 2);
    
    // Вычисляем значение полинома Эрмита H(t)
    return fa*h00 + dfa*(b-a)*h10 + fb*h01 + dfb*(b-a)*h11;
}

Graphic hermiteSpline(const Graphic& function, const Graphic& derivative, int N)
{
    double start = function.xVals[0];
    double finish = function.xVals[function.N-1];
    double dt = (finish - start) / (N-1);
    Graphic hermite({Vec(N), Vec(N), dt, N});
    double t = start;
    int i_H = 0;
    int i_F = 0;
    while (t<=finish || i_H<N){
        double a  = function.xVals[i_F];
        double b  = function.xVals[i_F+1];

        double fa = function.yVals[i_F];
        double fb = function.yVals[i_F+1];

        double dfa= derivative.yVals[i_F];
        double dfb= derivative.yVals[i_F+1];

        hermite.xVals[i_H] = t;
        hermite.yVals[i_H] = hermiteTerm(a, b, fa, fb, dfa, dfb, t);

        t += dt;
        i_H++;
        if (t>b)
            i_F++;
    }
    return hermite;
}


Graphic hermiteSpline(const Graphic& function, const Graphic& derivative, double dt)
{
    int N = ceil((function.xVals[function.N-1] - function.xVals[0]) / dt);
    return hermiteSpline(function, derivative, N);
}