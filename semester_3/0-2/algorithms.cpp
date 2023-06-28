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
            if (j != i && fabs(x - x_values[j]) > DBL_EPSILON) {
                term *= (x - x_values[j]) / (x_values[i] - x_values[j]);
            }
        }
        result += term;
    }
    return result;
}

Graphic lagrangeInterpol(const Graphic& func, double left, double right, int N)
{
    double dx = (right - left) / (N - 1);

    Vec xVals(N), yVals(N);
    for (int i = 0; i < N; ++i) {
        xVals[i] = left + i * dx;
        yVals[i] = lagrangeTerm(xVals[i], func.xVals, func.yVals);
    }
    return Graphic({xVals, yVals, dx, N});
}

Graphic lagrangeInterpol(const Graphic& func, double left, double right, double dx)
{
    int N = ceil((right - left) / dx);
    return lagrangeInterpol(func, left, right, N);
}



double hermitePolynomial(double a, double b, double fa, double fb, double dfa, double dfb, double t) {
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

Graphic hermiteSpline(const Graphic& function, const Graphic& derivative, double dt)
{
    Graphic hermite;
    hermite.N = 0;
    hermite.dx = dt < function.dx ? dt : function.dx;
    for (int i = 0; i < function.N - 1; ++i)
    {
        double a  = function.xVals[i];    // Начальная точка подотрезка
        double b  = function.xVals[i+1];  // Конечная точка подотрезка

        double fa = function.yVals[i];    // Значение функции в начальной точке
        double fb = function.yVals[i+1];  // Значение функции в конечной точке

        double dfa= derivative.yVals[i];  // Значение производной в начальной точке 
        double dfb= derivative.yVals[i+1];// Значение производной в конечной точке 

        double t = a;
        do {
            hermite.xVals.push_back(t);
            hermite.yVals.push_back(hermitePolynomial(a, b, fa, fb, dfa, dfb, t));
            hermite.N++;
            t += dt;
        } while (t < b);
    }
    return hermite;
}

Graphic hermiteSpline(const Graphic& function, const Graphic& derivative, int N)
{
    double dt = function.dx / N;
    return hermiteSpline(function, derivative, dt);
}