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

// Вычисление точек графика функции
Coords calculateGraphic(double (*f)(double, bool), const Vec& xLims, const int N)
{
    // Шаг дискретизации
    double dx = (xLims[1] - xLims[0]) / (N - 1);
    
    // Массивы координат точек
    Vec x_grid(N), y_prescise(N);

    // Заполнение
    for (int i = 0; i < N; ++i) {
        x_grid[i] = xLims[0] + i * dx;
        y_prescise[i] = f(x_grid[i], false);
    }
    return Coords(x_grid, y_prescise);
}


double lagrange(double x, const Vec& x_values, const Vec& y_values)
{
    double result = 0;
    for (size_t i = 0; i < x_values.size(); ++i) {
        double term = y_values[i];
        
        if(fabs(x - x_values[i]) < DBL_EPSILON) { // защита от деления на ноль
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