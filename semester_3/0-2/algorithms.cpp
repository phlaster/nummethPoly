#include "functions.hpp"


// Лагранжев член
double lagrangeTerm(double x, const Vec& x_values, const Vec& y_values)
{
    double result = 0;
    for (size_t i = 0; i < x_values.size(); ++i) {
        double term = y_values[i];
        if(fabs(x - x_values[i]) < DBL_EPSILON)
            return term;
        for (size_t j = 0; j < x_values.size(); ++j)
            if (j != i && fabs(x_values[i] - x_values[j]) > DBL_EPSILON)
                term *= (x - x_values[j]) / (x_values[i] - x_values[j]);
        result += term;
    }
    return result;
}

Graphic lagrangeInterpol(const Graphic& nodes, const Vec& grid)
{
    int nInterpol = grid.size();
    double left = grid[0];
    double right = grid[nInterpol-1];
    double dx = (right - left) / (nInterpol-1);
    Graphic lagrange({grid, Vec(nInterpol), dx, nInterpol});
    for (int i = 0; i < nInterpol; ++i) {
        lagrange.yVals[i] = lagrangeTerm(grid[i], nodes.xVals, nodes.yVals);
    }
    return lagrange;
}

pair<double, double> lagrange_uniform_single_value_with_error(
    double (*f)(double, bool),
    const double x,
    const Vec& lims,
    const int nNodes
){
    double left = lims[0];
    double right = lims[1];
    if (x< left || right < x)
        throw runtime_error("x не лежит внутри интервала интерполяции!\n");
    double y_exact = f(x, false);
    Vec uniformNodes_grid = uniformGrid(left, right, nNodes);
    Graphic uniformNodes = tabulateFunction(f, uniformNodes_grid);
    double y_appr = lagrangeTerm(x, uniformNodes.xVals, uniformNodes.yVals);
    double err = fabs(y_exact - y_appr);
    return make_pair(y_appr, err);
}

double hermiteTerm(double a, double b, double fa, double fb, double dfa, double dfb, double t)
{
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
    int i_F = 0;
    for (int i_H=0; i_H<N; i_H++){
        double a  = function.xVals[i_F];
        double b  = function.xVals[i_F+1];

        double fa = function.yVals[i_F];
        double fb = function.yVals[i_F+1];

        double dfa= derivative.yVals[i_F];
        double dfb= derivative.yVals[i_F+1];

        hermite.xVals[i_H] = t;
        hermite.yVals[i_H] = hermiteTerm(a, b, fa, fb, dfa, dfb, t);

        t = start + dt * (i_H+1);
        if (t>b)
            i_F++;
    }
    return hermite;
}
