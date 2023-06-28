#include "functions.hpp"

// Табуляция с известным количеством точек
Graphic tabulateFunction(double (*f)(double, bool),
                        double left,
                        double right,
                        int N)
{
    double dx = (right - left) / (N - 1);
    Vec x_grid(N), y_prescise(N);
    for (int i = 0; i < N; ++i) {
        x_grid[i] = left + i * dx;
        y_prescise[i] = f(x_grid[i], false);
    }
    return Graphic({x_grid, y_prescise, dx, N});
}

// Табуляция с выбранным шагом
Graphic tabulateFunction(double (*f)(double, bool),
                        double left,
                        double right,
                        double dx)
{
    int N = ceil((right - left) / dx);
    return tabulateFunction(f, left, right, N);
}


// Вычисление точек производной
Graphic tabulateDerivative(double (*f)(double, bool), const Graphic& nodes)
{
    Graphic derivative({nodes.xVals, Vec(), nodes.dx, nodes.N});
    for (auto x_i : derivative.xVals) 
    {
        derivative.yVals.push_back(f(x_i, true));
    }
    return derivative;
}

// Внутренние производные
Graphic truncatedTangents(const Graphic& gr)
{
    int n_tg = gr.N-2;
    Vec tangents(n_tg), x_grid;
    double dx = gr.dx;
    for (int i=1; i<=n_tg; i++)
    {
        tangents[i-1] = (gr.yVals[i+1]-gr.yVals[i-1])/(2*dx);
        x_grid.push_back(gr.xVals[i]);
    }
    return Graphic({x_grid, tangents, dx, n_tg});
}

// Добавляем внешние (по Лагранжу)
Graphic merge3(double left, const Graphic& main, double right){
    Vec xVals = {main.xVals[0] - main.dx};
    xVals.insert(xVals.end(), main.xVals.begin(), main.xVals.end());
    xVals.push_back(main.xVals[main.N-1] + main.dx);

    Vec yVals = {left};
    yVals.insert(yVals.end(), main.yVals.begin(), main.yVals.end());
    yVals.push_back(right);

    int N = main.N + 2;
    double dx = main.dx;

    return Graphic({xVals, yVals, dx, N});
}


Graphic tabulateDerivativeNum(double (*f)(double, bool), const Graphic& main, int LagrangePoints)
{
    double left = main.xVals[0];
    double right = main.xVals[main.N-1];
    double dx = main.dx;

    Graphic leftEnd = tabulateFunction(f, left-(LagrangePoints-1)*dx, left, LagrangePoints);
    double beyondLeft = lagrangeTerm(left-dx, leftEnd.xVals, leftEnd.yVals);

    Graphic rightEnd = tabulateFunction(f, right, right+(LagrangePoints-1)*dx, LagrangePoints);
    double beyondRight = lagrangeTerm(right+dx, rightEnd.xVals, rightEnd.yVals);

    Graphic extended = merge3(beyondLeft, main, beyondRight);

    return truncatedTangents(extended);
}
