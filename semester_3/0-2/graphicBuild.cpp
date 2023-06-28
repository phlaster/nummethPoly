#include "functions.hpp"

Graphic calculateGraphic(double (*f)(double, bool),
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

// Перегрузка
Graphic calculateGraphic(double (*f)(double, bool),
                        double left,
                        double right,
                        double dx)
{
    int N = ceil((right - left) / dx);
    return calculateGraphic(f, left, right, N);
}


// Вычисление точек производной
Graphic calculateDerivativeAnalitical(double (*f)(double, bool), const Vec& x_grid)
{
    Vec derivative(x_grid.size());
    double dx = x_grid[2]-x_grid[1];
    int N = x_grid.size();
    for (size_t i = 0; i < N; ++i) {
        derivative[i] = f(x_grid[i], true);
    }
    return Graphic({x_grid, derivative, dx, N});
}


Graphic tangents_trunc(const Graphic& gr)
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


Graphic calculateDerivativeNumerical(double (*f)(double, bool), const Graphic& main, int LagrangeDegree)
{
    double left = main.xVals[0];
    double right = main.xVals[main.N-1];
    double dx = main.dx;

    Graphic leftEnd = calculateGraphic(f, left-(LagrangeDegree-1)*dx, left, LagrangeDegree);
    double beyondLeft = lagrange(left-dx, leftEnd.xVals, leftEnd.yVals);

    Graphic rightEnd = calculateGraphic(f, right, right+(LagrangeDegree-1)*dx, LagrangeDegree);
    double beyondRight = lagrange(right+dx, rightEnd.xVals, rightEnd.yVals);


    return tangents_trunc(merge3(beyondLeft, main, beyondRight));
}
