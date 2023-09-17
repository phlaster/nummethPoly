#include "functions.hpp"

// Равномерная решётка
Vec uniformGrid(double a, double b, int n)
{
    Vec grid(n);
    double dx = (b - a) / (n - 1);
    for(int i = 0; i < n; i++)
        grid[i] = a + i * dx;
    return grid;
}

// Решётка Чебышева 
Vec chebyshevGrid(double a, double b, int n)
{
    Vec grid(n);
    for(int k = 0; k < n; k++) {
        double t = cos(M_PI * (k + 0.5) / n); 
        grid[n-k-1] = (a+b)/2 + t*(b-a)/2;
    }
    return grid;
}


// Табуляция на заданной решётке
Graphic tabulateFunction(double (*f)(double, bool), const Vec& grid)
{
    Graphic g;
    g.N = grid.size();
    g.xVals = grid;
    g.yVals = Vec(g.N);
    g.dx = fabs(grid[grid.size()-1] - grid[0])/(g.N-1);
    for (int i=0; i<g.N; i++)
        g.yVals[i] = f(grid[i], false);
    return g;
}

// Вычисление точек производной
Graphic tabulateDerivative(double (*f)(double, bool), const Graphic& g)
{
    Graphic derivative({g.xVals, Vec(), g.dx, g.N});
    for (auto x_i : derivative.xVals) 
        derivative.yVals.push_back(f(x_i, true));
    return derivative;
}

// Внутренние производные
Graphic truncatedTangents(const Graphic& g)
{
    Graphic t({Vec(g.N-2), Vec(g.N-2), g.dx, g.N-2});
    for (int i=0; i<t.N; i++)
    {
        t.xVals[i] = g.xVals[i+1];
        t.yVals[i] = (g.yVals[i+2]-g.yVals[i]) / (g.xVals[i+2]-g.xVals[i]);
    }
    return t;
}

// Добавляем внешние (по Лагранжу)
Graphic merge3(double left, const Graphic& g, double right){
    Graphic merged({Vec(g.N+1), Vec(g.N+1), g.dx, g.N+2});

    merged.xVals[0] = (g.xVals[0] - g.dx);
    merged.yVals[0] = (left);

    for (int i=0; i<g.N; i++)
    {
        merged.xVals[i+1] = g.xVals[i];
        merged.yVals[i+1] = g.yVals[i];
    }

    merged.xVals.push_back(g.xVals[g.N-1] + g.dx);
    merged.yVals.push_back(right);

    return merged;
}

// Взятие производных методом секущих
Graphic tabulateDerivativeNum(const Graphic& g, int LagrangePoints)
{
    double left = g.xVals[0];
    double right = g.xVals[g.N-1];
    double dx = g.dx;

    double beyondLeft = lagrangeTerm(
        left-dx,
        Vec(g.xVals.begin(), g.xVals.begin()+LagrangePoints),
        Vec(g.yVals.begin(), g.yVals.begin()+LagrangePoints)
    );

    double beyondRight = lagrangeTerm(
        right+dx,
        Vec(g.xVals.end()-LagrangePoints, g.xVals.end()),
        Vec(g.yVals.end()-LagrangePoints, g.yVals.end())
    );

    Graphic extended = merge3(beyondLeft, g, beyondRight);

    return truncatedTangents(extended);
}

// Внесение шума в график
Graphic deviate(const Graphic& g, double sigma)
{
    random_device randomSource;
    mt19937 rGen(randomSource());
    normal_distribution<double> normDist(0, sigma);
    Graphic dev({g.xVals, Vec(g.N), g.dx, g.N});
    for (int i=0; i< g.N; i++)
        dev.yVals[i] = g.yVals[i] + normDist(rGen);
    return dev;
}