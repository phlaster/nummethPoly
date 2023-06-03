/*
    Runge-Kutta method
    Variant 11, modified Euler method (middle-point)
    
    xy^2y` = x^2 + y^3  =>  y` = x/y^2 + y/x

    x in [1.1, 3]
    a = 1.1; b = 3

    y^3 = 3x^2(x-1)

*/

#include <iostream>
#include <cmath>

using namespace std;


double dif_y(double x, double y)
{
    return x/(y * y) + y/x;
}

double y_exact(double x)
{
    return cbrt(3 * x * x * (x-1));
}

double y_wave_next(double y_i, double h_i, double f)
{
    return y_i + h_i/2 * f;
}

double y_next(double y_i, double h_i, double f)
{
    return y_i + h_i * f;
}

void calculate(double a, double b, double h_i)

{
    double y_0 = y_exact(a);
    double h = b-a;
    double y_wav, y_i;
    double x_i = a;
    int n_steps = round(h/h_i);
    for (int i = 1; i <= n_steps; i++)
    {
        y_wav = y_wave_next(y_0, h_i, dif_y(y_0, y_0+h_i));
        y_i = y_next(y_0, h_i, dif_y(x_i + h_i/2, y_wav));
        y_0 = y_wav;
        y_0 = y_i;
        x_i = a + i * h_i;
        cout << h_i << ","<< i << ',' << x_i << ',' << y_i << ',' << abs (y_i - y_exact(x_i)) << '\n';

    }
}


int main(void)
{
    double a = 1.1, b = 3.0;
    double h_i = 1 * pow(10, -1); // Здесь меняем шаг
    
    calculate(a, b, h_i);

    h_i = 5 * pow(10, -2); // Здесь меняем шаг
    calculate(a, b, h_i);
    return 0;
}
