/*
    Runge-Kutta method
    Variant 11, modified Euler method (middle-point)
    
    xy^2y` = x^2 + y^3  =>  y` = x/y^2 + y/x

    x in [1.1, 3]
    a = 1.1; b = 3

    y^3 = 3x^2(x-1)

*/

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;


double dy_dx(double x, double y){
    return x/(y * y) + y/x;
}

double y(double x){
    return cbrt(3 * x * x * (x-1));
}

double y_wave_next(double y_i, double h_i, double f){
    return y_i + h_i/2 * f;
}

double y_next(double y_i, double h_i, double f){
    return y_i + h_i * f;
}

void calculate(double a, double b, double h_i)
{
    double y_0 = y(a);
    double h = b-a;
    double x_i = a;
    int n_steps = round(h/h_i);

    cout << "x_i,y_appr_i,y_i,err\n";
    for (int i = 0; i < n_steps; i++)
    {
        double y_wave = y_wave_next(y_0, h_i, dy_dx(x_i, y_0));
        double y_i = y_next(y_0, h_i, dy_dx(x_i + h_i/2, y_wave));
        y_0 = y_i;
        x_i = a + i * h_i;

        cout << fixed << setprecision(6);
        cout << x_i << ',' << y_i << ',' << y(x_i) << "," << abs (y_i - y(x_i)) << '\n';

    }
}


void calculate1_8(double a, double b)
{
    double h = b-a;
    double x_i = a;
    cout << "h_i,err\n";
    for (double p=-1; p>=-8; p--){
        double y_0 = y(a);
        double h_i = pow(10, p);
        int n_steps = round(h/h_i);
        double max_err=0;
        for (int i = 0; i < n_steps; i++){
            double y_wave = y_wave_next(y_0, h_i, dy_dx(x_i, y_0));
            double y_i = y_next(y_0, h_i, dy_dx(x_i + h_i/2, y_wave));
            y_0 = y_i;
            x_i = a + i * h_i;
            double err = fabs(y_i - y(x_i));
            max_err = err > max_err ? err : max_err;
        }
        cout << fixed << setprecision(2-p);
        cout << h_i << ',' << max_err << "\n";
    }
}

int main()
{
    double a = 1.1, b = 3.0;
    
    // cout << "h_i=0.1\n";
    // calculate(a, b, 0.1);

    // cout << "\n\nh_i=0.05\n";
    // calculate(a, b, 0.05);
    calculate1_8(a, b);

    return 0;
}
