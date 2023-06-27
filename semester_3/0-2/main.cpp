// Вариант 11 Б
/*
    f1(x) = ctg(x) + x^2
    f2(x) = x^5 - 3.2x^3 + 2.5x^2 - 7x + 1.5
*/
#include "functions.hpp"
// #include <fstream>
#include <iomanip>

const Vec F1_INTERVAL = {0.5, 3.0};
const Vec F2_INTERVAL = {-2.5, 2.3};


int main()
{
    int N = 6; // число узлов интерполяции
    double (*f)(double, bool) = f1;

    Vec interval = (f == f1) ? F1_INTERVAL : F2_INTERVAL;

    auto [x_grid, y_prescise] = calculateGraphic(f, interval, N);

    double x = 1;
    while (x <= 2)
    {
        double y = lagrange(x, x_grid, y_prescise);
        cout << fixed << setprecision(1) << x << "  "
             << scientific << setprecision(2) << fabs(f(x, false) - y)
             << endl;
        x += 0.1;
    }
    
    return 0;
}