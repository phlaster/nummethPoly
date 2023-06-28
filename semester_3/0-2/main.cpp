// Вариант 11 Б
/*
    f1(x) = ctg(x) + x^2
    f2(x) = x^5 - 3.2x^3 + 2.5x^2 - 7x + 1.5
*/
#include "functions.hpp"
// #include <fstream>
// #include <iomanip>

const Vec LIMS_1 = {0.5, 2.9};
const Vec LIMS_2 = {-2.5, 2.3};



int main()
{
    int N = 3; // Степень полинома Лагранжа
    double dx = 0.5;
    double (*f)(double, bool) = f1;
    Vec lims = (f==f1) ? LIMS_1 : LIMS_2;

    // С использованием полинома Лагранжа
    Graphic function = calculateGraphic(f, lims[0], lims[1], dx);
    Graphic derivative = calculateDerivativeNumerical(f, function, N);

    // cout.precision(3);
    cout << function.N << " значений функции.\n";
    cout << derivative.N << " значений её производной.\n";
    cout << " i    x       f(x)    df/dx    dL/dx\n";
    
    for (int i = 0; i<5; i++){
        cout << i+1 << "   " << fixed
            << function.xVals[i] << "   "
            << function.yVals[i] << "   "
            << f(function.xVals[i], true) << "   "
            << derivative.yVals[i] << "   "
            << endl;
    }
    cout << "<...>\n";
    for (int i= function.N -5; i<function.N ; i++){
        cout << i+1 << "   " << fixed
            << function.xVals[i] << "   "
            << function.yVals[i] << "   "
            << f(function.xVals[i], true) << "   "
            << derivative.yVals[i] << "   "
            << endl;
    }

    Graphic H = hermiteSpline(function, derivative, 0.1);
    cout << "Эрмитов сплайн:\nx:                 y:\n";
    for (int i=0; i<H.N; i++)
    {
        cout << i<< " " << H.xVals[i] << " " << H.yVals[i] << endl;
    }
    return 0;
}
