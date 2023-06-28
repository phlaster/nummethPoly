// Вариант 11 Б
/*
    f1(x) = ctg(x) + x^2
    f2(x) = x^5 - 3.2x^3 + 2.5x^2 - 7x + 1.5
*/
#include "functions.hpp"
#include <ios>
#include <fstream>
#include <iomanip>

const Vec LIMS_1 = {0.5, 2.9};
const Vec LIMS_2 = {-2.5, 2.3};

int task2(double (*f)(double, bool), int N, Str filename)
{
    Vec lims = (f==f1) ? LIMS_1 : LIMS_2; // Подходящие пределы

    // Небольшое количество точных узлов для дальнейших расчётов
    Graphic exactNodes = tabulateFunction(f, lims[0], lims[1], N);

    //Точная табуляция для сравнения
    Graphic exactFunction = tabulateFunction(f, lims[0], lims[1], 10*N);
    
    // Точные производные для сравнения
    Graphic exactDerivative = tabulateDerivative(f, exactFunction);
    
    // Производные, взятые численно (в соответствии с заданием)
    Graphic derivativeNum = tabulateDerivativeNum(f, exactNodes);

    // У приближений в 10 раз больше узлов
    Graphic interpolL = lagrangeInterpol(exactNodes, lims[0], lims[1], 10*N);
    Graphic interpolH = hermiteSpline(exactNodes, derivativeNum, 10);

    ofstream outStream(filename);
    if (!outStream.is_open()) {
        throw std::runtime_error("Не удалось открыть файл для записи!");
    }
    outStream << "n,x_exact,y_exact,dfdx_exact,Lagrange,errL,Hermit,errH\n";
    for (int i = 0; i<N*10; i++){
        double abs_errL = fabs(exactFunction.yVals[i] - interpolL.yVals[i]);
        double abs_errH = fabs(exactFunction.yVals[i] - interpolH.yVals[i]);

        int decimalsL = abs_errL < 1e-2 ? ceil(-log10(abs_errL)) : 2;
        int decimalsH = abs_errH < 1e-2 ? ceil(-log10(abs_errH)) : 2;
        
        outStream << i+1 << "," << fixed << setprecision(16)
            << exactFunction.xVals[i] << ","
            << exactFunction.yVals[i] << ","
            << exactDerivative.yVals[i] << ","
            << fixed << setprecision(decimalsL)
            << interpolL.yVals[i] << ","
            << scientific << setprecision(1) << abs_errL << ","
            << fixed << setprecision(decimalsH)
            << interpolH.yVals[i] << ","
            << scientific << setprecision(1) << abs_errH << "\n";
    }
    outStream.close();
    return 0;
}


int main()
{
    task2(f1, 6, "f1.csv");
    task2(f2, 6, "f2.csv");
    

    

    
    // for (int i = 0; i<5; i++){
    //     cout << i+1 << "   " << fixed
    //         << function.xVals[i] << "   "
    //         << function.yVals[i] << "   "
    //         << f(function.xVals[i], true) << "   "
    //         << derivative.yVals[i] << "   "
    //         << endl;
    // }
    // cout << "<...>\n";
    // for (int i= function.N -5; i<function.N ; i++){
    //     cout << i+1 << "   " << fixed
    //         << function.xVals[i] << "   "
    //         << function.yVals[i] << "   "
    //         << f(function.xVals[i], true) << "   "
    //         << derivative.yVals[i] << "   "
    //         << endl;
    // }

    // Graphic H = hermiteSpline(function, derivative, 0.1);
    // cout << "Эрмитов сплайн:\nx:                 y:\n";
    // for (int i=0; i<H.N; i++)
    // {
    //     cout << i<< " " << H.xVals[i] << " " << H.yVals[i] << endl;
    // }
    return 0;
}
