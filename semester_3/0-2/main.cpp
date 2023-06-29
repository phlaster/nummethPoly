// Вариант 11 Б
/*
    f1(x) = ctg(x) + x^2
    f2(x) = x^5 - 3.2x^3 + 2.5x^2 - 7x + 1.5
*/
#include "functions.hpp"
#include <ios>
#include <fstream>
#include <iomanip>

const Vec LIMS_1 = {0.5, 2.75};
const Vec LIMS_2 = {-2.4, 2.1};

int task2(double (*f)(double, bool), int nNodes, int nInterpol, Str filename)
{
    Vec lims = (f==f1) ? LIMS_1 : LIMS_2; // Подходящие пределы
    double left = lims[0];
    double right = lims[1];

    // Небольшое количество точных узлов для дальнейших расчётов
    Graphic exactNodes = tabulateFunction(f, left, right, nNodes);

    //Точная табуляция для сравнения
    Graphic exactFunction = tabulateFunction(f, left, right, nInterpol);
    
    // Точные производные для сравнения
    Graphic exactDerivative = tabulateDerivative(f, exactFunction);
    
    // Производные, взятые численно (в соответствии с заданием)
    Graphic derivativeNum = tabulateDerivativeNum(f, exactNodes);

    // У приближений больше узлов
    Graphic interpolL = lagrangeInterpol(exactNodes, left, right, nInterpol);
    Graphic interpolH = hermiteSpline(exactNodes, derivativeNum, nInterpol);

    ofstream outStream(filename);
    if (!outStream.is_open()) {
        throw std::runtime_error("Не удалось открыть файл для записи!");
    }

    outStream << "n,nodesX,nodesY,derNumY,x_exact,y_exact,dfdx_exact,LagrangeX,LagrangeY,errL,HermitX,HermitY,errH\n";
    for (int i = 0; i<nInterpol; i++){
        double abs_errL = fabs(exactFunction.yVals[i] - interpolL.yVals[i]);
        double abs_errH = fabs(exactFunction.yVals[i] - interpolH.yVals[i]);

        int decimalsL = abs_errL < 1e-2 ? ceil(-log10(abs_errL)) : 2;
        int decimalsH = abs_errH < 1e-2 ? ceil(-log10(abs_errH)) : 2;
        
        outStream << i+1 << ",";
        if (i>=nNodes) {
            outStream << ",,,";
        }
        else {
            outStream
            << exactNodes.xVals[i] << ","
            << exactNodes.yVals[i] << ","
            << derivativeNum.yVals[i] << ",";
        }
        
        outStream 
            // << fixed << setprecision(16)
            << exactFunction.xVals[i] << ","
            << exactFunction.yVals[i] << ","
            << exactDerivative.yVals[i] << ","
            // << fixed << setprecision(decimalsL)
            << interpolL.xVals[i] << ","
            << interpolL.yVals[i] << ","
            // << scientific << setprecision(1) 
            << abs_errL << ","
            // << fixed << setprecision(decimalsH)
            << interpolH.xVals[i] << ","
            << interpolH.yVals[i] << ","
            // << scientific << setprecision(1)
            << abs_errH << "\n";
    }


    outStream.close();
    return 0;
}


int main(int argc, char *argv[])
{
    
    int nNodes = stoi(argv[1]);
    int nInterpol = stoi(argv[2]);

    task2(f1, nNodes, nInterpol, "f1.csv");
    task2(f2, nNodes, nInterpol, "f2.csv");
    return 0;
}