// Вариант 11 Б
/*
    f1(x) = ctg(x) + x^2
    f2(x) = x^5 - 3.2x^3 + 2.5x^2 - 7x + 1.5
*/

#include "functions.hpp"

const Vec LIMS_1 = {0.5, 2.75};
const Vec LIMS_2 = {-2.4, 2.1};

int main()
{
    double t = 0.5;
    double dot = t * sum(LIMS_1);
    int nInterpol = 1;

    for (int N=4; N<=100; N++){
        int nNodes = N * nInterpol;

        Vec uniformNodes_grid = uniformGrid(LIMS_1[0], LIMS_1[1], nNodes);
        Vec uniformTab_grid = {dot};

        Graphic uniformNodes = tabulateFunction(f1, uniformNodes_grid);
        
        // Методом секущих вычисляем значение производной в узлах рабочей решётки
        Graphic uniformNodes_DerNum = tabulateDerivativeNum(uniformNodes); // ??????????????
        
        // Табулируем значения функции и её производной функции по подробной решётке
        Graphic uniformTab = tabulateFunction(f1, uniformTab_grid);
        Graphic uniformTab_Der = tabulateDerivative(f1, uniformTab);
        
        // Вычисляем полином Лагранжа и сплайн Эрмита
        Graphic uniformLagr = lagrangeInterpol(uniformNodes, uniformTab_grid);
        Graphic uniformHerm = hermiteSpline(uniformNodes, uniformNodes_DerNum, nInterpol);

        // Вычисляем профили ошибок для интерполяций
        Vec err_uniformLagr = errorProfile(uniformLagr.yVals, uniformTab.yVals);
        Vec err_uniformHerm = errorProfile(uniformHerm.yVals, uniformTab.yVals);

        cout << nNodes << endl;;
        cout << nInterpol << endl;;
        cout << mean(err_uniformLagr) << endl;;
        cout << mean(err_uniformHerm) << endl;;

        write2_3(
                "CSVs/task2-3_4-100_singledot_f1.csv",
                uniformNodes,
                uniformNodes_DerNum,

                uniformTab,
                uniformTab_Der,

                uniformLagr,
                uniformHerm,

                err_uniformLagr,
                err_uniformHerm
            );
    }
    
    return 0;
}
