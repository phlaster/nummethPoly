#include "functions.hpp"
#include <string>

void task2_3(double (*f)(double, bool), int minNodes, int maxNodes, int coeffMult, Vec lims, Str threadname, Buffer& buffer)
{

    double left = lims[0];
    double right = lims[1];

    buffer.append("nNodes,nInterpol,err_LagrUniform,err_Herm\n");

    for (int nNodes=minNodes; nNodes<=maxNodes; nNodes++)
    {   
        // Задаём количество узлов интерполяции
        int nInterpol = coeffMult * nNodes;

        // Задаём равномерные решётки для узлов и для подробно табулированной функции
        Vec uniformNodes_grid = uniformGrid(left, right, nNodes);
        Vec uniformTab_grid = uniformGrid(left, right, nInterpol);

        // Вычисляем значение функции в узлах рабочей и подробной решёток
        Graphic uniformNodes = tabulateFunction(f, uniformNodes_grid);
       
        // Методом секущих вычисляем значение производной в узлах рабочей решётки
        Graphic uniformNodes_DerNum = tabulateDerivativeNum(uniformNodes);
        
        // Табулируем значения функции и её производной функции по подробной решётке
        Graphic uniformTab = tabulateFunction(f, uniformTab_grid);
        Graphic uniformTab_Der = tabulateDerivative(f, uniformTab);
        
        // Вычисляем полином Лагранжа и сплайн Эрмита
        Graphic uniformLagr = lagrangeInterpol(uniformNodes, nInterpol);
        Graphic uniformHerm = hermiteSpline(uniformNodes, uniformNodes_DerNum, nInterpol);

        // Вычисляем профили ошибок для интерполяций
        Vec err_uniformLagr = errorProfile(uniformLagr.yVals, uniformTab.yVals);
        Vec err_uniformHerm = errorProfile(uniformHerm.yVals, uniformTab.yVals);

        // Запись в буфер
        buffer.append(nNodes);
        buffer.append(nInterpol);
        buffer.append(mean(err_uniformLagr));
        buffer.append(mean(err_uniformHerm), true);

        write2_3(
            "CSVs/task2-3_" + threadname + "_" + to_string(nNodes) + ".csv",
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
}

void task4(double (*f)(double, bool), int minNodes, int maxNodes, int coeffMult, Vec lims, Str threadname, Buffer& buffer)
{
    double left = lims[0];
    double right = lims[1];

    buffer.append("nNodes,nInterpol,err_LagrUniform,err_LagrCheb\n");

    for (int nNodes=minNodes; nNodes<=maxNodes; nNodes++)
    {   
        // Задаём количество узлов интерполяции
        int nInterpol = coeffMult * nNodes;

        // Задаём равномерные решётки для узлов и для подробно табулированной функции
        Vec uniformNodes_grid = uniformGrid(left, right, nNodes);
        Vec uniformTab_grid = uniformGrid(left, right, nInterpol);
        // Вычисляем значение функции в узлах
        Graphic uniformNodes = tabulateFunction(f, uniformNodes_grid);
        // Табулируем значение функции по подробной решётке
        Graphic uniformTab = tabulateFunction(f, uniformTab_grid);
        // Вычисляем полином Лагранжа на подробной решётке
        Graphic uniformLagr = lagrangeInterpol(uniformNodes, nInterpol);


        // Задаем решётки Чебышева для узлов и для подробно табулированной функции
        Vec chebNodes_grid = chebyshevGrid(left, right, nNodes);
        // Вычисляем значение функции в узлах решётки Чебышева
        Graphic chebNodes = tabulateFunction(f, chebNodes_grid);
        // Вычисляем полином Лагранжа на подробной решётке Чебышева
        Graphic chebLagr = lagrangeInterpol(chebNodes, nInterpol);
        Graphic chebLagrTab = tabulateFunction(f, chebLagr.xVals);


        // Вычисляем профили ошибок для равномерных решёток и для решётки Чебышева
        Vec err_uniformLagr = errorProfile(uniformLagr.yVals, uniformTab.yVals);
        Vec err_chebLagr = errorProfile(chebLagr.yVals, chebLagrTab.yVals);


        // Запись в буфер
        buffer.append(nNodes);
        buffer.append(nInterpol);
        buffer.append(mean(err_uniformLagr));
        buffer.append(mean(err_chebLagr), true);

        write4(
            "CSVs/task4_" + threadname + "_" + to_string(nNodes) + ".csv",
            uniformNodes,
            chebNodes,

            uniformTab,
            uniformLagr,
            chebLagr,

            err_uniformLagr,
            err_chebLagr
        );
    }
}

void task5(double (*f)(double, bool), int minNodes, int maxNodes, int coeffMult, double distModulo, Vec lims, Str threadname, Buffer& buffer)
{
    double left = lims[0];
    double right = lims[1];

    buffer.append("nNodes,nInterpol,err_LagrUniform_Noize,err_LagrCheb_Noize,err_Herm_Noize\n");

    for (int nNodes=minNodes; nNodes<=maxNodes; nNodes++)
    {   
        // Задаём количество узлов интерполяции
        int nInterpol = coeffMult * nNodes;

        // Задаём решётки узлов
        Vec uniformNodes_grid = uniformGrid(left, right, nNodes);
        Vec chebNodes_grid = chebyshevGrid(left, right, nNodes);
        Vec uniformTab_grid = uniformGrid(left, right, nInterpol);
       
        // Точные значения узлов
        Graphic uniformNodes = tabulateFunction(f, uniformNodes_grid);
        Graphic chebNodes = tabulateFunction(f, chebNodes_grid);

        // Добавляем шум в узлах
        Graphic uniformNodes_Noize = deviate(uniformNodes, distModulo);
        Graphic chebNodes_Noize = deviate(chebNodes, distModulo);

        // Вычисляем методом секущих производные в узлах с шумом
        Graphic uniformNodes_Noize_DerNum = tabulateDerivativeNum(uniformNodes_Noize);

        // Интерполируем на решётках с шумом
        Graphic uniformLagr_Noize = lagrangeInterpol(uniformNodes_Noize, nInterpol);
        Graphic uniformHerm_Noize = hermiteSpline(uniformNodes_Noize, uniformNodes_Noize_DerNum, nInterpol);
        Graphic chebLagr_Noize = lagrangeInterpol(chebNodes_Noize, nInterpol);

        // Табулируем точные значения функции и производной по подробной решётке для вычисления ошибок
        Graphic uniformTab = tabulateFunction(f, uniformTab_grid);
        Graphic uniformTab_Der = tabulateDerivative(f, uniformTab);

        // Вычисляем профили ошибок для решёток с шумом
        Vec err_uniformLagr_Noize = errorProfile(uniformLagr_Noize.yVals, uniformTab.yVals);
        Vec err_uniformHerm_Noize = errorProfile(uniformHerm_Noize.yVals, uniformTab.yVals);
        Vec err_chebLagr_Noize = errorProfile(chebLagr_Noize.yVals, uniformTab.yVals);

        // Запись в буфер
        buffer.append(nNodes);
        buffer.append(nInterpol);
        buffer.append(mean(err_uniformLagr_Noize));
        buffer.append(mean(err_chebLagr_Noize));
        buffer.append(mean(err_uniformHerm_Noize), true);

        write5(
            "CSVs/task5_" + threadname + "_" + to_string(nNodes)+"_"+to_string(distModulo)+"_.csv",
            uniformNodes_Noize,
            chebNodes_Noize,
            uniformNodes_Noize_DerNum,

            uniformTab,
            uniformTab_Der,

            uniformLagr_Noize,
            uniformHerm_Noize,
            chebLagr_Noize,

            err_uniformLagr_Noize,
            err_uniformHerm_Noize,
            err_chebLagr_Noize
        );
    }
}