// Вариант 11 Б
/*
    f1(x) = ctg(x) + x^2
    f2(x) = x^5 - 3.2x^3 + 2.5x^2 - 7x + 1.5
*/
#include "functions.hpp"
#include <string>
#include <thread>

const Vec LIMS_1 = {0.5, 2.75};
const Vec LIMS_2 = {-2.4, 2.1};



// void task2(double (*f)(double, bool), int minNodes, int maxNodes, bool toWrite)
// {
//     Vec lims;
//     Str threadname;
//     if (f==f1){
//         lims = LIMS_1;
//         threadname = "f1";
        
//     } else {
//         lims = LIMS_2;
//         threadname = "f2";
//     }
//     double left = lims[0];
//     double right = lims[1];

//     cout
//         << "Задание 2: Средняя ошибка интерполяции "
//         << threadname
//         << " полиномом Лагранжа и сплайном Эрмита:\n";
//     ;

//     if (!toWrite)
//     {
//         cout << "Узлов   точек     Лагранж       Эрмит\n";
//     }

//     for (int nNodes=minNodes; nNodes<=maxNodes; nNodes++)
//     {
//         int nInterpol = 100*nNodes;
//         if (!toWrite)
//             cout << nNodes << "       " << nInterpol << "       ";
//         ;

//         // Небольшое количество точных узлов для дальнейших расчётов
//         Graphic exactNodes = tabulateFunction(f, left, right, nNodes);

//         //Точная табуляция для сравнения
//         Graphic exactFunction = tabulateFunction(f, left, right, nInterpol);
        
//         // Точные производные для сравнения
//         Graphic exactDerivative = tabulateDerivative(f, exactFunction);
        
//         // Производные, взятые численно (в соответствии с заданием)
//         Graphic derivativeNum = tabulateDerivativeNum(f, exactNodes);

//         // У приближений больше узлов
//         Graphic interpolL = lagrangeInterpol(exactNodes, nInterpol);
//         Graphic interpolH = hermiteSpline(exactNodes, derivativeNum, nInterpol);

//         if (!toWrite){
//             cout
//                 << error(interpolL.yVals, exactFunction.yVals) << "       "
//                 << error(interpolH.yVals, exactFunction.yVals) << "\n";
//         }else{
//             task2Print(
//                 exactNodes,
//                 exactFunction,
//                 exactDerivative,
//                 derivativeNum,
//                 interpolL,
//                 interpolH,
//                 "CSVs/" + threadname+"_" + to_string(nNodes) + ".csv"
//             );
//         }
//     }
// }

// void task3(double (*f)(double, bool), int minNodes, int maxNodes, bool toWrite)
// {
//     Vec lims;
//     Str threadname;
//     if (f==f1){
//         lims = LIMS_1;
//         threadname = "f1";
        
//     } else {
//         lims = LIMS_2;
//         threadname = "f2";
//     }

//     double left = lims[0];
//     double right = lims[1];

//     cout
//         << "Задание 3: Для "
//         << threadname
//         << " сравнение сходимости полинома Лагранжа на равномерной и чебышевской сетках:\n";
//     ;

//     if (!toWrite)
//     {
//         cout << "Узлов   точек      равномерная        Чебышев\n";
//     }

//     for (int nNodes=minNodes; nNodes<=maxNodes; nNodes++)
//     {
//         int nInterpol = nNodes*100;
//         if (!toWrite)
//             cout << nNodes << "       " << nInterpol << "       ";
//         ;

//         Graphic uniformNodes = tabulateFunction(f, left, right, nNodes);
//         Graphic evenlyL = lagrangeInterpol(uniformNodes, nInterpol);
//         Graphic exactForEvenly = tabulateFunction(f, evenlyL.xVals);


//         Vec chebyshev = chebyshevGrid(left, right, nNodes);
//         Graphic chebyshevNodes = tabulateFunction(f, chebyshev);
//         Graphic chebyshevL = lagrangeInterpol(chebyshevNodes, nInterpol);
//         Graphic exactForChebyshev = tabulateFunction(f, chebyshevL.xVals);
        
//         if (!toWrite){
//             cout
//                 << error(evenlyL.yVals, exactForEvenly.yVals) << "    "
//                 << error(chebyshevL.yVals, exactForChebyshev.yVals) << "\n";
//         }else{
//             task3Print(
//                 uniformNodes,
//                 chebyshevNodes,
//                 exactForEvenly,
//                 exactForChebyshev,
//                 evenlyL,
//                 chebyshevL,
//                 nNodes,
//                 nInterpol,
//                 "CSVs/" + threadname + "_Cheb_" + to_string(nNodes) + ".csv"
//             );
//         }
//     }
// }

void task(double (*f)(double, bool), int minNodes, int maxNodes, int coeffMult, bool toWrite)
{
    Vec lims;
    Str threadname;
    if (f==f1){
        lims = LIMS_1;
        threadname = "f1";
        
    } else {
        lims = LIMS_2;
        threadname = "f2";
    }
    double left = lims[0];
    double right = lims[1];

    if (!toWrite)
        cout 
            << threadname+":\n"
            << "nNodes,nInterpol,errLagrUniform,errLagrCheb,errHerm\n";
    ;

    for (int nNodes=minNodes; nNodes<=maxNodes; nNodes++)
    {
        int nInterpol = coeffMult * nNodes;

        Graphic uniformNodes = tabulateFunction(f, left, right, nNodes);
        Graphic uniformNodes_DerNum = tabulateDerivativeNum(f, uniformNodes);
        Graphic uniformTab = tabulateFunction(f, left, right, nInterpol);
        Graphic uniformTab_Der = tabulateDerivative(f, uniformTab);
        Graphic uniformLagr = lagrangeInterpol(uniformNodes, nInterpol);
        Graphic uniformHerm = hermiteSpline(uniformNodes, uniformNodes_DerNum, nInterpol);

        Graphic chebNodes = tabulateFunction(f, chebyshevGrid(left, right, nNodes));
        Graphic chebLagr = lagrangeInterpol(chebNodes, nInterpol);
        Graphic chebTab = tabulateFunction(f, chebLagr.xVals);
        
        Vec err_uniformLagr = errorProfile(uniformLagr.yVals, uniformTab.yVals);
        Vec err_uniformHerm = errorProfile(uniformHerm.yVals, uniformTab.yVals);
        Vec err_chebLagr = errorProfile(chebLagr.yVals, chebTab.yVals);

        if (!toWrite)
            cout
                << nNodes << "," << nInterpol << ","
                << mean(err_uniformLagr) << ","
                << mean(err_uniformHerm) << ","
                << mean(err_chebLagr) << "\n";
        else
            print(
                "CSVs/" + threadname + "_" + to_string(nNodes) + ".csv",
                uniformNodes,
                uniformNodes_DerNum,
                chebNodes,
                uniformTab,
                uniformTab_Der,
                chebTab,
                uniformLagr,
                uniformHerm,
                chebLagr,
                err_uniformLagr,
                err_uniformHerm,
                err_chebLagr
            );
    }
    
}

int main(int argc, char *argv[])
{
    int minNodes = stoi(argv[1]);
    int maxNodes = stoi(argv[2]);
    assert(minNodes < maxNodes && "Сначала нижняя граница!\n");

    // if (argv[3]!=NULL)
    // {
    //     thread thread1(task2, f1, minNodes, maxNodes, true);
    //     thread thread2(task2, f2, minNodes, maxNodes, true);
    //     thread1.join();
    //     thread2.join();
    // } else
    // {
    //     task2(f1, minNodes, maxNodes, false);
    //     // task2(f2, minNodes, maxNodes, false);
    // }
    // cout << "Задание 2: Успешно!\n\n";

    // if (argv[3]!=NULL)
    // {
    //     thread thread1(task3, f1, minNodes, maxNodes, true);
    //     thread thread2(task3, f2, minNodes, maxNodes, true);
    //     thread1.join();
    //     thread2.join();
    // } else
    // {
    //     task3(f1, minNodes, maxNodes, false);
    //     // task3(f2, minNodes, maxNodes, false);
    // }
    assert(1==1==1);
    task(f1, minNodes, maxNodes, 100, true);
    // task(f2, minNodes, maxNodes, 100, false);
    // task3(f1, minNodes, maxNodes, false);
    
    // thread1 = thread(task3, f1, minNodes, maxNodes, name1, false);
    // thread2 = thread(task3, f2, minNodes, maxNodes, name2, false);

    // thread1.join();
	// thread2.join();
    cout << "Задание 3: Успешно!\n\n";

    return 0;
}