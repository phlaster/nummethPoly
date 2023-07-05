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


void task(double (*f)(double, bool), int minNodes, int maxNodes, int coeffMult, double distModulo, Buffer& buffer)
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

    buffer.append("nNodes,nInterpol,err_LagrUniform,err_LagrCheb,err_Herm,err_LagrUniformDev,err_LagrChebDev,err_HermDev\n");

    for (int nNodes=minNodes; nNodes<=maxNodes; nNodes++)
    {
        int nInterpol = coeffMult * nNodes;

        Graphic uniformNodes = tabulateFunction(f, left, right, nNodes);
        Graphic uniformNodes_DerNum = tabulateDerivativeNum(uniformNodes);
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



        Graphic uniformNodesDev = deviate(uniformNodes, distModulo);
        Graphic chebNodesDev = deviate(chebNodes, distModulo);
        Graphic uniformNodesDev_DerNum = tabulateDerivativeNum(uniformNodesDev);

        Graphic uniformLagrDev = lagrangeInterpol(uniformNodesDev, nInterpol);
        Graphic uniformHermDev = hermiteSpline(uniformNodesDev, uniformNodes_DerNum, nInterpol);
        Graphic chebLagrDev = lagrangeInterpol(chebNodesDev, nInterpol);

        Vec err_uniformLagrDev = errorProfile(uniformLagrDev.yVals, uniformTab.yVals);
        Vec err_uniformHermDev = errorProfile(uniformHermDev.yVals, uniformTab.yVals);
        Vec err_chebLagrDev = errorProfile(chebLagrDev.yVals, chebTab.yVals);


        buffer.append(nNodes);
        buffer.append(nInterpol);

        buffer.append(mean(err_uniformLagr));
        buffer.append(mean(err_uniformHerm));
        buffer.append(mean(err_chebLagr));
        
        buffer.append(mean(err_uniformLagrDev));
        buffer.append(mean(err_uniformHermDev));
        buffer.append(mean(err_chebLagrDev), true);

        write(
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
            err_chebLagr,
            uniformNodesDev,
            chebNodesDev,
            uniformNodesDev_DerNum,
            uniformLagrDev,
            uniformHermDev,
            chebLagrDev,
            err_uniformLagrDev,
            err_uniformHermDev,
            err_chebLagrDev
        );
    }
    
}

int main(int argc, char *argv[])
{
    int minNodes = stoi(argv[1]);
    int maxNodes = stoi(argv[2]);
    int interpolCoeff = stoi(argv[3]);
    double distModulo = fabs(stod(argv[4]));

    assert(minNodes < maxNodes && "Сначала нижняя граница!\n");
    Buffer b1, b2;
    thread thread1(task, f1, minNodes, maxNodes, interpolCoeff, distModulo, ref(b1));
    thread thread2(task, f2, minNodes, maxNodes, interpolCoeff, distModulo, ref(b2));
    thread1.join();
    thread2.join();

    ofstream outSummary_1(
        "CSVs/f1_summ_"+
        to_string(minNodes)+
        "-"+to_string(maxNodes)+
        "_"+to_string(distModulo)
        +"_.csv"
    );
    ofstream outSummary_2(
        "CSVs/f2_summ_"+
        to_string(minNodes)+
        "-"+to_string(maxNodes)+
        "_"+to_string(distModulo)
        +"_.csv"
    );

    if (!outSummary_1.is_open() || !outSummary_1.is_open())
        throw runtime_error("Не удалось открыть один из файлов для записи сводки!\n");
    outSummary_1 << b1.data;
    outSummary_2 << b2.data;
    outSummary_1.close();
    outSummary_2.close();

    cout << "Успешно!\n";
    return 0;
}