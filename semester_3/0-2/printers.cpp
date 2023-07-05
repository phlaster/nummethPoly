#include "functions.hpp"


void write(
    const Str& filename,

    const Graphic& uniformNodes,
    const Graphic& uniformNodes_DerNum,
    const Graphic& chebNodes,
    
    const Graphic& uniformTab,
    const Graphic& uniformTab_Der,
    const Graphic& chebTab,
    const Graphic& uniformLagr,
    const Graphic& uniformHerm,
    const Graphic& chebLagr,
    
    const Vec& err_uniformLagr,
    const Vec& err_uniformHerm,
    const Vec& err_chebLagr,
    //
    const Graphic& uniformNodesDev,
    const Graphic& chebNodesDev,
    const Graphic& uniformNodesDev_DerNum,
    
    const Graphic& uniformLagrDev,
    const Graphic& uniformHermDev,
    const Graphic& chebLagrDev,
    
    const Vec& err_uniformLagrDev,
    const Vec& err_uniformHermDev,
    const Vec& err_chebLagrDev
){
    int nNodes = uniformNodes.N;
    int nInterpol = uniformTab.N;

    ofstream outStream(filename);
    if (!outStream.is_open())
        throw runtime_error("Не удалось открыть файл для записи!\n");
    outStream // Заголовок .csv файла
        << "x_uniformNodes,y_uniformNodes,y_uniformNodes_DerNum,"
        << "y_uniformNodesDev,y_uniformNodesDev_DerNum,"
        << "x_chebNodes,y_chebNodes,y_chebNodesDev,"

        << "x_uniformTab,y_uniformTab,y_uniformTab_Der,"
        << "x_chebTab,y_chebTab,"
        << "y_uniformLagr,y_chebLagr,y_uniformHerm,"
        << "y_uniformLagrDev,y_chebLagrDev,y_uniformHermDev,"//
        << "err_uniformLagr,err_chebLagr,err_uniformHerm,"
        << "err_uniformLagrDev,err_chebLagrDev,err_uniformHermDev\n";//
    for (int i = 0; i<nInterpol; i++){
        if (i>=nNodes)
            outStream << ",,,,,,,,";
        else
            outStream
                << uniformNodes.xVals[i] << "," << uniformNodes.yVals[i] << "," << uniformNodes_DerNum.yVals[i] << ","
                << uniformNodesDev.yVals[i] << "," << uniformNodesDev_DerNum.yVals[i] << ","
                << chebNodes.xVals[i] << "," << chebNodes.yVals[i] << "," << chebNodesDev.yVals[i] << ",";
        outStream 
            << uniformTab.xVals[i] << "," << uniformTab.yVals[i] << "," << uniformTab_Der.yVals[i] << ","
            << chebTab.xVals[i] << "," << chebTab.yVals[i] << ","
            << uniformLagr.yVals[i] << "," << chebLagr.yVals[i] << "," << uniformHerm.yVals[i] << ","
            << uniformLagrDev.yVals[i] << "," << chebLagrDev.yVals[i] << "," << uniformHermDev.yVals[i] << ","
            << err_uniformLagr[i] << "," << err_chebLagr[i] << "," << err_uniformHerm[i] << ","
            << err_uniformLagrDev[i] << "," << err_chebLagrDev[i] << "," << err_uniformHermDev[i] << "\n";
    }
    outStream.close();
}