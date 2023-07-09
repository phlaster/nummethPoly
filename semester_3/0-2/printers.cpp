#include "functions.hpp"


void write(
    const Str& filename,

    const Graphic& uniformNodes,
    const Graphic& uniformNodes_DerNum,
    const Graphic& chebNodes,
    
    const Graphic& uniformTab,
    const Graphic& uniformTab_Noize,
    // const Graphic& chebTab,
    const Graphic& uniformLagr,
    const Graphic& uniformHerm,
    const Graphic& chebLagr,
    
    const Vec& err_uniformLagr,
    const Vec& err_uniformHerm,
    const Vec& err_chebLagr,

    const Graphic& uniformNodes_Noize,
    const Graphic& chebNodes_Noize,
    const Graphic& uniformNodesNoize_DerNum,
    
    const Graphic& uniformLagr_Noize,
    const Graphic& uniformHerm_Noize,
    const Graphic& chebLagr_Noize,
    
    const Vec& err_uniformLagr_Noize,
    const Vec& err_uniformHerm_Noize,
    const Vec& err_chebLagr_Noize
){
    int nNodes = uniformNodes.N;
    int nInterpol = uniformTab.N;

    ofstream outStream(filename);
    if (!outStream.is_open())
        throw runtime_error("Не удалось открыть файл для записи!\n");
    outStream // Заголовок .csv файла
        << "x_uniformNodes,y_uniformNodes,y_uniformNodes_DerNum,"
        << "y_uniformNodesNoize,y_uniformNodesNoize_DerNum,"
        << "x_chebNodes,y_chebNodes,y_chebNodes_Noize,"

        << "x_uniformTab,y_uniformTab,y_uniformTab_Der,"
        // << "x_chebTab,y_chebTab,"
        << "y_uniformLagr,y_chebLagr,y_uniformHerm,"
        << "y_uniformLagr_Noize,y_chebLagr_Noize,y_uniformHerm_Noize,"//
        << "err_uniformLagr,err_chebLagr,err_uniformHerm,"
        << "err_uniformLagr_Noize,err_chebLagr_Noize,err_uniformHerm_Noize\n";//
    for (int i = 0; i<nInterpol; i++){
        if (i>=nNodes)
            outStream << ",,,,,,,,";
        else
            outStream
                << uniformNodes.xVals[i] << "," << uniformNodes.yVals[i] << "," << uniformNodes_DerNum.yVals[i] << ","
                << uniformNodes_Noize.yVals[i] << "," << uniformNodesNoize_DerNum.yVals[i] << ","
                << chebNodes.xVals[i] << "," << chebNodes.yVals[i] << "," << chebNodes_Noize.yVals[i] << ",";
        outStream 
            << uniformTab.xVals[i] << "," << uniformTab.yVals[i] << "," << uniformTab_Noize.yVals[i] << ","
            // << chebTab.xVals[i] << "," << chebTab.yVals[i] << ","
            << uniformLagr.yVals[i] << "," << chebLagr.yVals[i] << "," << uniformHerm.yVals[i] << ","
            << uniformLagr_Noize.yVals[i] << "," << chebLagr_Noize.yVals[i] << "," << uniformHerm_Noize.yVals[i] << ","
            << err_uniformLagr[i] << "," << err_chebLagr[i] << "," << err_uniformHerm[i] << ","
            << err_uniformLagr_Noize[i] << "," << err_chebLagr_Noize[i] << "," << err_uniformHerm_Noize[i] << "\n";
    }
    outStream.close();
}

void write2_3(
    const Str& filename,

    const Graphic& uniformNodes,
    const Graphic& uniformNodes_DerNum,
    
    const Graphic& uniformTab,
    const Graphic& uniformTab_Der,

    const Graphic& uniformLagr,
    const Graphic& uniformHerm,
    
    const Vec& err_uniformLagr,
    const Vec& err_uniformHerm
){
    int nNodes = uniformNodes.N;
    int nInterpol = uniformTab.N;

    ofstream outStream(filename);
    if (!outStream.is_open())
        throw runtime_error("Не удалось открыть файл для записи!\n");
    outStream
        << "x_uniformNodes,y_uniformNodes,y_uniformNodes_DerNum,"
        << "x_uniformTab,y_uniformTab,y_uniformTab_Der,"
        << "y_uniformLagr,y_uniformHerm,"
        << "err_uniformLagr,err_uniformHerm\n";

    for (int i = 0; i<nInterpol; i++)
    {
        if (i>=nNodes)
            outStream << ",,,";
        else
            outStream
                << uniformNodes.xVals[i] << "," << uniformNodes.yVals[i] << "," << uniformNodes_DerNum.yVals[i] << ",";
        outStream 
            << uniformTab.xVals[i] << "," << uniformTab.yVals[i] << "," << uniformTab_Der.yVals[i] << ","
            << uniformLagr.yVals[i] << "," << uniformHerm.yVals[i] << ","
            << err_uniformLagr[i] << "," << err_uniformHerm[i] << "\n";
    }
    outStream.close();
}

void write4(
    const Str& filename,

    const Graphic& uniformNodes,
    const Graphic& chebNodes,
    
    const Graphic& uniformTab,
    const Graphic& uniformLagr,
    const Graphic& chebLagr,
    
    const Vec& err_uniformLagr,
    const Vec& err_chebLagr
){
    int nNodes = uniformNodes.N;
    int nInterpol = uniformTab.N;

    ofstream outStream(filename);
    if (!outStream.is_open())
        throw runtime_error("Не удалось открыть файл для записи!\n");
    outStream
        << "x_uniformNodes,y_uniformNodes,"
        << "x_chebNodes,y_chebNodes,"

        << "x_uniformTab,y_uniformTab,"
        << "y_uniformLagr,y_chebLagr,"

        << "err_uniformLagr,err_chebLagr\n";
    for (int i = 0; i<nInterpol; i++){
        if (i>=nNodes)
            outStream << ",,,,";
        else
            outStream
                << uniformNodes.xVals[i] << "," << uniformNodes.yVals[i] << ","
                << chebNodes.xVals[i] << "," << chebNodes.yVals[i] << ",";
        outStream 
            << uniformTab.xVals[i] << "," << uniformTab.yVals[i] << ","
            << uniformLagr.yVals[i] << "," << chebLagr.yVals[i] << ","
            << err_uniformLagr[i] << "," << err_chebLagr[i] << "\n";
    }
    outStream.close();
}

void write5(
    const Str& filename,
    const Graphic& uniformNodes_Noize,
    const Graphic& chebNodes_Noize,
    const Graphic& uniformNodes_Noize_DerNum,

    const Graphic& uniformTab,
    const Graphic& uniformTab_Der,

    const Graphic& uniformLagr_Noize,
    const Graphic& uniformHerm_Noize,
    const Graphic& chebLagr_Noize,

    const Vec& err_uniformLagr_Noize,
    const Vec& err_uniformHerm_Noize,
    const Vec& err_chebLagr_Noize
){
    int nNodes = uniformNodes_Noize.N;
    int nInterpol = uniformTab.N;

    ofstream outStream(filename);
    if (!outStream.is_open())
        throw runtime_error("Не удалось открыть файл для записи!\n");
    outStream
        << "x_uniformNodes_Noize,y_uniformNodes_Noize,y_uniformNodes_Noize_DerNum,"
        << "x_chebNodes_Noize,y_chebNodes_Noize,"

        << "x_uniformTab,y_uniformTab,y_uniformTab_Der"
        << "y_uniformLagr_Noize,y_uniformHerm_Noize,y_chebLagr_Noize,"
        << "err_uniformLagr_Noize,err_uniformHerm_Noize,err_chebLagr_Noize\n";
    for (int i = 0; i<nInterpol; i++){
        if (i>=nNodes)
            outStream << ",,,,,";
        else
            outStream
                << uniformNodes_Noize.xVals[i] << "," << uniformNodes_Noize.yVals[i] << "," << uniformNodes_Noize_DerNum.yVals[i] << ","
                << chebNodes_Noize.xVals[i] << "," << chebNodes_Noize.yVals[i] << ",";
        outStream 
            << uniformTab.xVals[i] << "," << uniformTab.yVals[i] << "," << uniformTab_Der.yVals[i] << ","
            << uniformLagr_Noize.yVals[i] << "," << uniformHerm_Noize.yVals[i] << "," << chebLagr_Noize.yVals[i] << ","
            << err_uniformLagr_Noize[i] << "," << err_uniformHerm_Noize[i] << "," << err_chebLagr_Noize[i] << "\n";
    }
    outStream.close();
}