// Вариант 8 Б
#include "functions.hpp"
#include <fstream>
#include <iomanip>

const Vec F1_INTERVAL = {0.5, 2.73205};
const Vec F2_INTERVAL = {0.6, 1.5};

const double X1 = 1.732050807568877;
const double X2 = 0.947404881615152;


void errorConverg(ans    (*method)(double (*f)(double), double, Vec, double),
                  double (*f)     (double),
                  string filename)
{
    ofstream outStream(filename);
    if (!outStream.is_open()) {
        throw std::runtime_error("Не удалось открыть файл для записи!");
    }

    outStream << "eps,err,root,steps\n";
    for(int i=1; i <=15; i++)
    {
        double presc = pow(10, -i);
        ans res;
        if (f==f1)
            res = method(f, X1, F1_INTERVAL, presc);
        else
            res = method(f, X2, F2_INTERVAL, presc);
        
        outStream << scientific << setprecision(0) << presc << ","
                  << scientific << setprecision(1) << res.absErr << ","
                  << fixed << setprecision(i+1) << res.root << ","
                  << res.nSteps << "\n";
    }
    outStream.close();
}


int main()
{
    errorConverg(bisection, f1, "bisect1.csv");
    errorConverg(bisection, f2, "bisect2.csv");

    errorConverg(newton, f1, "newton1.csv");
    errorConverg(newton, f2, "newton2.csv");
    return 0;
}