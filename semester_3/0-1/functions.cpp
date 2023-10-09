#include "functions.hpp"

double f1(double x){ return pow(x, 4) - pow(x, 3) -2*pow(x, 2) + 3*x - 3; }

double f1d1(double x){ return 4*pow(x, 3) -3*pow(x, 2) - 4*x + 3; }

double f1d2(double x){ return 12*pow(x, 2) - 6*x - 4; }

double f2(double x){ return 3 * exp(x) - 5*x - 3; }

double f2d1(double x){ return 3 * exp(x) - 5; }

double f2d2(double x){ return 3 * exp(x); }

int sign(double x){
    return x>0.0 ?  1 :
           x<0.0 ? -1 : 0;
}

void errorConverg(double exact, Vec interval,
                  ans    (*method)(double (*f)(double), double, Vec, double),
                  double (*f)     (double),
                  string filename)
{
    ofstream outStream(filename);
    if (!outStream.is_open()) {
        throw runtime_error("Не удалось открыть файл для записи!");
    }

    outStream << "eps,err,root,steps\n";
    for(int i=1; i <=15; i++)
    {
        double presc = pow(10, -i);
        ans res = method(f, exact, interval, presc);
        
        outStream << scientific << setprecision(0) << presc << ","
                  << scientific << setprecision(1) << res.absErr << ","
                  << fixed << setprecision(i+1) << res.root << ","
                  << res.nSteps << "\n";
    }
    outStream.close();
}

void brokenConverg(double exact, double x0, Vec interval, string filename){
    ofstream outStream(filename);
    if (!outStream.is_open()) {
        throw runtime_error("Не удалось открыть файл для записи!");
    }

    outStream << "eps,err,root,steps\n";
    for(int i=1; i <=15; i++)
    {
        double presc = pow(10, -i);
        ans res = newton_broken(exact, interval, x0, presc);
        
        outStream << scientific << setprecision(0) << presc << ","
                  << scientific << setprecision(1) << res.absErr << ","
                  << fixed << setprecision(i+1) << res.root << ","
                  << res.nSteps << "\n";
    }
    outStream.close();
}