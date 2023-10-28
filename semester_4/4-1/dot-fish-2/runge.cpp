#include "Runge.hpp"

double RK2_maxerror(double y0, double a, double b, double h, double (*func)(double, double), double (*etalon)(double)){
    double maxerr = 0;
    double x = a, y = y0;
    const double hh = h*0.5;
    while(x <= b){
        double _y = y + hh*func(x,y);
        y += h*func(x+hh, _y);
        x += h;
        double err = fabs(etalon(x) - y);
        maxerr = err > maxerr ? err : maxerr;
    }
    return maxerr;
}

vector<point> RK2(double y0, double a, double b, double h, double (*func)(double, double))
/*returns a vector of x,y pair values calculated using modified Euler (middle-point) method 
 *h - integration step, func - right hand of the equation for the first derivitive dy/dx*/{
    vector<point> Y = {{a,y0}};
    double x = a, y = y0;
    const double hh = h*0.5;
    while(x <= b){
        double _y = y + hh*func(x,y);
        y += h*func(x+hh, _y);
        x += h;
        Y.push_back({x,y});
    }
    return Y;
}

vector<point> RK2_adaptive(double y0, double a, double b, double h0, double p, double (*func)(double, double)){
    vector<point> Y = {{a,y0}};
    double x = a, y = y0;
    double h = h0;
    while(x <= b){
        while(1){
            double hhi = h*0.5, hhj = h*0.25;
            double _yi = y + hhi*func(x,y), _yj = y + hhj*func(x,y);
            double yi = y + h*func(x+hhi, _yi), yj = y + hhi*func(x+hhj, _yj);
            _yj = yj + hhj*func(x,yj);
            yj += hhi*func(x+hhj, _yj);
            //cout << "<" << fabs(yi - yj)/3.0 << ">\n";
            if(fabs(yi - yj)/3.0 < p){
                y = yi;
                break;
            }
            h = hhi;
        }
        x += h;
        Y.push_back({x,y});
    }
    return Y;
}