#include "Runge.hpp"

double RK2_maxerror(
    double y0, double a, double b, double h, double (*dydx)(double, double), double (*y_exact)(double)
){
    double maxerr = 0.0;
    double x = a, y = y0;
    const double hh = h * 0.5;
    while(x <= b){
        double _y = y + hh * dydx(x,y);
        y += h * dydx(x+hh, _y);
        x += h;
        double err = fabs(y_exact(x) - y);
        maxerr = err > maxerr ? err : maxerr;
    }
    return maxerr;
}


/*
    returns a vector of x,y pair values calculated using modified Euler (middle-point) method 
    *h - integration step,
    func - right hand of the equation for the first derivitive dy/dx
*/
vector<point> RK2(
    double y0, double a, double b, double h, double (*dydx)(double, double)
){
    vector<point> Y = {{a,y0}};
    double x = a, y = y0;
    const double hh = h*0.5;
    while(x <= b){
        double _y = y + hh * dydx(x,y);
        y += h * dydx(x+hh, _y);
        x += h;
        Y.push_back({x,y});
    }
    return Y;
}

vector<point> RK2_adaptive(
    double y0, double a, double b, double h0, double p, double (*dydx)(double, double)
){
    vector<point> Y = {{a,y0}};
    double x = a, y = y0;
    double h = h0;
    while(x <= b){
        while(1){
            double hhi = h * 0.5;
            double _yi = y + hhi * dydx(x,y);
            double  yi = y + h   * dydx(x+hhi, _yi);

            double hhj = h * 0.25;
            double _yj = y + hhj * dydx(x,y);
            double  yj = y + hhi * dydx(x+hhj, _yj);

            _yj = yj + hhj*dydx(x,yj);
            yj += hhi*dydx(x+hhj, _yj);
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