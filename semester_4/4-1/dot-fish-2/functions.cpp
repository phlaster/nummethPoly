#include "functions.hpp"

//exact (etalon) function value
double exact_y(double x){
    return cbrt(3*x*x*(x-1));
}

//returns the value of the first derivitive dy/dx
double dydx(double x, double y){
    //just in case
    if (y * x == 0){
        cerr << "X or Y equals to 0 => dy/dx -> Inf\nAborting the program :(\n";
        exit(1);
    }
    return x/y/y + y/x;
}

double error_max(vector<point> &vec){
    double maxerr = 0.0;
    for (auto p : vec){
        double err = fabs(p.y - exact_y(p.x));
        maxerr = (err>maxerr) ? err : maxerr;
    }
    return maxerr;
}

int max_error_index(vector<point> &vec){
    int index = 0;
    double maxerr = 0.0;
    for (int i = 0; i < vec.size(); i++){
        double err = fabs(vec[i].y - exact_y(vec[i].x));
        if(err>maxerr){
            index = i;
            maxerr = err;
        }
    }
    return index;
}