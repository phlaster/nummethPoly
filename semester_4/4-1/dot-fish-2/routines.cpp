#include "headers/routines.hpp"
#include "headers/functions.hpp"
#include "headers/runge.hpp"

//limits
const double a = 1.1, b = 3.0;

void Baza(void){
    cout << "executing \"Baza\"\n";

    for (double h : {0.1, 0.05}){
        ofstream stream("data_baza_"+to_string(h)+".csv", ofstream::trunc);
        cout << "calculating for h="<<h<<endl;

        vector<point> res = RK2(exact_y(a), a, b, h, &dydx);

        stream << "x,y_h="<<h<<",y_exact,err_h="<<h<<endl;
        for(size_t i = 0; i < res.size(); i++){
            stream
            << res[i].x << ","
            << res[i].y << ","
            << exact_y(res[i].x) << ","
            << fabs(res[i].y - exact_y(res[i].x)) << endl;
        }
        stream.close();
    }

    ofstream stream("data_baza_8dots.csv",ofstream::trunc);
    stream << "h,h_square,max_err,nsteps\n";
    for (int i = 1; i <= 8; i++){
        double h = pow(10,-i);
        cout << "calculating for h =" << h << endl;
        auto [y_n, maxerr, nsteps] = RK2_maxerror(exact_y(a), a, b, h, &dydx, &exact_y);
        stream << h << "," << h*h << "," << maxerr << "," << nsteps << "\n" << flush;
    }
    stream.close();

    cout << "Baza done!\n";
    return;
}

void Minimum(void){
    cout << "Minimum"<<endl; 
    ofstream stream("data_min.csv",ofstream::trunc);

    stream << "eps,N,max_err\n";
    for (int i = 1; i <= 12; i++){
        double eps = pow(10,-i);
        double h = 1.0;
        while(1){            
            auto [y_n, err_n, steps_n  ] = RK2_maxerror(exact_y(a), a, b, h,   &dydx, &exact_y);
            auto [y_2n,err_2n,nsteps_2n] = RK2_maxerror(exact_y(a), a, b, h/2, &dydx, &exact_y);

            if(fabs(y_2n-y_n)/3.0 <= eps){
                stream << eps << "," << nsteps_2n << "," << err_2n << endl;
                break;
            }
            h /= 2.0;
        }
    }
    stream.close();
    cout << "Minimum done!\n";
    return;
}

void Dostatochno(void){
    cout << "Dostatochno"<<endl; 
    ofstream stream("data_dost.csv",ofstream::trunc);

    stream << "eps,N,max_err\n";
    for (int i = 1; i <= 8; i++){
        double eps = pow(10,-i);
        cout << "calculating for eps = " << eps << "... " << flush;
        auto [y_n, err, nsteps]  = RK2_adaptive(exact_y(a), a, b, 0.1, eps, &dydx, &exact_y);
        cout << "and writing to file..." << endl;
        stream << eps << "," << nsteps << "," << err << endl;
    }
    stream.close();
    cout << "...all done!\n";
    return;
}
