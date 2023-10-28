#include "routines.hpp"

//limits
const double a = 1.1, b = 3.0;

void Baza(void){
    cout << "executing \"Baza\" part of the assigment\n";
    ofstream 
        stream("data_baza_a.csv", ofstream::trunc),
        stream2("data_baza_b.csv", ofstream::trunc);
    // .IA.
    cout << "calculating for h = 0.1 ...\n";
    vector<point> numsolution = RK2(exact_y(a), a, b, 0.1,&dydx);
    //
    cout << "writing to file ...\n";
    stream << "x,y_h=0.1,y_exact,err_h=0.1\n";
    for(int i = 0; i < numsolution.size(); i++){
        stream << numsolution[i].x << "," << numsolution[i].y << "," << exact_y(numsolution[i].x) << "," << fabs(numsolution[i].y - exact_y(numsolution[i].x)) << "\n";
    }
    stream.close();
    // .IB.
    cout << "calculating for h = 0.05 ...\n";
    numsolution = RK2(exact_y(a), a, b, 0.05,&dydx);
    //
    cout << "writing to file ...\n";
    stream2 << "x,y_h=0.05,y_exact,err_h=0.05\n";
    for(int i = 0; i < numsolution.size(); i++){
        stream2 << numsolution[i].x << "," << numsolution[i].y << "," << exact_y(numsolution[i].x) << "," << fabs(numsolution[i].y - exact_y(numsolution[i].x)) << "\n";
    }
    stream2.close();
    cout << "-------------------\n";
    //Построить зависимость (No2) ошибки от шага. На график нанести линию h 2 (почему?)
    ofstream stream3("data_baza_c.csv",ofstream::trunc);
    stream3 << "h,h_square,max_err\n";
    for (int i = 1; i < 9; i++){
        double h = pow(10,-i);
        cout << "calculating for h = " << h << "... ";
        double maxerr = RK2_maxerror(exact_y(a), a, b, h,&dydx, &exact_y);
        cout << "and writing to file...\n";
        stream3 << h << "," << h*h << "," << maxerr << "\n";
    }
    stream3.close();
    //-----------------------------
    cout << "...all done!\n";
    return;
}

void Minimum(void){
    cout << "executing \"Minimum\" part of the assigment\n"; 
    //Построить зависимость 3, 4
    ofstream stream3("data_min.csv",ofstream::trunc);
    stream3 << "p,N,max_err\n";
    for (int i = 1; i < 9; i++){
        double p = pow(10,-i);
        double h = 1.0;
        while(1){
            cout << "calculating for p = " << p << "... ";
            
            vector<point> numsolution  = RK2_adaptive(exact_y(a), a, b, h, p, &dydx);
            vector<point> numsolution2  = RK2_adaptive(exact_y(a), a, b, h*0.5, p, &dydx);

            //int index = max_error_index(numsolution);
            double yi = numsolution[numsolution.size()].y, yj = numsolution2[numsolution2.size()].y;
            if(fabs(yi - yj)/3.0 < p){
                cout << "and writing to file...\n";
                stream3 << p << "," << numsolution.size() << "," << error_max(numsolution) << "\n";
                break;
            }
            h = 0.5*h;
        }
        
    }
    stream3.close();
    //-----------------------------
    cout << "...all done!\n";
    return;
}

void Dostatochno(void){
    cout << "executing \"Dostatochno\" part of the assigment\n"; 
    //Построить зависимость 3, 4
    ofstream stream3("data_dost.csv",ofstream::trunc);
    stream3 << "p,N,max_err\n";
    for (int i = 1; i < 9; i++){
        double p = pow(10,-i);
        cout << "calculating for p = " << p << "... ";
        vector<point> numsolution  = RK2_adaptive(exact_y(a), a, b, 1.0, p, &dydx);
        double maxerr = error_max(numsolution);
        cout << "and writing to file...\n";
        stream3 << p << "," << numsolution.size() << "," << maxerr << "\n";
    }
    stream3.close();
    //-----------------------------
    cout << "...all done!\n";
    return;
}

