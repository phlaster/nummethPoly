#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>
#include "Runge.hpp"


double analytical_y(double x)
//analytical (etalon) function value
{
    return std::cbrt(3*x*x*(x-1));
}

double foo(double x, double y)
//returns the value of the first derivitive dy/dx 
{
    if (y == 0)
    //just in case
    {
        std::cout << "Y equals to 0 => dy/dx -> Inf\nAborting the program :(\n";
        exit(1);
    }
    if (x == 0)
    //just in case    
    {
        std::cout << "X equals to 0 => dy/dx -> Inf\nAborting the program :(\n";
        exit(1);
    }

    return x/y/y + y/x;
}

double error_max(std::vector<point> &vec)
{
    double maxerr = 0.0;
    for (auto p : vec)
    {
        double err = std::abs(p.y - analytical_y(p.x));
        maxerr = (err>maxerr) ? err : maxerr;
    }
    return maxerr;
}

int max_error_index(std::vector<point> &vec)
{
    int index = 0;
    double maxerr = 0.0;
    for (int i = 0; i < vec.size(); i++)
    {
        double err = std::abs(vec[i].y - analytical_y(vec[i].x));
        if(err>maxerr)
        {
            index = i;
            maxerr = err;
        }
    }
    return index;
}


//limits of integration:
const double a = 1.1, b = 3.0;

//Base part of the assignment (+0)
void Baza(void)
{
    std::cout << "executing \"Baza\" part of the assigment\n";
    std::ofstream 
        stream("data_baza_a.csv", std::ofstream::trunc),
        stream2("data_baza_b.csv", std::ofstream::trunc);
    // .IA.
    std::cout << "calculating for h = 0.1 ...\n";
    std::vector<point> numsolution = RK2(analytical_y(a), a, b, 0.1,&foo);
    //
    std::cout << "writing to file ...\n";
    stream << "x, y (h = 0.1), y_{analytical}, error (h = 0.1)\n";
    for(int i = 0; i < numsolution.size(); i++)
    {
        stream << numsolution[i].x << ", " << numsolution[i].y << ", " << analytical_y(numsolution[i].x) << ", " << std::abs(numsolution[i].y - analytical_y(numsolution[i].x)) << "\n";
    }
    stream.close();
    // .IB.
    std::cout << "calculating for h = 0.05 ...\n";
    numsolution = RK2(analytical_y(a), a, b, 0.05,&foo);
    //
    std::cout << "writing to file ...\n";
    stream2 << "x, y (h = 0.05), y_{analytical}, error (h = 0.05)\n";
    for(int i = 0; i < numsolution.size(); i++)
    {
        stream2 << numsolution[i].x << ", " << numsolution[i].y << ", " << analytical_y(numsolution[i].x) << ", " << std::abs(numsolution[i].y - analytical_y(numsolution[i].x)) << "\n";
    }
    stream2.close();
    std::cout << "-------------------\n";
    //Построить зависимость (No2) ошибки от шага. На график нанести линию h 2 (почему?)
    std::ofstream stream3("data_baza_c.csv",std::ofstream::trunc);
    stream3 << "h, h^2, max error\n";
    for (int i = 1; i < 9; i++)
    {
        double h = std::pow(10,-i);
        std::cout << "calculating for h = " << h << "... ";
        double maxerr = RK2_maxerror(analytical_y(a), a, b, h,&foo, &analytical_y);
        std::cout << "and writing to file...\n";
        stream3 << h << ", " << h*h << ", " << maxerr << "\n";
    }
    stream3.close();
    //-----------------------------
    std::cout << "...all done!\n";
    return;
}

//Minimum part of the assignment (+1)
void Minimum(void)
{
    std::cout << "executing \"Minimum\" part of the assigment\n"; 
    //Построить зависимость 3, 4
    std::ofstream stream3("data_min.csv",std::ofstream::trunc);
    stream3 << "p, N, max error\n";
    for (int i = 1; i < 9; i++)
    {
        double p = std::pow(10,-i);
        double h = 1.0;
        while(1)
        {
            std::cout << "calculating for p = " << p << "... ";
            std::vector<point> numsolution  = RK2_adaptive(analytical_y(a), a, b, h, p, &foo);
            std::vector<point> numsolution2  = RK2_adaptive(analytical_y(a), a, b, h*0.5, p, &foo);
            //int index = max_error_index(numsolution);
            double yi = numsolution[numsolution.size()].y, yj = numsolution2[numsolution2.size()].y;
            if(std::abs(yi - yj)/3.0 < p)
            {
                std::cout << "and writing to file...\n";
                stream3 << p << ", " << numsolution.size() << ", " << error_max(numsolution) << "\n";
                break;
            }
            h = 0.5*h;
        }
        
    }
    stream3.close();
    //-----------------------------
    std::cout << "...all done!\n";
    return;
}

//Dostatochno part of the assignment (+1)
void Dostatochno(void)
{
    std::cout << "executing \"Dostatochno\" part of the assigment\n"; 
    //Построить зависимость 3, 4
    std::ofstream stream3("data_dost.csv",std::ofstream::trunc);
    stream3 << "p, N, max error\n";
    for (int i = 1; i < 9; i++)
    {
        double p = std::pow(10,-i);
        std::cout << "calculating for p = " << p << "... ";
        std::vector<point> numsolution  = RK2_adaptive(analytical_y(a), a, b, 1.0, p, &foo);
        double maxerr = error_max(numsolution);
        std::cout << "and writing to file...\n";
        stream3 << p << ", " << numsolution.size() << ", " << maxerr << "\n";
    }
    stream3.close();
    //-----------------------------
    std::cout << "...all done!\n";
    return;
}

int main(int argc, char** argv)
{  
    if (argc == 1)
    {
        std::cout << "flags to execute tasks:\n";
        std::cout << "\t--all\n";
        std::cout << "\t--baza\n";
        std::cout << "\t--minimum / --min\n";
        std::cout << "\t--dostatochno / --dost\n";
        return 0;
    }
    bool 
        do_baza = false,
        do_minimum = false,
        do_dostatochno = false;
    for(int i = 1; i < argc; i++)
    {
        std::string s = argv[i];
        if (!s.compare("--all"))
        {
            Baza();
            Minimum();
            Dostatochno();
            return 0;
        }
        else if (!s.compare("--baza"))
        {
            do_baza = true;
        }
        else if (!(s.compare("--min")) || !(s.compare("--minimum")))
        {
            do_minimum = true;
        }
        else if (!(s.compare("--dostatochno")) || !(s.compare("--dost")))
        {
            do_dostatochno = true;
        }
    }
    if (do_baza)
        Baza();
    if (do_minimum)
        Minimum();
    if (do_dostatochno)
        Dostatochno();
    
    return 0;
}