#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

// Variant 23:
double polynom(double x)
{
    return  pow(x,5) -\
        3.2*pow(x,3) +\
        9.5*pow(x,2) -\
        7.0*x -\
        7.5;
}

// Antiderivative of given polynom
double F(double x)
{
    return  pow(x,6)/6.0 -\
        3.2*pow(x,4)/4.0 +\
        9.5*pow(x,3)/3.0 -\
        3.5*pow(x,2) -\
        7.5*x;
}

// Newton-Leibniz formula
double integralPrescise(double a=-0.7, double b=1.6)
{
    return F(b) - F(a);
}

//Numeric calculations within given segment
double integral(int n, double a=-0.7, double b=1.6)
{
    double I = 0.0;
    double left = a;
    double dx = (b-a)/n;
    double right = a + dx;
    for(int i=1; i <= n; i++)
        {
            I += dx * (polynom(left) + polynom(right)) / 2.0;
            left = right;
            right += dx;
        }
    return I;
}

void integralBase(int nSteps = 2000)
{
    double I_presc = integralPrescise();
    ofstream outfile("baseResult.csv");
    if (outfile.is_open()) {
        for(int n = 1; n<=nSteps; n++)
        {
            outfile << abs(I_presc-integral(n)) << endl;
        }
        outfile.close();
        cout << "baseResult.csv has been written!" << endl;
    }
    else cerr << "Couldn't write to file baseResult.csv!" << endl;
}


void integralRunge(int lowestPower = -10)
{
    double ETA = 1/3.0; // Method constant for trapez method
    ofstream outfile("rungeResult.csv");
    if (outfile.is_open()) {
        for(int eps_power = -1; eps_power >= lowestPower; eps_power--)
        {
            double delta = 1.0; // Initial *big* error
            int n = 1;
            double eps = pow(10, eps_power); //higher prescision at each iteration
            while(delta > eps) // Runge method
            {        
                delta = ETA*abs(integral(2*n) - integral(n));
                n *= 2;
            }
            outfile << n/2 << "," << delta << "," << eps << endl;
        }
        outfile.close();
        cout << "rungeResult.csv has been written!" << endl;
    }
    else cerr << "Couldn't write to file rungeResult.csv!" << endl;
}

int main(void)
{

    integralBase();
    integralRunge();
    return 0;
}
