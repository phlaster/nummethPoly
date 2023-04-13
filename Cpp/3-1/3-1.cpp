#include <iostream>
#include <cmath>

using namespace std;


// Variant 23:
double polynom(double x)
{
    return\
            pow(x, 5) -\
        3.2*pow(x, 3) +\
        9.5*pow(x,2) -\
        7.0*x -\
        7.5;
}

double F(double x)
{
    return\
            pow(x,6)/6.0 -\
        3.2*pow(x,4)/4.0 +\
        9.5*pow(x,3)/3.0 -\
        3.5*pow(x,2) -\
        7.5*x;
}

double integralPrescise(double a=-0.7, double b=1.6)
{
    return F(b) - F(a); // Newton-Leibniz formula
}


double integral(int n, double a=-0.7, double b=1.6)
{
    double I = 0.0;
    double left = a;
    double h = (b-a)/n;
    double right = a + h;
    
    for(int i=1; i <= n; i++)
        {
            I += h * (polynom(left) + polynom(right)) / 2.0;
            left = right;
            right += h;
        }
    return I;
}



void integralBase(int nSteps = 2000)
{

    double I_presc = integralPrescise();
    double result;
    for(int n = 1; n<=nSteps; n++)
    {
        result = abs(I_presc - integral(n));
        cout << result << endl;
    }
}


void integralRunge(int lowestPower = -10)
{
    double ETA = 1/3.0;
    double eps;
    double delta = 1.0;
    
    for(int eps_power = -1; eps_power > lowestPower; eps_power--)
    {
        delta = 1.0;
        int n = 1;
        eps = pow(10, eps_power);
        while(delta > eps)
            {        
                delta = ETA*abs(integral(2*n) - integral(n));
                n *= 2;
            }
        cout << n << ", " << delta << ", " << eps << endl;
    }
}

int main()
{
    // float N;
    // cout << "Enter number of steps:\n";
    // cin >> N;
    integralBase();
    integralRunge();
    // cout << "Prescise value = " << integralPrescise() << endl;
    // cout << "Integral = " << I << endl;
}
