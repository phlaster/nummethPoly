#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

// Variant 23:
double f(double x)
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
double integral_rectangular(int n_steps, double a=-0.7, double b=1.6) {
    double h = (b-a)/n_steps;
    double I = 0.0;
    for(int i=0; i < n_steps; i++){
        double left = a + i*h;
        I += f(left) + f(left + h);
    }
    return (h/2) * I;
}

double integral_trapez(int n_steps, double a=-0.7, double b=1.6){
    double h = (b-a)/n_steps;
    double I = (f(a) + f(b)) / 2;
    for(int i=1; i < n_steps; i++)
        I += f(a + i*h);
    return h * I;
}

pair<int, double> integral_trapez_runge(double eps, double a=-0.7, double b=1.6) {
    int n_steps = 1;
    double I1 = integral_trapez(n_steps, a, b);
    n_steps *= 2;
    double I2 = integral_trapez(n_steps, a, b);
    while (fabs(I2 - I1) / 3 > eps) {
        n_steps *= 2;
        I1 = I2;
        I2 = integral_trapez(n_steps, a, b);
    }
    return make_pair(n_steps, I2);
}

double integral_trapez_adaptive_runge(double eps, int& n_steps, double a=-0.7, double b=1.6) {
    double h = b-a;
    double I1 = integral_trapez(1, a, b);
    double I2 = integral_trapez(2, a, b);
    double delta = 1./3. * fabs(I2 - I1);

    if (delta <= eps * h / 2)
    {
        n_steps++;
        return I2;
    }else{
        return integral_trapez_adaptive_runge(eps, n_steps, a, (a+b)/2.) +\
               integral_trapez_adaptive_runge(eps, n_steps, (a+b)/2., b);
    }
}

void int_convergence_1(int maxiters){
    double I_r, I_t, I_true = integralPrescise();
    int di = 0;
    cout << "n_steps,err_rect,err_trapez,delta_err\n";
    for (int i = 2; i <= maxiters; i *= 2){
        I_r = integral_rectangular(i);
        I_t = integral_trapez(i);

        cout << i  << "," << fabs(I_r - I_true) << "," << fabs(I_t - I_true) << "," << fabs(I_r - I_t) << '\n';
    }
    cout << endl;
}

void int_convergence_2(int minpower=12){
    double I_true = integralPrescise();
    cout << "n_steps,eps,err\n";
    for (int i = 1; i <= minpower; i += 1){
        double eps = pow(10, -i);

        auto [nsteps, I] = integral_trapez_runge(eps);
        cout << nsteps  << "," << eps << "," << fabs(I - I_true) << '\n';
    }
    cout << endl;
}

void int_convergence_3(int minpower=12){
    double I_true = integralPrescise();
    cout << "n_steps,eps,err\n";
    for (int i = 1; i <= minpower; i += 1){
        double eps = pow(10, -i);
        int nsteps = 1;
        double I = integral_trapez_adaptive_runge(eps, nsteps);
        cout << nsteps  << "," << eps << "," << fabs(I - I_true) << '\n';
    }
    cout << endl;
}


int main(void)
{
    cout << "1:\n";
    int_convergence_1(pow(2, 25));

    cout << "\n\n2:\n";
    int_convergence_2();

    cout << "\n\n3:\n";
    int_convergence_3();
    return 0;
}
