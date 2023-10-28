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
