double integral_trapez(int n_steps, double a, double b){
    double h = (b-a)/n_steps;
    double I = (f(a) + f(b)) / 2;
    for(int i=1; i < n_steps; i++)
        I += f(a + i*h);
    return h * I;
}