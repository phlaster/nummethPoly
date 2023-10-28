double integral_rectangular(int n_steps, double a, double b) {
    double h = (b-a)/n_steps;
    double I = 0.0;
    for(int i=0; i < n_steps; i++){
        double left = a + i*h;
        I += f(left) + f(left + h);
    }
    return (h/2) * I;
}