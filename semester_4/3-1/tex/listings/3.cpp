pair<int, double> integral_trapez_runge(double eps, double a, double b) {
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