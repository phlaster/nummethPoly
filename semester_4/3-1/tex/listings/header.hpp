double integral_rectangular(int n_steps, double a, double b);
double integral_trapez(int n_steps, double a, double b);
pair<int, double> integral_trapez_runge(double eps, double a, double b);
double integral_trapez_adaptive_runge(double eps, int& n_steps, double a, double b);
