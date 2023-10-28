double dy_dx(double x, double y){
    return x/(y * y) + y/x;
}

double y(double x){
    return cbrt(3 * x * x * (x-1));
}

double y_wave_next(double y_i, double h_i, double f){
    return y_i + h_i/2 * f;
}

double y_next(double y_i, double h_i, double f){
    return y_i + h_i * f;
}

void calculate(double a, double b, double h_i){
    double y_0 = y(a);
    double h = b-a;
    double x_i = a;
    int n_steps = round(h/h_i);

    cout << "x_i,y_appr_i,y_i,err\n";
    for (int i = 0; i < n_steps; i++){
        double y_wave = y_wave_next(y_0, h_i, dy_dx(x_i, y_0));
        double y_i = y_next(y_0, h_i, dy_dx(x_i + h_i/2, y_wave));
        y_0 = y_i;
        x_i = a + i * h_i;

        cout << fixed << setprecision(6);
        cout << x_i << ',' << y_i << ',' << y(x_i) << "," << fabs (y_i - y(x_i)) << '\n';

    }
}
