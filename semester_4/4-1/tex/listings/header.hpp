//1
typedef struct point {
    double x;
    double y;
} point;
//2
double exact_y(const double x);
//3
double dydx(const double x, double y);
//4
vector<point> RK2(
    const double y0,
    const double a,
    const double b,
    const double h,
    double (*f)(double, double)
);
//5
tuple<double, double, unsigned> RK2_adaptive(
    const double y0,
    const double a,
    const double b,
    const double h_start,
    const double eps,
    double (*f)(double, double),
    double (*exact_solution)(double)
);
