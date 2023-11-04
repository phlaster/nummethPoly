#include "headers/runge.hpp"
#include <thread>

double eta(double (*f)(double, double), const point p, const double h){
    double h_div_2 = h/2.0;
    double _x = p.x + h_div_2;
    double _y = p.y + h_div_2 * f(p.x, p.y);
    return f(_x, _y);
}

pair<double, int> RK2_recursive(
    double (*f)(double, double),
    double (*exact)(double),
    const double y0,
    const double x0,
    const double b,
    const double eps
) {
    double h = b - x0;
    double h_div_2 = h/2.0;
    // I
    double eta_i = eta(f, {x0, y0}, h);
    double Delta_y_i = h * eta_i;
    double y_i = y0+ Delta_y_i;
    // I.1
    double eta_j1 = eta(f, {x0, y0}, h_div_2);
    double Delta_y_j1 = h_div_2 * eta_j1;
    double y_j1 = y0 + Delta_y_j1;
    double x_j1 =x0+h_div_2;
    double err_05 = fabs(exact(x_j1) - y_j1);
    // I.2
    double eta_j2 = eta(f, {x_j1, y_j1}, h_div_2);
    double Delta_y_j2 = h_div_2 * eta_j2;
    double y_j2 = y_j1 + Delta_y_j2;
    double x_j2 =x0+h;
    double err_10 = fabs(exact(x_j2) - y_j2);

    double R = fabs(y_j2 - y_i)/3.0;
    if (R <= eps){
        return make_pair(max(err_05, err_10), 2);
    } else {
        auto [l_maxerr, l_n] = RK2_recursive(f, exact, y0, x0, x_j1, eps);
        auto [r_maxerr, r_n] = RK2_recursive(f, exact, exact(x_j1), x_j1, b, eps); // !!!!!!!!

        return make_pair(max(l_maxerr, r_maxerr), l_n+r_n);
    }
}


vector<point> RK2(
    const double y0,
    const double a,
    const double b,
    const double h,
    double (*f)(double, double)
){
    if (h<=1e-7){
        throw invalid_argument("For this routine h>1e-7!");
    }
    const int n = round((b-a)/h);
    vector<point> Y(n);
    double x_i = a, y_i = y0;
    for (int i=0; i<n; i++){
        Y[i] = {x_i, y_i};
        double eta_i = eta(f, Y[i], h);
        double Delta_y_i = h * eta_i;
        x_i = a + (i+1)*h;
        y_i += Delta_y_i;
    }
    Y.push_back({x_i, y_i}); // n+1 points for n intervals
    return Y;
}

tuple<double, double, int> RK2_maxerror(
    const double y0,
    const double a,
    const double b,
    const double h,
    double (*f)(double, double),
    double (*y_exact)(double)
){
    const int n = round((b-a)/h) + 1;
    double x_i = a,
           y_i = y0,
           maxerr = 0.0;
    for (int i=0; i<n; i++){
        double eta_i = eta(f, {x_i,y_i}, h);
        double Delta_y_i = h * eta_i;
        x_i = a + (i+1)*h;
        y_i += Delta_y_i;

        double err = fabs(y_exact(x_i) - y_i);
        maxerr = err > maxerr ? err : maxerr;
    }
    return make_tuple(y_i, maxerr, n);
}






// tuple<double, double, ulong> RK2_adaptive(
//     const double y0,
//     const double a,
//     const double b,
//     const double h0,
//     const double eps,
//     double (*dydx)(double, double),
//     double (*y_exact)(double)
// ){
//     double x = a, y = y0;
//     double h = h0;
//     ulong i = 0;
//     double maxerr = 0.0;
//     while(x <= b){
//         while(1){
//             i++;
//             double hhi = h/2;
//             double _yi = y + hhi * dydx(x,y);
//             double  yi = y + h   * dydx(x+hhi, _yi);
//             double  hhj = h/4;
//             double _yj1 = y + hhj * dydx(x,y);
//             double  yj1 = y + hhi * dydx(x+hhj, _yj1);
//             double _yj2 = yj1 + hhj * dydx(x+hhi, yj1);
//             double  yj2 = yj1 + hhi * dydx(x+hhi, _yj2);
//             if(fabs(yi - yj2)/3.0 < eps){
//                 y = yj2;
//                 break;
//             }
//             h = hhi;
//         }
//         x += h;
//         double err = fabs(y_exact(x) - y);
//         maxerr = err > maxerr ? err : maxerr;
//     }
//     return make_tuple(y, maxerr, i);
// }

// tuple<double, double, ulong> RK2_adaptive(
//     const double y0,
//     const double a,
//     const double b,
//     const double h0,
//     const double eps,
//     double (*dydx)(double, double),
//     double (*y_exact)(double)
// ){
//     double x = a, y = y0;
//     double h = h0;
//     ulong i = 0;
//     double maxerr = 0.0;
//     while(x <= b){
//         while(1){
//             i++;
//             double hd2 = h/2;
//             double yi  = y   + h   * eta(dydx, {x,y},   h);
//             double yj1 = y   + hd2 * eta(dydx, {x,y},   hd2);
//             double yj2 = yj1 + hd2 * eta(dydx, {x+hd2,yj1}, hd2);

//             if(fabs(yi - yj2)/3.0 < eps){
//                 y = yj2;
//                 break;
//             }
//             h = hd2;
//         }
//         x += h;
        
//         double err = fabs(y_exact(x) - y);
//         maxerr = err > maxerr ? err : maxerr;
//         h = h0;
//     }
//     return make_tuple(y, maxerr, i);
// }

tuple<double, double, int> RK2_adaptive(
    const double y0,
    const double a,
    const double b,
    const double h0,
    const double eps,
    double (*f)(double, double),
    double (*y_exact)(double)
){
    int n = round((b-a)/h0) + 1,
        steps = 0;
    double x = a,
           y = y0,
           h = h0,
           maxerr = 0.0;
    for(int i=1;i<=n;i++){
        while(1){
            tuple<double, double, int> res;
            thread t([&res, y, x, h0, h, &f, &y_exact]() {
                res = RK2_maxerror(y, x, x+h0, h/2, f, y_exact);
            });
            auto [y1, err1, n1] = RK2_maxerror(y, x, x+h0, h, f, y_exact);
            t.join();
            auto [y2, err2, n2] = res;

            if(fabs(y2 - y1)/3.0 < eps){
                y = y2;
                maxerr = max(maxerr, err2);
                steps += n2;
                break;
            }
            h /= 2.0;
        }
        x = a + i*h0;
        h = h0;
    }
    return make_tuple(y, maxerr, steps);
}
