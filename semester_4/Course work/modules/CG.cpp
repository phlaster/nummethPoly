#include "headers/HEADER.hpp"

Vec solve_L(const spMtr& L, const Vec& b) {
    size_t n = b.size();
    Vec y(n, 0.0);
    for (size_t i = 0; i < n; i++) {
        double sum = 0.0;
        for (size_t j = 0; j < i; j++)
            sum += L.get(i, j) * y[j];
        y[i] = (b[i] - sum) / L.get(i, i);
    }
    return y;
}
Vec solve_U(const spMtr& U, const Vec& b) {
    size_t n = b.size();
    Vec x(n, 0.0);
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (size_t j = i + 1; j < n; j++)
            sum += U.get(i, j) * x[j];
        x[i] = (b[i] - sum) / U.get(i, i);
    }
    return x;
}
Vec solve_L_U(const spMtr& L, const spMtr& U, const Vec& b){
    Vec y = solve_L(L, b);
    Vec x = solve_U(U, y);
    return x;
}
Vec solve_L_U(const spMtr& L, const Vec& b){
    Vec y = solve_L(L, b);
    Vec x = solve_U(T(L), y);
    return x;
}


pair<Vec, int> pcg(const spMtr& A, const spMtr& L, const Vec& b, double eps, int maxIter) {
    Vec x_k = Vec(b.size(), 1);
    Vec r_k = b - A * x_k;
    Vec z_k = solve_L_U(L, r_k);
    Vec p_k = z_k;
    int k = 1;
    while (k < maxIter) {
        Vec q_k = A * p_k;
        double alpha_denom = p_k*q_k;
        double alpha = z_k*r_k / alpha_denom;
        x_k = x_k + alpha*p_k;
        Vec r_k_new = r_k - alpha*q_k;
        if (euclideanNorm(r_k_new) <= eps) {
            return make_pair(x_k, k);
        }
        Vec z_k_new = solve_L_U(L, r_k_new);
        
        double beta = (z_k_new * r_k_new) / (z_k * r_k);
        p_k = z_k_new + beta * p_k;
        
        r_k = r_k_new;
        z_k = z_k_new;
        ++k;
    }
    cerr << "Метод не сошёлся за " << maxIter << " шагов!" << endl;
    return make_pair(x_k, -1);
}
pair<Vec, int> cg(const spMtr& A, const Vec& b, double eps, int maxIter) {
    Vec x_k = Vec(b.size(), 1);
    Vec r_k = b - A * x_k;
    Vec p_k = r_k;
    int k = 1;
    while (k < maxIter) {
        Vec q_k = A * p_k;
        double alpha_denom = p_k*q_k;
        double alpha = r_k*r_k / alpha_denom;
        x_k = x_k + alpha*p_k;
        Vec r_k_new = r_k - alpha*q_k;
        if (euclideanNorm(r_k_new) <= eps) {
            return make_pair(x_k, k);
        }        
        double beta = (r_k_new * r_k_new) / (r_k * r_k);
        p_k = r_k_new + beta * p_k;
        
        r_k = r_k_new;
        ++k;
    }
    cerr << "Метод не сошёлся за " << maxIter << " шагов!" << endl;
    return make_pair(x_k, -1);
}