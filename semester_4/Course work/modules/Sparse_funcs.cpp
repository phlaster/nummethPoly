#include "headers/HEADER.hpp"
#include "headers/Sparse.hpp"

spMtr normalize(const spMtr& M){
    if (M.cols != 1)
        throw invalid_argument("Sparse: only column matricies for normalization!");
    double norm = euclideanNorm(M);
    if (norm==0.0)
        throw invalid_argument("Can't norm zero sparse vector!");
    return M/norm;
}

pair<double, double> mean_std(const spMtr& A) {
    size_t n = A.rows;
    size_t m = A.cols;
    double sum = 0.0;
    double sum_sq = 0.0;

    for (auto const& pair : A.data) {
        double value = pair.second;
        sum += value;
        sum_sq += value * value;
    }
    double mean = sum / (n * m);
    double variance = (sum_sq - (sum * sum)/n/m)/(n*m-1);
    double stdv = sqrt(variance);
    return make_pair(mean, stdv);
}

spMtr T(const spMtr& M){
    spMtr res(M.cols, M.rows);
    for (size_t i = 0; i < M.rows; ++i)
        for(size_t j = 0; j < M.cols; ++j){
            double val = M.get(i, j);
            res.set(val, j, i);
        }
    return res;
}

spMtr E(size_t size, bool sparse){
    spMtr res(size, size);
    if (sparse){
        for (size_t i=0; i<size; i++)
            res.set(1.0, i, i);
        return res;
    }
    return res;
}


spMtr householder(const spMtr& M){
    if (M.cols != 1)
        throw invalid_argument("Only column matricies for sparse Householder!");

    spMtr W = normalize(M);
    spMtr H = E(M.rows, true) - 2 * W * T(W);
    return H;
}

spMtr generateRndSymPos(int n, double cond, double sparsity){
    spMtr w0 = spMtr(n, 1, sparsity);
    spMtr H = householder(w0);
    spMtr D = E(n, true);

    Vec diag_elems = randVec(n-1, 1, cond);
    for (int i=0; i<n-1; i++)
        D.set(diag_elems[i], i,i);
    D.set(cond, 0, 0);

    spMtr A = T(H) * D * H;
    return A;
}

spMtr erase_above_diag(spMtr A, bool below){
    if (A.cols != A.rows)
        throw invalid_argument("Eraser applies only to sqr mtrx!");
    
    size_t n = A.rows;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i+1; j < n; ++j) {
            if (below) {
                A.set(0.0, j, i);
            } else {
                A.set(0.0, i, j);
            }
        }
    }
    return A;
}

spMtr chol(spMtr A, const double threshold) {
    size_t n = A.rows;
    for (size_t i = 0; i < n; i++) {
        double S = A.get(i, i);
        for (size_t ip = 0; ip < i; ip++) {
            S -= A.get(i, ip) * A.get(i, ip);
        }

        A.set(sqrt(S), i, i);

        for (size_t j = i + 1; j < n; j++) {
            double S = A.get(j, i);
            for (size_t ip = 0; ip < i; ip++) {
                S -= A.get(i, ip) * A.get(j, ip);
            }
            double to_set = S / A.get(i, i);
            // if (fabs(to_set) > threshold) {
                A.set(to_set, j, i);
            // } else {
            //     A.set(0.0, j, i); 
            // }
        }
    }
    return erase_above_diag(A);
}

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

pair<Vec, int> pcg(const spMtr& A, const Vec& b, double eps, int maxIter){
    Vec x_k = Vec(b.size(), 1);
    Vec r = b - A*x_k;
    
    Vec p = r;
    int k = 0;
    while (k <= maxIter){
        Vec q = A * p;
        double pq_denom = p*q;
        double alpha = r*p / pq_denom; 
        x_k = x_k + alpha*p;
        r = r - alpha*q;
        if (euclideanNorm(r) <= eps)
            return make_pair(x_k, k);
        double beta = r*q / pq_denom;
        p = r - beta*p;
        k++;
    }
    return make_pair(x_k, -1);
}

// Overload for preconditioner
// pair<Vec, int> pcg(const spMtr& A, const Vec& b, const spMtr& C, double eps, int maxIter){
//     Vec x_k = Vec(b.size(), 0);
//     Vec r = b;
//     int k = 0;
//     spMtr M = C*T(C);
//     Vec p, z0, r0;
//     while (euclideanNorm(r) > eps){
//         auto [z_k, n_z] = pcg(M, r, eps, maxIter); // solving for preconditioner
//         k++;
//         if (k==1){
//             p = z_k;
//         } else {
//             double beta = r*z_k / (r0*z0);
//             p = z_k + beta*p;
//         }
//         double alpha = r*z_k / (A*p * p);
//         r0 = r;
//         z0 = z_k;
//         x_k = x_k + alpha*p;
//         r = r - alpha*A*p;
//     }
//     return make_pair(x_k, k);
// }

pair<Vec, int> pcg(const spMtr& A, const Vec& b, const spMtr& Chol, double eps, int maxIter){

    Vec x_k = Vec(b.size(), 1);
    Vec r = b - A*x_k;
    Vec z = solve_L_U(Chol, T(Chol), r);
    Vec p = z;
    int k = 0;
    while (k <= maxIter){
        k++;
        Vec q = A*p;
        double denom = p*q;
        double alpha = z*q/denom;
        x_k = x_k + alpha*p;
        Vec r_k = r - alpha*q;
        if (euclideanNorm(r_k) <= eps)
            break;
        Vec z_k = solve_L_U(Chol, T(Chol), r_k);
        double beta = r_k*z_k / (r*z);
        p = r_k + beta*p;

        z = z_k;
        r = r_k;
    }
    return make_pair(x_k, k);
}

