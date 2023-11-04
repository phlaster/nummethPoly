#include "headers/HEADER.hpp"


Mtr T(const Mtr& M){
    size_t rows = M.size();
    size_t cols = M[0].size();

    Mtr res(cols, Vec(rows));

    for (size_t i = 0; i < rows; ++i)
        for(size_t j = 0; j < cols; ++j)
            res[j][i] = M[i][j];
    return res;
}
Mtr E(size_t size){
    Mtr m(size, Vec(size, 0.0));
    for (size_t i=0; i<size; i++)
        m[i][i] = 1.0;
    return m;
}
Mtr toCol(const Vec& v){
    Mtr res(size(v), {0.0});
    for (size_t i=0; i< size(v); i++)
        res[i][0] = v[i];
    return res;
}

Mtr toRow(const Vec& v){
    return Mtr(1, v);
}

Mtr householder(const Vec& V){
    size_t n = size(V);
    Vec W = normalize(V);
    Mtr H = E(n) - 2*toCol(W)*toRow(W);
    return H;
}

Mtr generateRndSymPos(size_t n, double cond){
    Vec w0 = randVec(n);
    Mtr H = householder(w0);
    Mtr D = E(n);

    Vec diag_elems = randVec(n-1, 1, cond);
    for (size_t i=0; i<n-1; i++)
        D[i][i] = diag_elems[i];
    D[0][0] = cond;

    Mtr A = T(H) * D * H;
    return A;
}

Mtr generateRndSymPos(const Vec& V, const Vec& conds){
    size_t n = V.size();
    Mtr H = householder(V);
    Mtr D = E(n);

    for (size_t i=0; i<n; i++)
        D[i][i] = conds[i];

    Mtr A = T(H) * D * H;
    return A;
}

Mtr to_dense(const spMtr& SP){
    size_t rows = SP.rows;
    size_t cols = SP.cols;
    Mtr res(rows, Vec(cols, 0.0));
    for (size_t i = 0; i < rows; i++){
        for (size_t j = 0; j < cols; j++){
            res[i][j] = SP.get(i,j);
        }
    }
    return res;
}

pair<double, double> mean_std(const Mtr& A) {
    size_t n = A.size();
    size_t m = A[0].size();
    double sum = 0.0;
    double sum_sq = 0.0;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            sum += A[i][j];
            sum_sq += A[i][j] * A[i][j];
        }
    }
    double mean = sum / (n * m);
    double variance = (sum_sq - (sum * sum)/n/m)/(n*m-1);
    double stdv = sqrt(variance);
    return make_pair(mean, stdv);
}

