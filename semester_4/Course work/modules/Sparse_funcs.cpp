#include "headers/HEADER.hpp"
#include "headers/Sparse.hpp"
#include <cmath>


void spy(const spMtr& sA){
    cout << "non-zero: " << sA.valCounter << endl;
    if (sA.rows>45){
        cerr << "Matrix is too wide to display" << endl;
        return;
    }
    for (size_t i = 0; i < sA.rows; i++){
        for (size_t j = 0; j < sA.cols; j++){
            double value = sA.get(i,j);
            if (value == 0.0){
                cout << "   ";
            }
            else
                cout << "███";
        }
        cout << "\n";
    }
}

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

double maximum(const spMtr& A){
    double maximum = -INFINITY;
    for (const auto& [key, val] : A.data){
        maximum = max(maximum, val);
    }
    return maximum;
}
double maxabs(const spMtr& A){
    double maximum = 0.0;
    for (const auto& [key, val] : A.data){
        maximum = max(maximum, fabs(val));
    }
    return maximum;
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

spMtr generateRndSymPos(int n, double cond, double density){
    spMtr w0 = spMtr(n, 1, 1);
    spMtr H = householder(w0);
    spMtr D = E(n, true);

    Vec diag_elems = randVec(n-1, 1, cond);
    for (int i=0; i<n-1; i++)
        D.set(diag_elems[i], i,i);
    D.set(cond, 0, 0);

    spMtr A = T(H) * D * H;
    return sparsen(A, 1-density);
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

spMtr sparsen(spMtr M, double prob) {
    if (prob < 0.0 || prob > 1.0) {
        cerr << "Error: Probability must be in the range [0, 1]." << endl;
        exit(1);
    }
    if (M.rows!=M.cols){
        cerr << "Square matricies only!" << endl;
        exit(1);
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> rdist(0.0, 1.0);

    for (size_t i = 0; i < M.rows; ++i) {
        for (size_t j = i + 1; j < M.cols; ++j) {
            if (rdist(gen) < prob) {
                M.set(0.0, i, j);
                M.set(0.0, j, i);
            }
        }
    }
    return M;
}


