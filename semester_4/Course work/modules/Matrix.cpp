#include "HEADER.hpp"
#include <stdexcept>

Vec randVec(int n, double lower, double upper){
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> distrib(lower, upper);
    
    Vec V(n);
    for (int i = 0; i < n; ++i)
        V[i] = distrib(gen);
    return V;
}


Mtr T(const Mtr& M){
    int rows = M.size();
    int cols = M[0].size();

    Mtr res(cols, Vec(rows));

    for (int i = 0; i < rows; ++i)
        for(int j = 0; j < cols; ++j)
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

Mtr generateRndSymPos(int n, double cond){
    Vec w0 = randVec(n);
    Mtr H = householder(w0);
    Mtr D = E(n);

    Vec diag_elems = randVec(n-1, 1, cond);
    for (int i=0; i<n-1; i++)
        D[i][i] = diag_elems[i];
    D[0][0] = cond;

    Mtr A = T(H) * D * H;
    return A;
}




Mtr to_dense(const spMtr& SP){
    int rows = SP.rows;
    int cols = SP.cols;
    Mtr res(rows, Vec(cols, 0.0));
    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            res[i][j] = SP.get(i,j);
        }
    }
    return res;
}







spMtr T(const spMtr& M){
    spMtr res(M.cols, M.rows);
    for (int i = 0; i < M.rows; ++i)
        for(int j = 0; j < M.cols; ++j){
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

void erase_above_diag(spMtr& A, bool below){
    if (A.cols != A.rows)
        throw invalid_argument("Eraser applies only to sqr mtrx!");
    
    int n = A.rows;
    for (int i = 0; i < n; ++i) {
        for (int j = i+1; j < n; ++j) {
            if (below) {
                A.set(0.0, j, i);
            } else {
                A.set(0.0, i, j);
            }
        }
    }
}

void chol(spMtr& A) {
    int n = A.rows;
    for (int i = 0; i < n; ++i) {
        double S = A.get(i, i);

        for (int ip = 0; ip < i; ++ip) {
            S -= A.get(i, ip) * A.get(i, ip);
        }
        double root = sqrt(S);
        A.set(root, i, i);        
        for (int j = i + 1; j < n; ++j) {
            double S = A.get(j,i);
            for (int ip = 0; ip < i; ++ip) {
                S -= A.get(i, ip) * A.get(j, ip);
            }
            double to_set = S / A.get(i, i);
            A.set(to_set, j, i);
        }
    }
    erase_above_diag(A);
}


void incomp_chol_zero_tol(spMtr& A, const double threshold) {
    int N = A.rows;
    
    for (int I = 0; I < N; ++I) {
        double S = A.get(I,I);

        for (int IP = 0; IP < I; ++IP) {
            if (A.get(I,IP) != 0.0 && A.get(IP,IP) != 0.0 &&  fabs(A.get(I, IP) / A.get(IP, IP)) > threshold) {
                S -= A.get(I,IP) * A.get(I,IP) / A.get(IP,IP);
            }
        }

        if (S <= 0) {
            A.set(0.0, I,I);
        } else {
            A.set(sqrt(S), I,I);
        }
        
        for (int J = I + 1; J < N; ++J) {
            double S = A.get(J,I);

            for (int IP = 0; IP < I; ++IP) {
                if (A.get(I,IP) != 0.0 && A.get(IP,IP) != 0.0 && fabs(A.get(J, IP) / A.get(IP, IP)) > threshold) {
                    S -= A.get(J,IP) * A.get(I,IP) / A.get(IP,IP);
                }
            }
            
            if (A.get(I,I) != 0.0) {
                A.set(S / A.get(I,I), J,I);
            } else {
                A.set(0.0, J,I);
            }
        }
    }
    erase_above_diag(A);
}



void ___incomp_chol_tol(spMtr& A, const double threshold) {
    int N = A.rows;
    for (int I = 0; I < N; I++) {
        double S = A.get(I, I);

        for (int IP = 0; IP < I; IP++) {
            S -= A.get(I, IP) * A.get(I, IP);
        }
        A.set(sqrt(S), I,I);
        for (int J = I + 1; J < N; J++) {
            double S = A.get(J, I);
            for (int IP = 0; IP < I; IP++) {
                S -= A.get(I, IP) * A.get(J, IP);
            }
            A.set(S / A.get(I, I), J, I);
            if (fabs(A.get(J, I)) < threshold) {
                A.set(0.0, J, I);
            }
        }
    }
    erase_above_diag(A);
}