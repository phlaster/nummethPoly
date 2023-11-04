#include "headers/HEADER.hpp"

Vec randVec(size_t n, double lower, double upper){
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> distrib(lower, upper);
    
    Vec V(n);
    for (size_t i = 0; i < n; ++i)
        V[i] = distrib(gen);
    return V;
}


double euclideanNorm(const Vec& V){
    return sqrt(V*V);
}

double euclideanNorm(const spMtr& M){
    if (M.cols != 1)
        throw invalid_argument("Sparse: only column matricies for Eucledian norm!");

    Vec transposed = T(M).get(0);

    return euclideanNorm(transposed);
}

Vec normalize(const Vec& V){
    double norm = euclideanNorm(V);
    if (norm==0.0)
        throw invalid_argument("Can't norm zero vector!");
    return V/norm;
}

