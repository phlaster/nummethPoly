#include "HEADER.hpp"
#include "Sparse.hpp"

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

spMtr normalize(const spMtr& M){
    if (M.cols != 1)
        throw invalid_argument("Sparse: only column matricies for normalization!");
    double norm = euclideanNorm(M);
    if (norm==0.0)
        throw invalid_argument("Can't norm zero sparse vector!");
    return M/norm;
}