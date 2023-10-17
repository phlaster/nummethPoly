#include "LinearAlgebra.hpp"

bool issquare(const Mtr& M){
    return M.size() == M[0].size();
}
bool isapprox(const Vec& v1, const Vec& v2, double tol) {
    if (v1.size() != v2.size())
        return false;
    for (size_t i = 0; i < v1.size(); ++i)
        if (fabs(v1[i] - v2[i]) > tol)
            return false;
    return true;
}
bool isapprox(const Mtr& m1, const Mtr& m2, double tol) {
    if (m1.size() != m2.size() || m1[0].size() != m2[0].size())
        return false;
    for (size_t i = 0; i < m1.size(); ++i)
        if (!isapprox(m1[i], m2[i], tol))
            return false;
    return true; 
}
bool issingleval(const Mtr& M){
    return issquare(M) && M.size() == 1;
}