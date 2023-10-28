#include "LinearAlgebra.hpp"

// pair<int, double> naive_PM(const Mtr& A, double delta){
//     size_t N = size(A);
//     Vec x_k, x = fill(1, N);
//     Vec lambda_k, lambda = fill(0, N);

//     int nsteps = 1;
//     for(;;) {
//         x_k = mul(A, x);
//         lambda_k = div(x_k, x, delta);
//         if (lambda_k.size() == 0){ // Если не осталось ни одного x_i, вернуть lambda с предыдущего шага
//             cerr << "Все x_i в знаменателе меньше пороговых значений, используем \n";
//             return make_pair(nsteps, mean(lambda));
//         }
//         if (fabs(mean(lambda) - mean(lambda_k)) <= delta) // euclideanNorm(sum(lambda, lambda_k, 1, -1))
//             return make_pair(nsteps, mean(lambda_k));
//         x = x_k;
//         lambda = lambda_k;
//         nsteps++;
//     }
// }
pair<int, double> naive_PM(const Mtr& A, double delta){
    size_t N = size(A), nsteps = 1;
    Vec x(N, 1.0), lambda(N, 0.0);
    for(;;) {
        Vec x_k = mul(A, x);
        Vec lambda_k = div(x_k, x, delta);
        if (lambda_k.size() == 0){
            cerr << "All x_i in denominator are less than delta, return lambda from prev. step \n";
            return make_pair(nsteps, mean(lambda)); }
        if (fabs(mean(lambda) - mean(lambda_k)) <= delta)
            return make_pair(nsteps, mean(lambda_k));
        x = x_k;  lambda = lambda_k;
        nsteps++;
    }
}

pair<int, double> normed_PM(const Mtr& A, double delta){
    size_t N = size(A), nsteps = 1;
    Vec y(N, 1), x = normalize(y), lambda(N, 0.0);
    for(;;) {
        Vec y = mul(A, x);
        Vec x_k = normalize(y);
        Vec lambda_k = div(y, x, 0.1*delta);
        if (lambda_k.size() == 0){
            cerr << "All x_i in denominator are less than delta, return lambda from prev. step \n";
            return make_pair(nsteps, mean(lambda)); }
        if (fabs(mean(lambda) - mean(lambda_k)) <= delta)
            return make_pair(nsteps, mean(lambda_k));
        lambda = lambda_k;  x = x_k;
        nsteps++;
    }
}

// pair<int, double> normed_PM(const Mtr& A, double delta){
//     size_t N = size(A);
//     Vec y = fill(1, N);
//     Vec x_k, x = normalize(y);
//     Vec lambda_k, lambda = fill(0, N);

//     int nsteps = 1;
//     for(;;) {
//         y = mul(A, x);
//         x_k = normalize(y);

//         lambda_k = div(y, x, delta);
//         if (lambda_k.size() == 0){
//             cerr << "Все x_i в знаменателе меньше пороговых значений, используем lambda с предыдущего шага\n";
//             return make_pair(nsteps, mean(lambda));
//         }
        
//         if (fabs(mean(lambda) - mean(lambda_k)) <= delta)
//             return make_pair(nsteps, mean(lambda_k));
//         lambda = lambda_k;
//         x = x_k;
//         nsteps++;
//     }
// }