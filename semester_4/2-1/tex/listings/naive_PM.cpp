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
