pair<int, double> normed_PM(const Mtr& A, double delta){
    size_t N = size(A), nsteps = 1;
    Vec y(N, 1), x = normalize(y), lambda(N, 0.0);
    for(;;) {
        Vec y = mul(A, x);
        Vec x_k = normalize(y);
        Vec lambda_k = div(y, x, delta);
        if (lambda_k.size() == 0){
            cerr << "All x_i in denominator are less than delta, return lambda from prev. step \n";
            return make_pair(nsteps, mean(lambda)); }
        if (fabs(mean(lambda) - mean(lambda_k)) <= delta)
            return make_pair(nsteps, mean(lambda_k));
        lambda = lambda_k;  x = x_k;
        nsteps++;
    }
}
