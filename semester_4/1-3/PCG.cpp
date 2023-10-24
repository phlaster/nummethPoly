#include "Functions.hpp"

pair<Vec, int> conjugateGradientMethod(const Mtr& A, const Vec& b, double eps, int maxIter, bool verbose){
    Vec x_k = Vec(b.size(), 1);
    Vec r = vecSum(1, b, -1, multiplyMatrixVector(A, x_k));
    Vec p = r;
    int k = 0;

    if (verbose){
        cout << "Before iterations:\n";
        cout << "A:\n"; print(A);
        cout << "b = "; print(b);
        cout << "k = " << k << endl;
        cout << "x_i = "; print(x_k);
        cout << "r = "; print(r);
        cout << "p = "; print(p);
        cout << "\n\nIterations:\n\n";
    }

    while (k <= maxIter){
        Vec q = multiplyMatrixVector(A, p);
        double pq_denom = dotProduct(p, q);
        double alpha = dotProduct(r, p) / pq_denom; 
        x_k = vecSum(1, x_k, alpha, p);
        r = vecSum(1, r, -alpha, q);
        if (euclideanNorm(r) <= eps)
            return make_pair(x_k, k);
        
        
        double beta = dotProduct(r, q) / pq_denom;
        p = vecSum(1, r, -beta, p);

        if (verbose){
            cout << "k = " << k << endl;
            cout << "alpha = " << alpha << endl;
            cout << "beta = " << beta << endl;
            cout << "x_i = "; print(x_k);
            cout << "r = "; print(r);
            cout << "||r||= " << euclideanNorm(r) << endl;
            cout << "p = "; print(p);
            cout << "q = "; print(q);
            cout << endl;
        }
        k++;
    }
    return make_pair(x_k, -1);
}

void PCG(const Mtr& A, const Vec& b, int lowestDeg, int maxIter, bool verbose){
    cout << "eps,err,iters,n="<<b.size()<<"\n";
    for (int p = -1; p >= lowestDeg; p--){
        double eps = pow(10, p);
        auto [x, n_iters] = conjugateGradientMethod(A, b, eps, maxIter, verbose);
        if (n_iters == -1){
            cout << eps << " и далее " <<">"<< maxIter <<" итераций.\n";
            break;
        } else{
            Vec r = vecSum(1, b, -1, multiplyMatrixVector(A, x));
            cout << eps << "," << euclideanNorm(r) << "," << n_iters <<"\n";
        }
    }
}

pair<double, double> PCG_multi(int n_repeats, int mtr_size, double eps, int maxIter, double cond){
    vInt sums(n_repeats);
    double max_err = 0;
    for (int i=0; i<n_repeats; i++){

        Mtr A = cond < 1 ? generateRndSymPos(mtr_size) : generateRndSymPos(mtr_size, cond);
        Vec x = generateRandomVector(mtr_size);
        Vec b = multiplyMatrixVector(A, x);
        auto [_x_, n_iters] = conjugateGradientMethod(A, b, eps, maxIter, false);
        if (n_iters == -1){
            cerr << "Iteration limit exceeded!\n";
            break;
        } else{
            double err_norm = euclideanNorm(vecSum(1, x, -1, _x_));
            max_err = err_norm > max_err ? err_norm : max_err;
            sums[i] = n_iters;
        }
    }
    return make_pair(mean(sums), max_err);
}

void PCG_multi_convergence(int n_repeats, int maxIter, double cond){
    cout << "eps,n_10,err_10,n_30,err_30,n_50,err_50,n_100,err_100,n_200,err_200,n_500,err_500,nrepeats="<<n_repeats<<",cond="<<cond<<endl;
    for (double eps : {1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14,1e-15}){
        cout << eps << ",";
        for (int size : {10,30,50,100,200,500}){
            auto [mean_iters, max_err] = PCG_multi(n_repeats, size, eps, maxIter, cond);
            cout << mean_iters << "," << max_err << ",";
        }
        cout << endl;
    }
}
