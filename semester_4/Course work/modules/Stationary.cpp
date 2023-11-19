#include "headers/HEADER.hpp"

spMtr block5diag(size_t m, size_t n) {
    size_t total_size = m * n;
    spMtr matrix(total_size, total_size);

    for (size_t block = 0; block < n; ++block) {
        size_t block_start = block * m;

        for (size_t row = 0; row < m; ++row) {
            size_t col = block_start + row;

            matrix.set(4, col, col);

            if (row > 0) { // Up
                matrix.set(-1, col - 1, col);
            }
            if (row < m - 1) { // Down
                matrix.set(-1, col + 1, col);
            }
            if (col > block_start) { // Left
                matrix.set(-1, col, col - 1);
            }
            if (col < block_start + m - 1) { // Right
                matrix.set(-1, col, col + 1);
            }
        }

        if (block < n - 1) {
            for (size_t row = 0; row < m; ++row) {
                size_t col = (block + 1) * m + row;
                matrix.set(-1, col - m, col);
            }
        }
        if (block > 0) {
            for (size_t row = 0; row < m; ++row) {
                size_t col = (block - 1) * m + row;
                matrix.set(-1, col + m, col);
            }
        }
    }
    return matrix;
}

spMtr block5diag(size_t n){
    return block5diag(n,n);
}

Vec estimate_b(size_t n) {
    const double PI = acos(-1);
    const double h = 1.0/(n+1);
    
    auto f = [PI](double x, double y) {
        return 2 * PI * PI * sin(PI * x) * cos(PI * y);
    };
    
    Vec b(n * n);
    int k = 0;
    for (size_t i = 1; i <= n; ++i) {
        for (size_t j = 1; j <= n; ++j) {
            double x = i * h;
            double y = j * h;
            b[k] = h * h * -f(x,y);

            if (j == 1) {
                b[k] -= sin(PI * x);
            }
            if (j == n) {
                b[k] += sin(PI * x);
            }
            k++;
        }
    }
    return b;
}

Vec exact_x(double n){
    const double PI = acos(-1);
    const double h = 1.0/(n+1);
    
    auto u = [PI](double x, double y) {
        return sin(PI * x) * cos(PI * y);
    };
    
    Vec exact(n * n);
    int k = 0;
    for (size_t i = 1; i <= n; ++i) {
        for (size_t j = n; j >= 1; j--) {
            double x = i * h;
            double y = j * h;

            exact[k] = u(x,y);
            k++;
        }
    }
    return exact;
}