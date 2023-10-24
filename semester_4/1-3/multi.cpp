#include "Functions.hpp"
#include <cstdlib>


int main(int argc, char **argv) {
    const int MAXITERS = 5000;
    
    int n_repeats = atoi(argv[1]);
    double cond = atof(argv[2]);

    PCG_multi_convergence(n_repeats, MAXITERS, cond);
    return 0;
}