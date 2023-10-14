// Вариант 11 Б
/*
    f1(x) = ctg(x) + x^2
    f2(x) = x^5 - 3.2x^3 + 2.5x^2 - 7x + 1.5
*/

#include "functions.hpp"

const Vec LIMS_1 = {0.5, 2.75};
const Vec LIMS_2 = {-2.4, 2.1};
const double rand_val = 0.9718075809237104;

void error_progression(
    int nodes_min,
    int nodes_max,
    double rand
){
    double a_1 = LIMS_1[0], b_1 = LIMS_1[1]; 
    double rand_x_1 = a_1 + rand*(b_1-a_1)/2;

    double a_2 = LIMS_2[0], b_2 = LIMS_2[1]; 
    double rand_x_2 = a_2 + rand*(b_2-a_2)/2;

    cout << "nNodes,err_1,err_2\n";
    for (int nNodes=nodes_min; nNodes<=nodes_max; nNodes++){
        auto [y_val_1, err_1] = lagrange_uniform_single_value_with_error(f1, rand_x_1, LIMS_1, nNodes);
        auto [y_val_2, err_2] = lagrange_uniform_single_value_with_error(f2, rand_x_2, LIMS_2, nNodes);
        cout << nNodes << "," << err_1 << "," << err_2 << endl;
    }
}

int main()
{
    error_progression(4, 100, rand_val);
    return 0;
}
