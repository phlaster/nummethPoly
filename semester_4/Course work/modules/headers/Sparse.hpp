#ifndef SPARSE_HPP
#define SPARSE_HPP

#include "Types.hpp"

using keypair = pair<size_t,size_t>;

struct spMtr{
    size_t rows;
    size_t cols;
    size_t zeroCounter;
    size_t valCounter;
    map<keypair, double> data;


    spMtr(size_t r, size_t c);
    spMtr(size_t r, size_t c, double sparsity);
    spMtr(Mtr matrix, double eps=1e-16);
    ~spMtr();
    void print();
    double sparsity();
    void rotate90Clockwise();

    
    double get(size_t rIndex, size_t cIndex) const;
    Vec get(size_t rIndex) const;
    void set(double value, size_t rIndex, size_t cIndex);
};


#endif
