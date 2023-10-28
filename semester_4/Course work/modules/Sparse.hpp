#ifndef SPARSE_HPP
#define SPARSE_HPP

#include "Types.hpp"

struct spMtr{
    int rows;
    int cols;
    int zeroCounter;
    int valCounter;
    map<string,double> data;


    spMtr(int r, int c);
    spMtr(int r, int c, double sparsity);
    spMtr(Mtr matrix, double eps=1e-16);
    ~spMtr();
    void print();
    double sparsity();
    void rotate90Clockwise();

    
    double get(int rIndex, int cIndex) const;
    Vec get(int rIndex) const;
    void set(double value, int rIndex, int cIndex);
};

#endif
