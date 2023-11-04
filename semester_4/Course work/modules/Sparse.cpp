#include "headers/Sparse.hpp"

spMtr::spMtr( size_t r, size_t c) : rows(r), cols(c), zeroCounter(r*c), valCounter(0){
}

spMtr::spMtr( size_t r, size_t c, double sparsity):rows(r), cols(c), zeroCounter(r*c), valCounter(0){
    if (sparsity<0. || sparsity>1.)
        throw invalid_argument("Sparsity value s must be 0<=s<=1");


    for( size_t i = 0; i < rows; ++i ){
        for( size_t j = 0; j < cols; ++j ){
            if ((double)rand()/RAND_MAX > sparsity){
            }
            else {
                keypair key = make_pair(i, j);
                data[key] = (double)rand()/RAND_MAX;
                valCounter++; zeroCounter--;
            }
        }
    }
}

spMtr::spMtr(Mtr M, double eps){
    rows = M.size();
    cols = M[0].size();

    for (size_t i = 0; i < rows; ++i){
        for (size_t j = 0; j < cols; ++j){
            if (fabs(M[i][j]) <= eps){
                zeroCounter++;
            } else {
                keypair key = make_pair(i, j);
                data[key] = M[i][j];
                valCounter++;
            }
        }
    }
}

spMtr::~spMtr(){
}

double spMtr::sparsity(){
    return (double)zeroCounter/(rows*cols);
}

double spMtr::get(size_t rIndex, size_t cIndex) const {
    if( rIndex >= rows || cIndex >= cols )
        throw out_of_range("Indices out of range");
    keypair key = make_pair(rIndex, cIndex);
    return data.count(key) ? data.at(key) : 0.0;
}

Vec spMtr::get(size_t rIndex) const {
    if( rIndex >= rows )
        throw out_of_range("Row index out of range");
    Vec res(cols);
    for (size_t i=0; i<cols; i++){
        keypair key = make_pair(rIndex, i);
        res[i] = data.count(key) ? data.at(key) : 0.0;
    }
    return res;
}

void spMtr::set(double value, size_t rIndex, size_t cIndex){
    if( rIndex >= rows || cIndex >= cols )
        throw out_of_range("Indices out of range");
    
    keypair key = make_pair(rIndex, cIndex);
    if (data.count(key) && fabs(value) > 0.0){
        data[key] = value;
    } else if (data.count(key) && fabs(value) == 0.0){
        data.erase(key);
        zeroCounter++;
        valCounter--;
    } else if (!data.count(key) && fabs(value) > 0.0){
        data[key] = value;
        zeroCounter--;
        valCounter++;
    }
}