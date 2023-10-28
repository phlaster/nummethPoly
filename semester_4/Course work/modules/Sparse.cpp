#include "Sparse.hpp"
#include <algorithm>
#include <cstdlib>
#include <stdexcept>


spMtr::spMtr( int r, int c) : rows(r), cols(c), zeroCounter(r*c), valCounter(0){
}

spMtr::spMtr( int r, int c, double sparsity):rows(r), cols(c), zeroCounter(r*c), valCounter(0){
    if (sparsity<0. || sparsity>1.)
        throw invalid_argument("Sparsity value s must be 0<=s<=1");


    for( int i = 0; i < rows; ++i ){
        for( int j = 0; j < cols; ++j ){
            if ((double)rand()/RAND_MAX > sparsity){
            }
            else {
                string key = to_string(i) + "," + to_string(j);
                data[key] = (double)rand()/RAND_MAX;
                valCounter++; zeroCounter--;
            }
        }
    }
}

spMtr::spMtr(Mtr M, double eps){
    rows = M.size();
    cols = M[0].size();

    for (int i = 0; i < rows; ++i){
        for (int j = 0; j < cols; ++j){
            if (fabs(M[i][j]) <= eps){
                zeroCounter++;
            } else {
                string key = to_string(i) + "," + to_string(j);
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

void spMtr::rotate90Clockwise(){
    map<string,double> temp;
    string oldKey;
    string newKey;

    //Transpose matrix first
    for(int r = 0; r < rows; r++)
    {
      for(int c = r; c < cols; c++)
      {
        oldKey = to_string(r) + "," + to_string(c);
        newKey = to_string(c) + "," + to_string(r);

        temp[newKey] = data[oldKey];

        if( oldKey != newKey && data.count(newKey))
            temp[oldKey] = data[newKey];
      }
    }

    // Assign the temp map to our data map
    data = temp;

    // Matrix is transposed now so rows and cols should be swaped
    swap(rows, cols);

    // Reverse elements of matrix on row order
    for(int r = 0; r < rows; r++)
    {
        for(int c =0; c < cols/2; c++)
        {
            oldKey = to_string(r) + "," + to_string(c);
            newKey = to_string(r) + "," + to_string(cols-c-1);

            data[newKey] = temp[oldKey];

            if( oldKey != newKey && temp.count(newKey))
                data[oldKey] = temp[newKey];
        }
    }
}

double spMtr::get(int rIndex, int cIndex) const {
    if( rIndex >= rows || cIndex >= cols )
        throw out_of_range("Indices out of range");
    string key = to_string(rIndex) + "," + to_string(cIndex);
    return data.count(key) ? data.at(key) : 0.0;
}

Vec spMtr::get(int rIndex) const {
    if( rIndex >= rows )
        throw out_of_range("Row index out of range");
    Vec res(cols);
    for (int i=0; i<cols; i++){
        string key = to_string(rIndex) + "," + to_string(i);
        res[i] = data.count(key) ? data.at(key) : 0.0;
    }
    return res;
}

void spMtr::set(double value, int rIndex, int cIndex){
    if( rIndex >= rows || cIndex >= cols )
        throw out_of_range("Indices out of range");
    
    string key = to_string(rIndex) + "," + to_string(cIndex);
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