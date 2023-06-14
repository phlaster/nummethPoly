#include <iostream>
#include <cassert>
#include <cmath>
#include "../functions.hpp"

using namespace std;

int main() {

    assert(scalarMultiply(0.5, {3, 4, 5})[0] == 1.5);
    assert(scalarMultiply(0.5, {3, 4, 5})[1] == 2.0);
    assert(scalarMultiply(0.5, {3, 4, 5})[2] == 2.5);

	Matrix A = {{-4,-7,-6, -3},
	             {1,-1, 1, 3},
	             {-2,-5, 2, -1},
                 {-3, 4,-8, -4}};
	Vector x = {-5, 6, 0, 1};

	assert(multiplyMatrixVector(A,x)[0] == -25 );
	assert(multiplyMatrixVector(A,x)[1] == -8 );
	assert(multiplyMatrixVector(A,x)[2] == -21 );
    assert(multiplyMatrixVector(A,x)[3] == 35 );

	Vector v1 = {10 ,20 ,30}; double c1=0.5; 
	Vector v2 ={40 ,50 ,60}; double c2=-2;

	assert(vecSum(c1,v1,c2,v2)[0]==-75); 
    assert(vecSum(c1,v1,c2,v2)[1]==-90);    
    assert(vecSum(c1,v1,c2,v2)[2]==-105);

	Vector u = {1.0, 2.0, -3.0};
    Vector v = {-4.5, 6.7, 8};

	assert(dotProduct(u,v) == -15.1 ); 

	Vector p = {3,-4,12};

	assert(euclideanNorm(p)==13);
  
    cout << "Тесты успешно пройдены!" << endl;

    return 0;
}