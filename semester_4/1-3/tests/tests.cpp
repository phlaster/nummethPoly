#include <iostream>
#include <cassert>
#include <cmath>
#include "../functions.hpp"

using namespace std;

int main() {

    // Test scalarMultiply function
    assert(scalarMultiply(0.5, {3, 4, 5})[0] == 1.5);
    assert(scalarMultiply(0.5, {3, 4, 5})[1] == 2.0);
    assert(scalarMultiply(0.5, {3, 4, 5})[2] == 2.5);

    // Test multiplyMatrixVector function
	Matrix A = {{{-4,-7,-6, -3},
	             {1,-1, 1, 3},
	             {-2,-5, 2, -1},
                 {-3, 4,-8, -4}}};
	Vector x = {-5, 6, 0, 1};

	assert(multiplyMatrixVector(A,x)[0] == -25 );
	assert(multiplyMatrixVector(A,x)[1] == -8 );
	assert(multiplyMatrixVector(A,x)[2] == -21 );
    assert(multiplyMatrixVector(A,x)[3] == 35 );

	// Test vector sum function (vecSum)
	Vector a = {10 ,20 ,30}; double b=0.5; 
	Vector c ={40 ,50 ,60}; double d=-2;

	assert(vecSum(b,a,d,c)[0]==-75); 
    assert(vecSum(b,a,d,c)[1]==-90);    
    assert(vecSum(b,a,d,c)[2]==-105);

	//Test dotProduct function
	Vector u = {1.0, 2.0, -3.0};
    Vector v = {-4.5, 6.7, 8};

	assert(dotProduct(u,v) == -15.1 ); 

	//Test euclideanNorm function
	Vector p = {3,-4,12};

	assert(euclideanNorm(p)==13);
  
    cout << "All tests passed successfully!" << endl;

    return 0;
}