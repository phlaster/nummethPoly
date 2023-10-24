#include <cassert>
#include "Functions.hpp"

using namespace std;

int main() {

    assert(scalarMultiply(0.5, {3, 4, 5})[0] == 1.5);
    assert(scalarMultiply(0.5, {3, 4, 5})[1] == 2.0);
    assert(scalarMultiply(0.5, {3, 4, 5})[2] == 2.5);

	Mtr A = {
		{-4,-7,-6,-3},
	    { 1,-1, 1, 3},
		{-2,-5, 2,-1},
        {-3, 4,-8,-4}
	};
	Vec x = {-5, 6, 0, 1};

	assert(multiplyMatrixVector(A,x)[0] == -25 );
	assert(multiplyMatrixVector(A,x)[1] == -8 );
	assert(multiplyMatrixVector(A,x)[2] == -21 );
    assert(multiplyMatrixVector(A,x)[3] == 35 );

	Vec v1 = {10 ,20 ,30}; double c1=0.5; 
	Vec v2 ={40 ,50 ,60}; double c2=-2;

	assert(vecSum(c1,v1,c2,v2)[0]==-75); 
    assert(vecSum(c1,v1,c2,v2)[1]==-90);    
    assert(vecSum(c1,v1,c2,v2)[2]==-105);

	Vec u = {1.0, 2.0, -3.0};
    Vec v = {-4.5, 6.7, 8};

	assert(dotProduct(u,v) == -15.1 ); 

	Vec p = {3,-4,12};

	assert(euclideanNorm(p)==13);
  
    cout << "Тесты успешно пройдены!" << endl;


	Vec _x_ = generateRandomVector(4);
	print(_x_);

	Mtr RSP = generateRndSymPos(15, 10);
	print(RSP);

    return 0;
}