#include "Functions.hpp"


void print(const Vec& v){
    cout << "[";
    for (auto item : v)
        cout << item << ", ";
    cout << "\b\b]\n\n";
}

void print(const vInt& v){
    cout << "[";
    for (auto item : v)
        cout << item << ", ";
    cout << "\b\b]\n\n";
}

void print(const Mtr& A){
    for (auto row : A){
        for (auto elem : row)
            cout << setw(11) << elem << " ";
        cout << "\n";
    }
    cout << "\n";
}

void print(const Mtr& A, const Vec& x, const Vec& b){
    int n = A.size();
    
    cout << "A:" << setw(12*n+1) << "|x:" << setw(12) << "  |b:" << '\n';

	for(int i=0; i<n; ++i){
		for(int j=0; j<n; ++j)
			cout << setw(11) << A[i][j] << " ";
		cout << "|"<< setw(11) << x[i] << "|" << setw(11) << b[i] << "\n";
	}
    cout << "\n";
}