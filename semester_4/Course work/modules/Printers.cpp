#include "headers/Printers.hpp"
#include <fstream>

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
            cout << setw(12) << elem << " ";
        cout << "\n";
    }
    cout << "\n";
}
void print(const Mtr& A, const Vec& x, const Vec& b){
    size_t n = A.size();
    
    cout << "A:" << setw(13*n+1) << "|x:" << setw(13) << "  |b:" << '\n';

	for(size_t i=0; i<n; ++i){
		for(size_t j=0; j<n; ++j)
			cout << setw(12) << A[i][j] << " ";
		cout << "|"<< setw(12) << x[i] << "|" << setw(12) << b[i] << "\n";
	}
    cout << "\n";
}
void print(const Mtr& A, const vector<string>& S, const Vec& b){
    size_t n = A.size();
    
    cout << "A:" << setw(13*n+1) << "|x:" << setw(13) << "  |b:" << '\n';

	for(size_t i=0; i<n; ++i){
		for(size_t j=0; j<n; ++j)
			cout << setw(12) << A[i][j] << " ";
		cout << "|"<< setw(12) << S[i] << "|" << setw(12) << b[i] << "\n";
	}
    cout << "\n";
}
void print(const double value) {
    cout << value << "\n\n";
}
void print(const bool b) {
    cout << (b ? "true" : "false") << "\n";
}
void print(const char* s) {
    cout << s << "\n";
}

void print(const spMtr& sA){
    for (size_t i = 0; i < sA.rows; i++){
        for (size_t j = 0; j < sA.cols; j++){
            double value = sA.get(i,j);
            if (value == 0.0){
                cout << "           . ";
            }
            else
                cout << setw(12) << value << " ";
        }
        cout << "\n";
    }
    cout << "\n";
}

void save(const spMtr& A, const string& fname){
    ofstream stream(fname, ofstream::trunc);
    for (int i=0; i<A.rows; i++){
        for (int j=0; j<A.cols; j++)
            stream << A.get(i, j) << " ";
        stream << "\n";
    }
}
