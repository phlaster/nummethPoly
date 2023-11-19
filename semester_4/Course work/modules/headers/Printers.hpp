#ifndef PRINTERS_HPP
#define PRINTERS_HPP

#include "Sparse.hpp"

void print(const Vec& v);
void print(const vInt& v);
void print(const Mtr& A);
void print(const Mtr& A, const Vec& x, const Vec& b);
void print(const Mtr& A, const vector<string>& S, const Vec& b);
void print(const double value);
void print(const bool b);
void print(const char* s);
void print(const spMtr& sA);
void save(const spMtr& sA, const string& fname);

#endif