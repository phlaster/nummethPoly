#ifndef OPERATORS_HPP
#define OPERATORS_HPP
#include "Types.hpp"
#include "Sparse.hpp"

Vec operator+(Vec a, const Vec& b);
Vec operator-(Vec a);
Vec operator-(Vec a, const Vec& b);
Vec operator*(const double c, Vec a);
// Dot product
double operator*(const Vec& a, const Vec& b);
Vec operator/(Vec a, const double c);



Mtr operator+(Mtr A, const Mtr& B);
Mtr operator-(Mtr A);
Mtr operator-(Mtr A, const Mtr& B);
Mtr operator*(const double c, Mtr A);
Vec operator*(const Mtr& A, const Vec& x);
// Mtr mul
Mtr operator*(const Mtr& A, const Mtr& B);
Mtr operator/(Mtr A, const double c);



spMtr operator+(spMtr A, const spMtr& B);
spMtr operator-(spMtr A);
spMtr operator-(const spMtr& A, const spMtr& B);
spMtr operator*(const double c, spMtr A);
Vec operator*(const spMtr& A, const Vec& x);
// Mtr mul
spMtr operator*(const spMtr& A, const spMtr& B);
spMtr operator/(spMtr A, const double c);

//
void operator+=(spMtr& A, const spMtr& B);



#endif