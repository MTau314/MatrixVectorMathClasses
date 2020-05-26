#ifndef MVMATH_H
#define MVMATH_h

#include "matrix.h"

namespace VMath
{
	float magnitude(const Vector& v);
	float dot(const Vector& v, const Vector& u);
	Vector cross(const Vector& u, const Vector& v);
	Vector proj(const Vector& v, Vector& onto);
}

namespace MMath
{
	int cofRank(const Matrix& A, int m, int n);
	int augRank(const Matrix& A);
	double det(const Matrix& A);
	Matrix transpose(const Matrix& A);
	Matrix adj(const Matrix& A);
	Matrix inv(const Matrix& A);
	Matrix multiply(const Matrix& A, const Matrix& B);
	Matrix multiply(const Matrix& A, Vector& v);
	Matrix selfMultiply(const Matrix& A, int exp);
}

#endif // !MATH_H

