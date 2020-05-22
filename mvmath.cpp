#include "mvmath.h"

/* Vector Math functions */
float VMath::magnitude(const Vector& v)
{
	float sum{};
	for (int i{}; i < v.getSize(); i++)
		sum += (v.at(i) * v.at(i));

	return sqrt(sum); //similar idea to pythagorean thm.
}

float VMath::dot(const Vector& v, const Vector& u)
{
	assert(v.getSize() == u.getSize() && "Not the same size");

	float dotSum{};
	for (int i{}; i < v.getSize(); i++)
		dotSum += (v.at(i) * u.at(i));

	return dotSum;
}

Vector VMath::cross(const Vector& u, const Vector& v)
{
	assert(v.getSize() == 3 && u.getSize() == 3 && "Not a 3D Vector");

	return Vector{ (u.at(1) * v.at(2)) - (u.at(2) * v.at(1)),
				   (u.at(2) * v.at(0)) - (u.at(0) * v.at(2)),
				   (u.at(0) * v.at(1)) - (u.at(1) * v.at(0)) };
}

Vector VMath::proj(const Vector& v, Vector& onto)
{
	float multi{ (dot(v,onto)) / (magnitude(onto) * magnitude(onto)) };
	Vector temp{ onto * multi };
	return Vector{ temp };
}

/* Matrix Math functions */

int MMath::cofRank(const Matrix& A, int m, int n)
{
	Matrix cofM{ A };
	cofM.resize(m, n); //narrow to wanted size ignoring the most-right hand side
	cofM.rowReduce();

	int rank{};
	for (int rlead{}; rlead < m; rlead++)
		for (int clead{}; clead < n; clead++) //checl entire row if non-zero
			if (cofM.at(rlead, clead) != 0)
			{
				++rank;
				break;
			}

	return rank;
}

int MMath::augRank(const Matrix& A)
{
	Matrix augM{ A };
	augM.rowReduce();

	int rank{};
	for (int rlead{}; rlead < A.getRow(); rlead++)
		for (int clead{}; clead < A.getCol(); clead++)
			if (augM.at(rlead, clead) != 0)
			{
				++rank;
				break;
			}

	return rank;
}

double MMath::det(const Matrix& A)
{
	if (A.getCol() == 2 && A.getRow() == 2)
		return (A.at(0, 0) * A.at(1, 1) - A.at(0, 1) * A.at(1, 0)); // ad - bc
	if (A.getCol() == 1 && A.getRow() == 1)
		return A.at(0, 0);

	double dDet{};

	for (int i{}; i < 1; i++)
		for (int j{}; j < A.getCol(); j++)
		{
			Matrix D{ A };
			dDet += det(D.remove(i, j)) * A.at(i, j) * pow(-1, i + j);
		}

	return dDet;
}

Matrix MMath::transpose(const Matrix& A)
{
	Matrix T{ A.getCol(), A.getRow() };
	for (int i{}; i < A.getRow(); i++)
		for (int j{}; j < A.getCol(); j++)
			T.get(j, i) = A.get(i, j); //turn i-j elements into j-i elements

	return Matrix{ T };
}

Matrix MMath::adj(const Matrix& A)
{
	Matrix C{ A.getRow(), A.getCol() }; //cofactor matrix
	Matrix Acpy{ A };
	for (int i{}; i < A.getRow(); i++)
		for (int j{}; j < A.getCol(); j++)
		{
			Matrix D{ Acpy.remove(i,j) }; //matrix with the ith row and jth column removed per cofactor expansion
			C.get(i, j) = det(D) * pow(-1, i + j);
		}

	return Matrix{ transpose(C) };
}

Matrix MMath::inv(const Matrix& A)
{
	assert(A.getRow() != A.getCol() || det(A) != 0 && "Inverse does not exist");
	return Matrix{ adj(A) * (1.0f / det(A)) };
}

Matrix MMath::multiply(const Matrix& A, const Matrix& B) // mxn and nxp will result in an mxp matrix 
{
	assert(A.getCol() == B.getRow() && "Dimensions must match");
	Matrix ABProduct{ A.getRow(), B.getCol() }; //mxp matrix

	for (int i{}; i < A.getRow(); i++) //row cursor on product
	{
		int j{}; //col cursor on product
		for (int ii{}; ii < B.getCol(); ii++) //row cursor on BT(col cursor on B)
		{
			for (int jj{}; jj < B.getRow(); jj++) //col cursor on BT(row cursor on B)
			{
				ABProduct.get(i, j) += (A.at(i, jj) * B.at(jj, ii));
			}
			++j; //force to next product col
		}
		--j; //undo
	}
	return Matrix{ ABProduct };
}

Matrix MMath::multiply(const Matrix& A, Vector& v)
{
	assert(A.getCol() == v.getSize() && "Dimensions must match");

	Matrix MVProduct{ A.getRow(), 1 }; //dimension will change to be rowCount by 1(colCount of vector)

	for (int i{}; i < A.getRow(); i++)
		for (int j{}; j < A.getCol(); j++)
			MVProduct(i, 0) += v.at(j) * A.at(i, j);

	return Matrix{ MVProduct };
}

Matrix MMath::selfMultiply(const Matrix& A, int exp)
{
	Matrix AProduct{ A };

	for (int i{ 1 }; i < exp; i++)
		AProduct = multiply(AProduct, A);

	return Matrix{ AProduct };
}