/* Matthew Louigi Cabie Ong 2020 */
#include "matrix.h"



/* Utility functions */
Matrix& Matrix::resize(int row, int col) //resizes current matrix to row and col specified
{
	if (row == m_rows && col == m_cols)
		return *this;

	Matrix temp{ row, col };
	temp.setAll(0);

	int rowIndex = (m_rows > row) ? row : m_rows;
	int colIndex = (m_cols > col) ? col : m_cols;

	for (int i{}; i < rowIndex; i++)
		for (int j{}; j < colIndex; j++)
			temp.m_matrix[i][j] = m_matrix[i][j];

	*this = temp;
	return *this;
}

Matrix Matrix::remove(int row, int col) //removes a specific row and col
{
	Matrix temp{ m_rows - 1, m_cols - 1 };

	int ii{};
	for (int i{}; i < m_rows && ii < temp.m_rows; i++)
	{
		int jj{};

		if (i == row)
			continue;

		for (int j{}; j < m_cols && jj < temp.m_cols; j++)
		{
			if (j == col)
				++j;

			temp.m_matrix[ii][jj] = m_matrix[i][j];
			++jj;
		}
		++ii;
	}

	return temp;
}

void Matrix::setAll(float value)
{
	for (int i{}; i < m_rows; i++)
		for (int j{}; j < m_cols; j++)
			m_matrix[i][j] = (float)value;
}

void Matrix::eye()
{
	assert(m_rows == m_cols && "Not square");
	for (int i{}; i < m_rows; i++)
		m_matrix[i][i] = 1;
}

int Matrix::cofRank(int m, int n) const
{
	Matrix cofM{ *this };
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

int Matrix::augRank() const
{
	Matrix augM{ *this };
	augM.rowReduce();

	int rank{};
	for (int rlead{}; rlead < m_rows; rlead++)
		for (int clead{}; clead < m_cols; clead++)
			if (augM.at(rlead, clead) != 0)
			{
				++rank;
				break;
			}

	return rank;
}

double Matrix::det() const
{
	if (m_cols == 2 && m_rows == 2)
		return (m_matrix[0][0] * m_matrix[1][1] - m_matrix[0][1] * m_matrix[1][0]); // ad - bc
	if (m_cols == 1 && m_rows == 1)
		return m_matrix[0][0];

	double dDet{};
	for (int i{}; i < 1; i++)
		for (int j{}; j < m_cols; j++)
		{
			Matrix D{ *this };
			dDet += D.remove(i, j).det() * m_matrix[i][j] * pow(-1, i + j);
		}

	return dDet;
}

Matrix& Matrix::eApply(float (*fcn)(float))
{
	for (int i{}; i < m_rows; i++)
		for (int j{}; j < m_cols; j++)
			m_matrix[i][j] = fcn(m_matrix[i][j]);
	return *this;
}

/* Matrix Operations for Rows and Cols */
Matrix& Matrix::rowSwap(int row1, int row2)
{
	for (int j{}; j < m_cols; j++)
	{
		float temp = m_matrix[row1][j];
		m_matrix[row1][j] = m_matrix[row2][j];
		m_matrix[row2][j] = temp;
	}

	return *this;
}

Matrix& Matrix::rowAdd(int row1, int row2, float k)
{
	for (int j{}; j < m_cols; j++)
		m_matrix[row2][j] += (m_matrix[row1][j] * k);

	return *this;
}

Matrix& Matrix::rowMult(int row, float k)
{
	for (int j{}; j < m_cols; j++)
		m_matrix[row][j] *= k;

	return *this;
}

Matrix& Matrix::colSwap(int col1, int col2)
{
	for (int i{}; i < m_rows; i++)
	{
		float temp = m_matrix[i][col1];
		m_matrix[i][col1] = m_matrix[i][col2];
		m_matrix[i][col2] = temp;
	}

	return *this;
}

Matrix& Matrix::colAdd(int col1, int col2, float k)
{
	for (int i{}; i < m_rows; i++)
		m_matrix[i][col2] += (m_matrix[i][col1] * k);

	return *this;
}

Matrix& Matrix::colMult(int col, float k)
{
	for (int i{}; i < m_rows; i++)
		m_matrix[i][col] *= k;

	return *this;
}

Matrix& Matrix::rowReduce() /*	courtesy of rosettacode/wikipedia pseudocode	*/
{
	int lead{ 0 }; //index of leading variable

	for (int row{}; row < this->getRow(); ++row) //will iterate through all rows of the matrix
	{
		if (lead > this->getCol() - 1) //will stop if reached the end of the matrix (col)
			return *this;

		int i{ row }; //cursor index

		while (this->at(i, lead) == 0) // looks for a non-zero leading variable
		{
			++i;
			if (i > this->getRow() - 1)  //resets if reached end of row
			{
				i = row;
				++lead;

				if (lead > this->getCol() - 1) //will stop if reached the end of the col
					return *this;
			}
		}

		this->rowSwap(i, row); //will swap for an available non-zero leading variable determined by the search algorithm
		this->rowMult(row, 1 / (this->at(row, lead))); //will make elements on pivot column into 1's

		for (int j{}; j < this->getRow(); ++j)  //will make elements on pivot column into 0's making it into a true pivot column
			if (j != row)
				this->rowAdd(row, j, -1 * this->at(j, lead));

	}

	return *this;
}

Matrix Matrix::transpose()
{
	Matrix T{ *this };
	for (int i{}; i < m_rows; i++)
		for (int j{}; j < m_cols; j++)
			T(j, i) = m_matrix[i][j]; //turn i-j elements into j-i elements

	return T;
}

/* operator */
std::ostream& operator<<(std::ostream& out, const Matrix& mat)
{
	static const double TOL{ 1e-12 }, ZERO{ 0 };
	out << std::fixed << std::setprecision(4);
	for (int i{}; i < mat.m_rows; i++)
	{
		for (int j{}; j < mat.m_cols; j++)
		{
			if (std::abs(mat.m_matrix[i][j] - ZERO) < TOL)
				mat.m_matrix[i][j] = std::abs(mat.m_matrix[i][j]);
			out << std::setw(10) << mat.m_matrix[i][j] << ' ';
		}
		out << std::endl;
	}

	return out;
}

Matrix Matrix::operator*(const float k)
{
	Matrix virtualProduct{ *this };
	for (int i{}; i < m_rows; i++)
		for (int j{}; j < m_rows; j++)
			virtualProduct.m_matrix[i][j] *= k;

	return virtualProduct;
}

Matrix Matrix::operator/(const float k)
{
	return (*this * (1 / k));
}

Matrix Matrix::operator+(const Matrix& m)
{
	Matrix virtualSum{ *this };
	for (int i{}; i < m_rows; i++)
		for (int j{}; j < m_rows; j++)
			virtualSum.m_matrix[i][j] += m.m_matrix[i][j];

	return virtualSum;
}

Matrix Matrix::operator-(const Matrix& m)
{
	Matrix virtualDiff{ *this };
	for (int i{}; i < m_rows; i++)
		for (int j{}; j < m_rows; j++)
			virtualDiff.m_matrix[i][j] -= m.m_matrix[i][j];

	return virtualDiff;
}

Matrix Matrix::operator-()
{
	return Matrix{*this * -1 };
}

float& Matrix::operator()(int i, int j) { return m_matrix[i][j]; }

/* library functions */
Matrix adj(const Matrix& A)
{
	Matrix C{ A.getRow(), A.getCol() }; //cofactor matrix
	Matrix Acpy{ A };
	for (int i{}; i < A.getRow(); i++)
		for (int j{}; j < A.getCol(); j++)
		{
			Matrix D{ Acpy.remove(i,j) }; //matrix with the ith row and jth column removed per cofactor expansion
			C.get(i, j) = D.det() * pow(-1, i + j);
		}

	return C.transpose();
}

Matrix inv(const Matrix& A)
{
	assert(A.getRow() != A.getCol() || A.det() != 0 && "Inverse does not exist");
	return Matrix{ adj(A) * (1.0f / A.det()) };
}

Matrix multiply(const Matrix& A, const Matrix& B) // mxn and nxp will result in an mxp matrix 
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
	return ABProduct;
}

Matrix multiply(const Matrix& A, const Vector& v)
{
	assert(A.getCol() == v.getSize() && "Dimensions must match");

	Matrix MVProduct{ A.getRow(), 1 }; //dimension will change to be rowCount by 1(colCount of vector)

	for (int i{}; i < A.getRow(); i++)
		for (int j{}; j < A.getCol(); j++)
			MVProduct(i, 0) += v.at(j) * A.at(i, j);

	return MVProduct;
}

Matrix power(const Matrix& A, int exp)
{
	Matrix AProduct{ A };

	for (int i{ 1 }; i < exp; i++)
		AProduct = multiply(AProduct, A);

	return AProduct;
}

/* Functions */
void Vector::setAll(int value)
{
	for (int i{}; i < m_size; i++)
		m_vector[i] = static_cast<float>(value);
}

float Vector::magnitude()
{
	float sum{};
	for (int i{}; i < m_size; i++)
		sum += (m_vector[i] * m_vector[i]);

	return sqrt(sum); //similar idea to pythagorean thm.
}

/* Operators */
std::ostream& operator<<(std::ostream& out, const Vector& v)
{
	for (int i{}; i < v.m_size; i++)
		out << v.m_vector[i] << std::endl;

	return out;
}

Vector Vector::operator*(const float& k)
{
	Vector virtualProduct{ *this };
	for (int i{}; i < m_size; i++)
		virtualProduct.m_vector[i] *= k;

	return virtualProduct;
}

Vector Vector::operator+(const Vector& v)
{
	Vector virtualSum{ *this };
	assert(m_size == v.m_size && "vectors do not match in size");

	for (int i{}; i < m_size; i++)
		virtualSum.m_vector[i] += v.m_vector[i];

	return virtualSum;
}

Vector Vector::operator-(const Vector& v)
{
	Vector virtualDiff{ *this };
	assert(m_size == v.m_size && "vectors do not match in size");

	for (int i{}; i < m_size; i++)
		virtualDiff.m_vector[i] -= v.m_vector[i];

	return virtualDiff;
}

float dot(const Vector& v, const Vector& u)
{
	assert(v.getSize() == u.getSize() && "Not the same size");

	float dotSum{};
	for (int i{}; i < v.getSize(); i++)
		dotSum += (v.at(i) * u.at(i));

	return dotSum;
}

Vector cross(const Vector& u, const Vector& v)
{
	assert(v.getSize() == 3 && u.getSize() == 3 && "Not a 3D Vector");

	return Vector{ (u.at(1) * v.at(2)) - (u.at(2) * v.at(1)),
				   (u.at(2) * v.at(0)) - (u.at(0) * v.at(2)),
				   (u.at(0) * v.at(1)) - (u.at(1) * v.at(0)) };
}

Vector proj(const Vector& v, Vector& onto)
{
	float multi{ (dot(v,onto)) / (onto.magnitude() * onto.magnitude()) };
	Vector temp{ onto * multi };
	return temp;
}