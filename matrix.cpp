/* Matthew Louigi Cabie Ong 2020 */
#include "matrix.h"

/* ----- Vector functions ----- */
/* Vector Utility Functions */
void Vector::setAll(float value)
{
	for (vSize_t i{}; i < m_size; i++)
		m_vector[i] = value;
}

float Vector::magnitude() const
{
	float sum{};
	for (vSize_t i{}; i < m_size; i++)
		sum += (m_vector[i] * m_vector[i]);
	return sqrt(sum); //similar idea to pythagorean thm.
}

/* Vector Operators */
Vector Vector::operator-()
{
	return *this * -1;
}

std::ostream& operator<<(std::ostream& out, const Vector& v)
{
	for (vSize_t i{}; i < v.getSize(); i++)
		out << v[i] << std::endl;
	return out;
}

Vector operator*(const Vector& v, float k)
{
	Vector prodVec{ v };
	for (vSize_t i{}; i < prodVec.getSize(); i++)
		prodVec[i] *= k;
	return prodVec;
}
Vector operator*(float k, const Vector& v)
{
	return v * k;
}

Vector operator/(float k, const Vector& v)
{
	return v * (1 / k);
}
Vector operator/(const Vector& v, float k)
{
	return v * (1 / k);
}

Vector operator+(const Vector& lhs, const Vector& rhs)
{
	assert(lhs.getSize() == rhs.getSize() && "vectors do not match in size");

	Vector sumVec{ lhs };
	for (vSize_t i{}; i < lhs.getSize(); i++)
		sumVec[i] += rhs[i];
	return sumVec;
}

Vector operator-(const Vector& lhs, const Vector& rhs)
{
	assert(lhs.getSize() == rhs.getSize() && "vectors do not match in size");

	Vector diffVec{ lhs };
	for (vSize_t i{}; i < lhs.getSize(); i++)
		diffVec[i] -= rhs[i];
	return diffVec;
}

/* Vector Math Functions */
float dot(const Vector& v, const Vector& u)
{
	assert(v.getSize() == u.getSize() && "Not the same size");

	float dotSum{};
	for (vSize_t i{}; i < v.getSize(); i++)
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

Vector proj(const Vector& v, const Vector& onto)
{
	return onto * (dot(v, onto)) / (onto.magnitude() * onto.magnitude());
}

/* ----- Matrix functions ----- */
/* Utility functions */
Matrix Matrix::resize(mSize_t row, mSize_t col) const //resizes current matrix to row and col specified
{
	if (row == m_rows && col == m_cols)
		return *this;

	Matrix temp{ row, col };

	mSize_t rowIndex = (m_rows > row) ? row : m_rows;
	mSize_t colIndex = (m_cols > col) ? col : m_cols;

	for (mSize_t i{}; i < rowIndex; i++)
		for (mSize_t j{}; j < colIndex; j++)
			temp.m_matrix[i][j] = m_matrix[i][j];
	return temp;
}

Matrix Matrix::remove(mSize_t row, mSize_t col) const //removes a specific row and col
{
	Matrix temp{ m_rows - 1, m_cols - 1 };
	mSize_t tempRowIndex{};
	for (mSize_t i{}; i < m_rows && tempRowIndex < temp.m_rows; i++)
	{
		if (i == row)
			continue;

		mSize_t tempColIndex{};
		for (mSize_t j{}; j < m_cols && tempColIndex < temp.m_cols; j++)
		{
			if (j == col)
				j++;

			temp.m_matrix[tempRowIndex][tempColIndex] = m_matrix[i][j];
			tempColIndex++;
		}

		tempRowIndex++;
	}

	return temp;
}

void Matrix::setAll(float value)
{
	for (mSize_t i{}; i < m_rows; i++)
		for (mSize_t j{}; j < m_cols; j++)
			m_matrix[i][j] = value;
}

void Matrix::eye()
{
	assert(m_rows == m_cols && "Not square");
	for (mSize_t i{}; i < m_rows; i++)
		m_matrix[i][i] = 1.0f;
}

int Matrix::cofRank(mSize_t rowLimit, mSize_t colLimit) const
{
	//user specifies the size of the coefficient matrix to compute
	Matrix cofM{ this->resize(rowLimit, colLimit) }; //narrow to wanted size ignoring the vector matrix
	cofM.rowReduce();

	int rank{};
	for (mSize_t rlead{}; rlead < rowLimit; rlead++)
	{
		bool isNonZero{ true };
		for (mSize_t clead{}; clead < colLimit && isNonZero; clead++) //check entire row if non-zero
		{
			if (cofM.at(rlead, clead) != 0)
			{
				++rank;
				isNonZero = false;
			}
		}

	}
	
	return rank;
}

int Matrix::augRank() const
{		
	return cofRank(m_rows, m_cols);
}

double Matrix::det() const
{
	if (m_cols == 2 && m_rows == 2)
		return (m_matrix[0][0] * m_matrix[1][1] - m_matrix[0][1] * m_matrix[1][0]); // ad - bc
	if (m_cols == 1 && m_rows == 1)
		return m_matrix[0][0];

	double dDet{};
	for (mSize_t j{}; j < m_cols; j++)
	{
		Matrix D{ *this };
		dDet += D.remove(1, j).det() * static_cast<double>(m_matrix[1][j]) * pow(-1, 1 + j);
	}

	return dDet;
}

Matrix& Matrix::eApply(float (*fcn)(float))
{
	for (mSize_t i{}; i < m_rows; i++)
		for (mSize_t j{}; j < m_cols; j++)
			m_matrix[i][j] = fcn(m_matrix[i][j]);
	return *this;
}

void readMatrix(std::ifstream& fileIn, Matrix& matrix)
{
	assert(fileIn && "File could not be opened"); //failsafe if user does not check prior to function call
	matrix.setAll(0);
	for (mSize_t rowCounter{}; rowCounter < matrix.getRow(); rowCounter++)
		for (mSize_t colCounter{}; colCounter < matrix.getCol(); colCounter++)
			fileIn >> matrix.at(rowCounter, colCounter);
}

/* Matrix Operations for Rows and Cols */
Matrix& Matrix::rowSwap(int row1, int row2)
{
	for (mSize_t j{}; j < m_cols; j++)
	{
		float temp = m_matrix[row1][j];
		m_matrix[row1][j] = m_matrix[row2][j];
		m_matrix[row2][j] = temp;
	}

	return *this;
}

Matrix& Matrix::rowAdd(int reference, int destination, float k)
{
	for (mSize_t j{}; j < m_cols; j++)
		m_matrix[destination][j] += (m_matrix[reference][j] * k);

	return *this;
}

Matrix& Matrix::rowMult(int row, float k)
{
	for (mSize_t j{}; j < m_cols; j++)
		m_matrix[row][j] *= k;

	return *this;
}

Matrix& Matrix::colSwap(int col1, int col2)
{
	for (mSize_t i{}; i < m_rows; i++)
	{
		float temp = m_matrix[i][col1];
		m_matrix[i][col1] = m_matrix[i][col2];
		m_matrix[i][col2] = temp;
	}

	return *this;
}

Matrix& Matrix::colAdd(int reference, int destination, float k)
{
	for (mSize_t i{}; i < m_rows; i++)
		m_matrix[i][destination] += (m_matrix[i][reference] * k);

	return *this;
}

Matrix& Matrix::colMult(int col, float k)
{
	for (mSize_t i{}; i < m_rows; i++)
		m_matrix[i][col] *= k;

	return *this;
}

Matrix& Matrix::rowReduce() /*	courtesy of rosettacode/wikipedia pseudocode	*/
{
	mSize_t lead{ 0 }; //index of leading variable
	for (mSize_t row{}; row < m_rows; ++row) //will iterate through all rows of the matrix
	{
		if (lead >= m_cols) //will stop if reached the end of the matrix (col)
			return *this;

		mSize_t i{ row }; //cursor index
		while (m_matrix[i][lead] == 0) // looks for a non-zero leading variable
		{
			i++;
			if (i >= m_rows)  //resets if reached end of row
			{
				i = row;
				lead++;
				if (lead >= m_cols) //will stop if reached the end of the col
					return *this;
			}
		}

		rowSwap(i, row); //will swap for an available non-zero leading variable determined by the search algorithm
		rowMult(row, 1 / (m_matrix[row][lead])); //will make elements on pivot column into 1's
		for (mSize_t j{}; j < m_rows; ++j)  //will make elements on pivot column into 0's making it into a true pivot column
			if (j != row)
				rowAdd(row, j, -m_matrix[j][lead]);
	}

	return *this;
}

Matrix Matrix::transpose()
{
	Matrix T{ *this };
	for (mSize_t i{}; i < m_rows; i++)
		for (mSize_t j{}; j < m_cols; j++)
			T.m_matrix[j][i] = m_matrix[i][j]; //turn i-j elements into j-i elements
	return T;
}

/* Matrix Operators */
std::ostream& operator<<(std::ostream& out, const Matrix& mat)
{
	static const double TOL{ 1e-6 }, ZERO{ 0 };
	out << std::fixed << std::setprecision(4);
	for (mSize_t i{}; i < mat.getRow(); i++)
	{
		for (mSize_t j{}; j < mat.getCol(); j++)
		{
			float element{ mat.at(i,j) };
			if (std::fabs(mat.at(i, j) - ZERO) < TOL)
				element = std::fabs(mat.at(i, j));
			out << std::setw(10) << element << ' ';
		}
		out << std::endl;
	}
	return out;
}

// multiplication and division
Matrix operator*(const Matrix& m, float k)
{
	Matrix prodMat{ m };
	for (mSize_t i{}; i < m.getRow(); i++)
		for (mSize_t j{}; j < m.getCol(); j++)
			prodMat.at(i,j) *= k;
	return prodMat;
}
Matrix operator*(float k, const Matrix& m)
{
	return m * k;
}

Matrix operator*(const Matrix& lhs, const Matrix& rhs)
{
	assert(lhs.getCol() == rhs.getRow() && "Dimensions must match");
	Matrix ABProduct{ lhs.getRow(), rhs.getCol() }; //lhs_row by rhs_col matrix

	for (mSize_t i{}; i < lhs.getRow(); i++) //row cursor on product
	{
		mSize_t j{}; //col cursor on product
		for (mSize_t ii{}; ii < rhs.getCol(); ii++) //row cursor on BT(col cursor on B)
		{
			for (mSize_t jj{}; jj < rhs.getRow(); jj++) //col cursor on BT(row cursor on B)
			{
				ABProduct.at(i, j) += (lhs.at(i,jj) * rhs.at(jj, ii));
			}
			++j; //force to next product col
		}
		--j; //undo
	}
	return ABProduct;
}
Matrix operator*(const Matrix& m, const Vector& v)
{
	assert(m.getCol() == v.getSize() && "Dimensions must match");
	Matrix MVProduct{ m.getRow(), 1 }; //dimension will change to be rowCount by 1(colCount of vector)
	for (mSize_t i{}; i < m.getRow(); i++)
		for (mSize_t j{}; j < m.getCol(); j++)
			MVProduct(i, 0) += v.at(j) * m.at(i,j);
	return MVProduct;
}

Matrix operator/(const Matrix& m, float k)
{
	return m * (1 / k);
}
Matrix operator/(float k, const Matrix& m)
{
	return m * (1 / k);
}

// addition and subtraction
Matrix operator+(const Matrix& lhs, const Matrix& rhs)
{
	assert(lhs.getCol() == rhs.getCol() && lhs.getRow() == rhs.getRow() && "Matrices does not have matching dimensions");
	Matrix sumMat{ lhs };
	for (mSize_t i{}; i < lhs.getRow(); i++)
		for (mSize_t j{}; j < lhs.getCol(); j++)
			sumMat.at(i,j) += rhs.at(i,j);
	return sumMat;
}

Matrix operator-(const Matrix& lhs, const Matrix& rhs)
{
	assert(lhs.getCol() == rhs.getCol() && lhs.getRow() == rhs.getRow() && "Matrices does not have matching dimensions");
	Matrix diffMat{ lhs };
	for (mSize_t i{}; i < lhs.getRow(); i++)
		for (mSize_t j{}; j < lhs.getCol(); j++)
			diffMat.at(i, j) -= rhs.at(i, j);
	return diffMat;
}

// misc operators
const float& Matrix::operator()(int row, int col) const { return m_matrix[row][col]; }
float& Matrix::operator()(int row, int col) { return m_matrix[row][col]; }

const float* Matrix::operator[](int index) const { return m_matrix[index]; }
float* Matrix::operator[](int index) { return m_matrix[index]; }

Matrix Matrix::operator-()
{
	return *this * -1;
}

/* Matrix Math Functions */
Matrix adj(const Matrix& A)
{
	Matrix C{ A.getRow(), A.getCol() }; //cofactor matrix
	Matrix Acpy{ A };
	for (mSize_t i{}; i < A.getRow(); i++)
		for (mSize_t j{}; j < A.getCol(); j++)
		{
			Matrix D{ Acpy.remove(i,j) }; //matrix with the ith row and jth column removed per cofactor expansion
			C(i, j) = D.det() * pow(-1, i + j);
		}

	return C.transpose();
}

Matrix inv(const Matrix& A)
{
	assert(A.getRow() != A.getCol() || A.det() != 0 && "Inverse does not exist");
	return Matrix{ adj(A) * (1.0f / A.det()) };
}

Matrix power(const Matrix& A, int exp)
{
	Matrix AProduct{ A };
	for (mSize_t i{ 1 }; i < static_cast<mSize_t>(exp); i++)
		AProduct = AProduct * AProduct;
	return AProduct;
}