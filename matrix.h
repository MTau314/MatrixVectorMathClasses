#ifndef MATRIX_H
#define MATRIX_H

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <fstream>

using vSize_t = std::size_t;
class Vector
{
private:
	float* m_vector{ nullptr };
	vSize_t m_size{ 0 };

public:
	Vector() = delete;
	Vector(const std::initializer_list<float>& elementsList) : m_size{ elementsList.size() }
	{
		m_vector = new float[static_cast<int>(elementsList.size())];
		int i{ 0 };
		for (const auto& element : elementsList)
		{
			m_vector[i] = element;
			i++;
		}
	}

	Vector(const Vector& srcVector) : m_size{ srcVector.m_size }
	{
		m_vector = new float[srcVector.m_size];
		for (vSize_t i{}; i < srcVector.m_size; i++)
			m_vector[i] = srcVector.m_vector[i];
	}

	Vector& operator=(const Vector& srcVector)
	{
		if (this == &srcVector)
			return *this;

		if (m_vector)
			delete[] m_vector;

		m_size = srcVector.m_size;
		m_vector = new float[srcVector.m_size];
		for (vSize_t i{}; i < srcVector.m_size; i++)
			m_vector[i] = srcVector.m_vector[i];

		return *this;
	}

	~Vector()
	{
		delete[] m_vector;
	}

public:
	/* Functions */
	float& at(int index) { return m_vector[index]; }
	const float& at(int index) const { return m_vector[index]; }
	const vSize_t& getSize() const { return m_size; }
	void setAll(float value);
	float magnitude() const;

	/* Member Operators */
	float& operator[](int index) { return m_vector[index]; }
	const float& operator[](int index) const { return m_vector[index]; }
	Vector operator-();
};

/* Non-member Operators for Vectors */
std::ostream& operator<<(std::ostream& out, const Vector& v);
Vector operator*(const Vector& v, float k);
Vector operator*(float k, const Vector& v);
Vector operator/(const Vector& v, float k);
Vector operator/(float k, const Vector& v);
Vector operator+(const Vector& lhs, const Vector& rhs);
Vector operator-(const Vector& lhs, const Vector& rhs);

/* Math functions for Vectors */
float dot(const Vector& v, const Vector& u);
Vector cross(const Vector& u, const Vector& v);
Vector proj(const Vector& v, const Vector& onto);

using mSize_t = std::size_t;
class Matrix
{
private:
	mSize_t m_rows{ 0 };
	mSize_t m_cols{ 0 };
	float** m_matrix{ nullptr };

/* Move semantics */
public:
	Matrix(Matrix&& srcMat) noexcept
		: m_matrix{ srcMat.m_matrix }, m_rows{ srcMat.m_rows }, m_cols{ srcMat.m_cols }
	{
		srcMat.m_rows = 0;
		srcMat.m_cols = 0;
		srcMat.m_matrix = nullptr;
	}

	Matrix& operator=(Matrix&& srcMat) noexcept
	{
		if (&srcMat == this)
			return *this;

		m_rows = srcMat.m_rows;
		m_cols = srcMat.m_cols;

		//delete existing resources
		for (mSize_t i{}; i < m_rows; i++)
			delete[] m_matrix[i];
		delete[] m_matrix;
		
		m_matrix = srcMat.m_matrix;

		srcMat.m_matrix = nullptr;
		srcMat.m_rows = 0;
		srcMat.m_cols = 0;

		return *this;
	}

/* Copy semantics */
public:
	Matrix(const Matrix& srcMatrix) : Matrix(srcMatrix.m_rows, srcMatrix.m_cols)
	{
		for (mSize_t i{}; i < m_rows; i++)
			for (mSize_t j{}; j < m_cols; j++)
				m_matrix[i][j] = srcMatrix.m_matrix[i][j];
	}

	Matrix& operator=(const Matrix& srcMatrix)
	{
		if (this == &srcMatrix)
			return *this;
		if (m_matrix) //if data exists already, clear it
		{
			//deallocate array of columns
			for (mSize_t i{}; i < m_rows; i++)
				delete[] m_matrix[i];

			//deallocate array of rows
			delete[] m_matrix;
		}

		m_rows = srcMatrix.m_rows;
		m_cols = srcMatrix.m_cols;

		m_matrix = new float* [srcMatrix.m_rows];

		for (mSize_t i{}; i < srcMatrix.m_rows; i++)
			m_matrix[i] = new float[srcMatrix.m_cols];

		for (mSize_t i{}; i < srcMatrix.m_rows; i++)
			for (mSize_t j{}; j < srcMatrix.m_cols; j++)
				m_matrix[i][j] = srcMatrix.m_matrix[i][j];

		return *this;
	}

/* Main Functors */
public:
	Matrix() = delete;
	Matrix(mSize_t m, mSize_t n) : m_rows{ m }, m_cols{ n }
	{
		//create an array of pointers as indices of a row
		m_matrix = new float*[m];
		//each row has an array of indices for columns
		for (mSize_t i{}; i < m_rows; i++)
			m_matrix[i] = new float[n];
		setAll(0);
	}

	//constructor with vectors
	Matrix(const std::initializer_list<Vector>& list) : m_cols { list.size() }
	{
		for (const auto& v : list)
		{
			mSize_t size{ v.getSize() };
			if (size <= m_rows)
				break;
			else if (size > m_rows)
				m_rows = size;
		}

		m_matrix = new float*[m_rows];
		for (mSize_t i{}; i < m_rows; i++)
			m_matrix[i] = new float[m_cols];

		mSize_t j{};
		for (const auto& vector : list)
		{
			if (j == m_cols)
				break;
			for (mSize_t i{}; i < m_rows; i++)
				m_matrix[i][j] = vector[i];
			++j;
		}
	}

	//augmented matrix constructor
	Matrix(const Matrix& cofMatrix, const Vector& b) : Matrix(cofMatrix)
	{
		assert(m_rows == b.getSize() && "Sizes do not match");
		resize(m_rows, m_cols + 1);

		for (mSize_t row{}; row < m_rows; row++)
		{
			if (row >= b.getSize())
			{
				m_matrix[row][m_cols - 1] = 0;
				break;
			}
			m_matrix[row][m_cols - 1] = b.at(row);
		}
	}

	//assignment by list for mxn declarations
	Matrix& operator=(const std::initializer_list<float>& list)
	{
		mSize_t i{}, j{};
		for (const auto& element : list)
		{
			if (j == m_cols)
			{
				i++;
				j = 0;
			}

			if (i == m_rows)
				break;

			m_matrix[i][j] = element;
			j++;
		}

		return *this;
	}

	//destructor
	~Matrix()
	{
		if (m_matrix)
		{
			//deallocate array of columns
			for (mSize_t i{}; i < m_rows; i++)
				delete[] m_matrix[i];
			//deallocate array of rows
			delete[] m_matrix;
		}
	}

public:
	/* Utility */
	Matrix resize(mSize_t row, mSize_t col) const;
	Matrix remove(mSize_t row, mSize_t col) const;
	void setAll(float value);
	void eye();

	int cofRank(mSize_t rowLimit, mSize_t colLimit) const;
	int augRank() const;
	double det() const;
	Matrix& eApply(float (*fcn)(float));

	float& at(int row, int col) { return m_matrix[row][col]; }
	const float& at(int row, int col) const { return m_matrix[row][col]; }
	const mSize_t& getRow() const { return m_rows; }
	const mSize_t& getCol() const { return m_cols; }
	const mSize_t dim() const { return m_rows * m_cols; }

	/* Matrix Operations for Row and Columns */
	Matrix& rowSwap(int row1, int row2);
	Matrix& rowMult(int row, float k);
	Matrix& rowAdd(int reference, int destination, float k = 1.0f);
	Matrix& colSwap(int col1, int col2);
	Matrix& colMult(int col, float k);
	Matrix& colAdd(int reference, int destination, float k = 1.0f);
	Matrix& rowReduce();
	Matrix transpose();

	/* Operators */
	const float& operator()(int row, int col) const;
	float& operator()(int row, int col);
	float* operator[](int index);
	const float* operator[](int index) const;
	Matrix operator-();
};

/* Non-member Operators for Matrices */
std::ostream& operator<<(std::ostream& out, const Matrix& mat);
Matrix operator*(const Matrix& m, float k);
Matrix operator*(float k, const Matrix& m);

Matrix operator*(const Matrix& lhs, const Matrix& rhs);
Matrix operator*(const Matrix& m, const Vector& v);

Matrix operator/(const Matrix& m, float k);
Matrix operator/(float k, const Matrix& m);

Matrix operator+(const Matrix& lhs, const Matrix& rhs);
Matrix operator-(const Matrix& lhs, const Matrix& rhs);

/* Functions */
Matrix adj(const Matrix& A);
Matrix inv(const Matrix& A);
Matrix power(const Matrix& A, int exp);
void readMatrix(std::ifstream& fileIn, Matrix& matrix);

#endif // !MATRIX_H
