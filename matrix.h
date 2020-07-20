#ifndef MATRIX_H
#define MATRIX_H

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <fstream>

class Vector
{
private:
	float* m_vector;
	int m_size;

public:
	Vector() = delete;
	Vector(const std::initializer_list<float>& elementsList) : m_size{ (int)elementsList.size() }
	{
		m_vector = new float[(int)elementsList.size()];

		int i{ 0 };
		for (const auto element : elementsList)
		{
			m_vector[i] = element;
			++i;
		}
	}

	Vector(Vector& srcVector) : m_size{ srcVector.m_size }
	{
		m_vector = new float[srcVector.m_size];

		for (int i{}; i < srcVector.m_size; i++)
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

		for (int i{}; i < srcVector.m_size; i++)
			m_vector[i] = srcVector.m_vector[i];

		return *this;
	}

	~Vector()
	{
		delete[] m_vector;
	}

public:
	/* Functions */
	float& get(int index) { return m_vector[index]; }
	float at(int index) const { return m_vector[index]; }
	const int getSize() const { return m_size; }
	void setAll(int value);
	float magnitude();

	/* Operators */
	friend std::ostream& operator<<(std::ostream& out, const Vector& v);
	Vector operator*(const float& k);
	Vector operator+(const Vector& v);
	Vector operator-(const Vector& v);
};

float dot(const Vector& v, const Vector& u);
Vector cross(const Vector& u, const Vector& v);
Vector proj(const Vector& v, Vector& onto);

class Matrix
{
private:
	int m_rows;
	int m_cols;
	float** m_matrix;

public:
	Matrix() = delete;
	Matrix(int m, int n) : m_rows{ m }, m_cols{ n }
	{
		//create an array of pointers as indices of a row
		m_matrix = new float* [m];

		//each row has an array of indices for columns
		for (int i{}; i < m_rows; i++)
			m_matrix[i] = new float[n];
		this->setAll(0);
	}

	//constructor with vectors
	Matrix(const std::initializer_list<Vector>& list) : m_cols { static_cast<int>(list.size()) }
	{
		for (const auto& v : list)
		{
			int size{ v.getSize() };
			if (size <= m_rows)
				break;
			else if (size > m_rows)
				m_rows = size;
		}

		m_matrix = new float* [m_rows];
		for (int i{}; i < m_rows; i++)
			m_matrix[i] = new float[m_cols];

		int j{};
		for (const auto& vector : list)
		{
			if (j == m_cols)
				break;

			for (int i{}; i < m_rows; i++)
				m_matrix[i][j] = vector.at(i);
			++j;
		}
	}

	//copy constructor
	Matrix(const Matrix& srcMatrix) : m_rows{ srcMatrix.m_rows }, m_cols{ srcMatrix.m_cols }
	{
		m_matrix = new float* [m_rows];

		for (int i{}; i < m_rows; i++)
			m_matrix[i] = new float[m_cols];

		for (int i{}; i < m_rows; i++)
			for (int j{}; j < m_cols; j++)
				m_matrix[i][j] = srcMatrix.m_matrix[i][j];
	}

	//augmented matrix constructor
	Matrix(const Matrix& cofMatrix, const Vector& b) : Matrix(cofMatrix)
	{
		assert(m_rows == b.getSize() && "Sizes do not match");

		this->resize(m_rows, m_cols + 1);

		for (int row{}; row < m_rows; row++)
		{
			if (row >= b.getSize())
			{
				m_matrix[row][m_cols - 1] = 0;
				break;
			}

			m_matrix[row][m_cols - 1] = b.at(row);
		}

	}

	//assignment operator
	Matrix& operator=(const Matrix& srcMatrix)
	{
		if (this == &srcMatrix)
			return *this;

		if (m_matrix) //if data exists already, clear it
		{
			//deallocate array of columns
			for (int i{}; i < m_rows; i++)
				delete[] m_matrix[i];

			//deallocate array of rows
			delete[] m_matrix;
		}

		m_rows = srcMatrix.m_rows;
		m_cols = srcMatrix.m_cols;

		m_matrix = new float* [srcMatrix.m_rows];

		for (int i{}; i < srcMatrix.m_rows; i++)
			m_matrix[i] = new float[srcMatrix.m_cols];

		for (int i{}; i < srcMatrix.m_rows; i++)
			for (int j{}; j < srcMatrix.m_cols; j++)
				m_matrix[i][j] = srcMatrix.m_matrix[i][j];

		return *this;
	}

	//assignment by list for mxn declarations
	Matrix& operator=(std::initializer_list<float> list)
	{
		int i{}, j{};
		for (const auto element : list)
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
		//deallocate array of columns
		for (int i{}; i < m_rows; i++)
			delete[] m_matrix[i];

		//deallocate array of rows
		delete[] m_matrix;
	}

public:
	/* Utility */
	Matrix& resize(int row, int col);
	Matrix remove(int row, int col);
	void setAll(float value);
	void eye();

	int cofRank(int m, int n) const;
	int augRank() const;
	double det() const;

	Matrix& eApply(float (*fcn)(float));

	/* Matrix Operations for Row and Columns */
	Matrix& rowSwap(int row1, int row2);
	Matrix& rowMult(int row, float k);
	Matrix& rowAdd(int row1, int row2, float k = 1);
	Matrix& colSwap(int col1, int col2);
	Matrix& colMult(int col, float k);
	Matrix& colAdd(int col1, int col2, float k = 1);

	Matrix& rowReduce();
	Matrix transpose();

	/* Getters and setters */
	float& get(int row, int col) const { return m_matrix[row][col]; }
	float at(int row, int col) const { return m_matrix[row][col]; }
	const int getRow() const { return m_rows; }
	const int getCol() const { return m_cols; }
	int dim() const { return m_rows * m_cols; }
	Vector size() { return Vector{ float(m_rows), float(m_cols) }; }

	/* Operators */
	friend std::ostream& operator<<(std::ostream& out, const Matrix& mat);
	Matrix operator*(const float k);
	Matrix operator/(const float k);
	Matrix operator+(const Matrix& m);
	Matrix operator-(const Matrix& m);
	Matrix operator-();
	float& operator()(int i, int j);
};

Matrix adj(const Matrix& A);
Matrix inv(const Matrix& A);
Matrix multiply(const Matrix& A, const Matrix& B);
Matrix multiply(const Matrix& A, const Vector& v);
Matrix power(const Matrix& A, int exp);

#endif // !MATRIX_H
