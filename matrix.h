/* Matthew Louigi Cabie Ong 2020 */
#ifndef MATRIX_H
#define MATRIX_H

#include "vector.h"

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

	/* Matrix Operations for Row and Columns */
	Matrix& rowSwap(int row1, int row2);
	Matrix& rowMult(int row, float k);
	Matrix& rowAdd(int row1, int row2, float k = 1);
	Matrix& colSwap(int col1, int col2);
	Matrix& colMult(int col, float k);
	Matrix& colAdd(int col1, int col2, float k = 1);

	Matrix& rowReduce();

	/* Getters and setters */
	float& get(int row, int col) const { return m_matrix[row][col]; }
	float at(int row, int col) const { return m_matrix[row][col]; }
	const int getRow() const { return m_rows; }
	const int getCol() const { return m_cols; }
	inline int dim() const { return m_rows * m_cols; }

	/* Operators */
	friend std::ostream& operator<<(std::ostream& out, const Matrix& mat);
	Matrix operator*(const float k);
	Matrix operator/(const float k);
	Matrix operator+(const Matrix& m);
	Matrix operator-(const Matrix& m);
	float& operator()(int i, int j);
};

#endif // !MATRIX_H
