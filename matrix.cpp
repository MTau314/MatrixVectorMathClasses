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

	return Matrix{ temp };
}

void Matrix::setAll(float value)
{
	for (int i{}; i < m_rows; i++)
		for (int j{}; j < m_cols; j++)
			m_matrix[i][j] = (float)value;
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

/* operator */
std::ostream& operator<<(std::ostream& out, const Matrix& mat)
{
	for (int i{}; i < mat.m_rows; i++)
	{
		std::cout << std::setw(2);
		for (int j{}; j < mat.m_cols; j++)
		{
			if (std::abs(mat.m_matrix[i][j] - 0.0) < 1e-12)
				mat.m_matrix[i][j] = std::abs(mat.m_matrix[i][j]);
			out << mat.m_matrix[i][j] << std::setw(10);
		}

		std::cout << std::setw(1);
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

	return Matrix{ virtualProduct };
}

Matrix Matrix::operator/(const float k)
{
	return *this * (1 / k);
}

Matrix Matrix::operator+(const Matrix& m)
{
	Matrix virtualSum{ *this };
	for (int i{}; i < m_rows; i++)
		for (int j{}; j < m_rows; j++)
			virtualSum.m_matrix[i][j] += m.m_matrix[i][j];

	return Matrix{ virtualSum };
}

Matrix Matrix::operator-(const Matrix& m)
{
	Matrix virtualDiff{ *this };
	for (int i{}; i < m_rows; i++)
		for (int j{}; j < m_rows; j++)
			virtualDiff.m_matrix[i][j] -= m.m_matrix[i][j];

	return Matrix{ virtualDiff };
}

float& Matrix::operator()(int i, int j) { return m_matrix[i][j]; }
