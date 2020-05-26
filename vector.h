/* Matthew Louigi Cabie Ong 2020 */
#ifndef VECTOR_H
#define VECTOR_H

#include <cassert>
#include <cmath>
#include <initializer_list>
#include <iomanip>
#include <iostream>

class Vector
{
private:
	float* m_vector;
	int m_size;

public:
	Vector() = delete;
	Vector(std::initializer_list<float> elementsList) : m_size{ (int)elementsList.size() }
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

	/* Operators */
	friend std::ostream& operator<<(std::ostream& out, const Vector& v);
	Vector operator*(const float& k);
	Vector operator+(const Vector& v);
	Vector operator-(const Vector& v);
};

#endif // !VECTOR_H
