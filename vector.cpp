#include "vector.h"

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