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