#ifndef SHARE_VECTOR_H
#define SHARE_VECTOR_H

#include <vector>

using namespace std;

class Vector {
private:
    vector<double> _vec;
public:
    size_t n;

    double &operator[](size_t idx);

    const double &operator[](size_t idx) const;

    Vector operator-(const Vector &rhs) const;

    Vector operator+(const Vector &rhs) const;

    Vector operator*(const double d) const;

    double operator*(const Vector &rhs) const;

    Vector operator/(const double d) const;

    void setDimension(size_t n);

    void set0();

    double length() const;

    static Vector get0(size_t n);

    double getNorm() const;
};

#endif // SHARE_VECTOR_H
