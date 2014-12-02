#include "vector.h"

#include <cmath>

double &Vector::operator[](size_t idx) {
    return _vec[idx];
}

void Vector::setDimension(size_t n) {
    this->n = n;
    _vec.resize(n);
}

void Vector::set0() {
    for (size_t i = 0; i < this->n; ++i) {
        _vec[i] = 0;
    }
}

Vector Vector::get0(size_t n) {
    Vector tmp;
    tmp.setDimension(n);
    tmp.set0();
    return tmp;
}

const double &Vector::operator[](size_t idx) const {
    return _vec[idx];
}

Vector Vector::operator-(Vector const &rhs) const {
    Vector tmp;
    tmp.setDimension(this->n);
    for (size_t i = 0; i < this->n; ++i) {
        tmp[i] = this->_vec[i] - rhs[i];
    }
    return tmp;
}

double Vector::getNorm() const {
    double sum = 0;
    for (int i = 0; i < n; ++i) {
        sum += abs(_vec[i]);
    }
    return sum;
}

Vector Vector::operator*(const double d) const {
    Vector tmp;
    tmp.setDimension(this->n);
    for (size_t i = 0; i < this->n; ++i) tmp[i] = _vec[i] * d;
    return tmp;
}

Vector Vector::operator+(const Vector &rhs) const {
    Vector tmp;
    tmp.setDimension(this->n);
    for (size_t i = 0; i < this->n; ++i) {
        tmp[i] = this->_vec[i] + rhs[i];
    }
    return tmp;
}

double Vector::length() const {
    double res = 0;
    for (size_t i = 0; i < this->n; ++i) {
        res += _vec[i] * _vec[i];
    }
    return sqrt(res);
}

Vector Vector::operator/(const double d) const {
    return (*this) * (1 / d);
}

double Vector::operator*(const Vector &rhs) const {
    double res = 0;
    for (size_t i = 0; i < this->n; ++i) res += _vec[i] * rhs[i];
    return res;
}
