#include "vector.h"

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
