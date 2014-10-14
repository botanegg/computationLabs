#include "matrix.h"

std::vector<double> &Matrix::operator[](size_t idx) {
    return _m[idx];
}

Matrix Matrix::getTransparent() {
    Matrix nm;

    nm.setDimension(this->m, this->n);

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            nm[i][j] = _m[j][i];
        }
    }
    return nm;
}

void Matrix::setDimension(size_t n, size_t m) {
    this->m = m;
    this->n = n;

    _m.resize(m);

    for (size_t j = 0; j < m; ++j) {
        _m[j].resize(n);
    }
}

Matrix Matrix::operator*(Matrix &rhs) {
    Matrix nm;
    nm.setDimension(this->m, rhs.n);

    for (size_t i = 0; i < this->m; ++i) {
        for (size_t j = 0; j < rhs.n; ++j) {
            double sum = 0;
            for (size_t k = 0; k < rhs.m; ++k) {
                sum += this->_m[i][k] * rhs[k][j];
            }
            nm[i][j] = sum;
        }
    }

    return nm;
}

void Matrix::setE() {
    for (size_t i = 0; i < n; ++i) {
        _m[i][i] = 1;
    }
}

Matrix Matrix::getE(size_t n, size_t m) {
    Matrix tmp;
    tmp.setDimension(n, m);
    tmp.setE();

    return tmp;
}

Matrix Matrix::operator*(double d) {
    Matrix nm;
    nm.setDimension(this->n, this->m);

    for (size_t i = 0; i < this->m; ++i) {
        for (size_t j = 0; j < this->n; ++j) {
            nm[i][j] = this->_m[i][j] * d;
        }
    }

    return nm;
}

Vector Matrix::operator*(Vector &rhs) {
    Vector n;
    n.setDimension(this->n);

    for (size_t i = 0; i < this->m; ++i) {
        double sum = 0;
        for (size_t k = 0; k < rhs.n; ++k) {
            sum += this->_m[i][k] * rhs[k];
        }
        n[i] = sum;
    }

    return n;
}

const vector<double> &Matrix::operator[](size_t idx) const {
    return _m[idx];
}
