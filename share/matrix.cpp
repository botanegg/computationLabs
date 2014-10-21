#include "matrix.h"

#include <cmath>
#include <iostream>

std::vector<double> &Matrix::operator[](size_t idx) {
    return _m[idx];
}

Matrix Matrix::getTransparent() const {
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

Matrix Matrix::operator*(const Matrix &rhs) const {
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

Matrix Matrix::operator*(const double d) const {
    Matrix nm;
    nm.setDimension(this->n, this->m);

    for (size_t i = 0; i < this->m; ++i) {
        for (size_t j = 0; j < this->n; ++j) {
            nm[i][j] = this->_m[i][j] * d;
        }
    }

    return nm;
}

Vector Matrix::operator*(const Vector &rhs) const {
    Vector n = Vector::get0(this->n);

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

Matrix Matrix::getInverse() const {
    Matrix src(*this);
    Matrix mat = Matrix::getE(n, m);

    size_t size = this->n;

    /* COMPUTE */
    //compute matrix[][] ROT method
    //forward

    for (size_t i = 0; i < size; ++i) {
        for (size_t j = i + 1; j < size; ++j) {
            double a = src[i][i]; //a11 ... (diagonal element)
            double b = src[j][i]; //a21 ... (next row element)
            double c = a / sqrt(a * a + b * b);
            double s = b / sqrt(a * a + b * b);
            for (size_t k = i; k < size; ++k) {
                double t1 = src[i][k]; //first row element
                double t2 = src[j][k]; //next row element
                src[i][k] = c * t1 + s * t2;
                src[j][k] = -s * t1 + c * t2;

                double h1 = mat[i][k]; //first row element
                double h2 = mat[j][k]; //next row element
                mat[i][k] = c * h1 + s * h2;
                mat[j][k] = -s * h1 + c * h2;
            }
        }
    }

    for (size_t i = size - 1; i < size; --i) { //warning sign trick
        double d = 1.0 / src[i][i];

        for (size_t j = 0; j < size; ++j) {
            src[i][j] *= d;
            mat[i][j] *= d;
        }

        for (size_t j = i - 1; j < size; --j) {
            if (src[j][i] != 0) {
                double h = src[j][i];
                for (size_t k = 0; k < size; ++k) {
                    src[j][k] = src[j][k] - src[i][k] * h;
                    mat[j][k] = mat[j][k] - mat[i][k] * h;
                }
            }
        }
    }

    return mat;
}

double Matrix::getNorm() const {
    double max = 0;
    for (size_t i = 0; i < n; ++i) {
        double sum = 0;
        for (size_t j = 0; j < n; ++j) {
            sum += abs(_m[j][i]);
        }
        if (sum > max) max = sum;
    }
    return max;
}

Matrix Matrix::operator+(Matrix const &rhs) const {
    Matrix nm;
    nm.setDimension(this->n, this->m);

    for (size_t i = 0; i < this->m; ++i) {
        for (size_t j = 0; j < this->n; ++j) {
            nm[i][j] = this->_m[i][j] + rhs._m[i][j];
        }
    }

    return nm;
}
