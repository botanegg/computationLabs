#ifndef SHARE_MATRIX_H
#define SHARE_MATRIX_H

#include <vector>
#include "vector.h"

using namespace std;

class Matrix {
private:
    vector<vector<double> > _m;
public:
    size_t n, m;

    vector<double> &operator[](size_t idx);

    const vector<double> &operator[](size_t idx) const;

    Matrix operator*(const Matrix &rhs) const;

    Matrix operator+(const Matrix &rhs) const;

    Matrix operator*(const double d) const;

    Vector operator*(const Vector &rhs) const;

    Matrix getTransparent() const;

    void setDimension(size_t n, size_t m);

    void setE();

    static Matrix getE(size_t n, size_t m);

    Matrix getInverse() const;

    double getNorm() const;
};

#endif // SHARE_MATRIX_H
