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

    Matrix operator*(Matrix &rhs);

    Matrix operator*(double d);

    Vector operator*(Vector &rhs);

    Matrix getTransparent();

    void setDimension(size_t n, size_t m);

    void setE();

    static Matrix getE(size_t n, size_t m);
};

#endif // SHARE_MATRIX_H