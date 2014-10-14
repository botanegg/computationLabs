#include "utils.h"

#include <iostream>

void Utils::printMatrix(const Matrix &m) {
    size_t mn_d1 = m.n - 1;
    for (int i = 0; i < m.m; ++i) {
        for (int j = 0; j < m.n; ++j) {
            cout << m[i][j];
            if (j != mn_d1) cout << ' ';
        }
        cout << endl;
    }
}

void Utils::printVector(const Vector &v) {
    size_t vn_d1 = v.n - 1;
    for (int i = 0; i < v.n; ++i) {
        cout << v[i];
        if (i != vn_d1) cout << ' ';
    }
    cout << endl;
}

Vector Utils::solveSystemLU(const Matrix &_A, const Vector &_b) {
    Matrix A = _A;
    Vector b = _b;

    size_t size = b.n;

    size_t size_d1 = size - 1;

    Matrix L, U;
    L.setDimension(size, size);
    U.setDimension(size, size);

    /* COMPUTE */
    //compute LU method
    //forward


    //compute L U matrix
    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < size; ++j) {
            double sum = 0;
            for (int k = 0; k < i; ++k) {
                sum += L[i][k] * U[k][j];
            }

            if (i <= j) {
                //compute U matrix
                U[i][j] = A[i][j] - sum;
            } else {
                //compute L matrix
                L[i][j] = (A[i][j] - sum) / U[j][j];
            }

            if (i == j) L[i][j] = 1;
        }
    }

    //Ly = b
    //Ux = y

    //compute y[]
    //forward L

    Vector y;
    y.setDimension(size);

    for (size_t i = 0; i < size; ++i) {
        double sum = 0;
        for (size_t j = 0; j < i; ++j)
            sum += L[i][j] * y[j];
        double L_ii_mul_x = b[i] - sum;
        y[i] = L_ii_mul_x;
    }

    //compute x[]
    //backward U

    Vector x;
    x.setDimension(size);

    for (size_t i = size_d1; i < size; --i) {
        double sum = 0;
        for (size_t j = i + 1; j < size; ++j)
            sum += U[i][j] * x[j];
        double U_ii_mul_x = y[i] - sum;
        x[i] = U_ii_mul_x / U[i][i];
    }
    return x;
}

Vector Utils::computeResidual(const Matrix &_A, const Vector &_x, const Vector &_b) {
    return _A*_x - _b;
}
