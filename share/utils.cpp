#include "utils.h"

#include <iostream>
#include <cmath>

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
    return _A * _x - _b;
}

Vector Utils::solveSystemRotation(const Matrix &_A, const Vector &_b) {
    Matrix A = _A;
    Vector b = _b;

    size_t size = b.n;

    size_t size_d1 = size - 1;

    /* COMPUTE */
    //compute matrix[][] ROT method
    //forward

    for (int i = 0; i < size; ++i) {
        for (int j = i + 1; j < size; ++j) {
            double at = A[i][i]; //a11 ... (diagonal element)
            double bt = A[j][i]; //a21 ... (next row element)
            double c = at / sqrt(at * at + bt * bt);
            double s = bt / sqrt(at * at + bt * bt);
            for (int k = i; k < size; ++k) {
                double t1 = A[i][k]; //first row element
                double t2 = A[j][k]; //next row element
                A[i][k] = c * t1 + s * t2;
                A[j][k] = -s * t1 + c * t2;
            }

            double t1 = b[j]; //first row element
            double t2 = b[j]; //next row element
            b[j] = c * t1 + s * t2;
            b[j] = -s * t1 + c * t2;

        }
    }

    //compute x[]
    //backward

    Vector x;
    x.setDimension(size);

    for (size_t i = size_d1; i < size; --i) {
        double sum = 0;
        for (size_t j = i + 1; j < size; ++j)
            sum += A[i][j] * x[j];
        double A_ii_mul_x = b[i] - sum;
        x[i] = A_ii_mul_x / A[i][i];
    }

    return x;
}

double Utils::computeCondition(const Matrix &_A) {
    Matrix A_1 = _A.getInverse();
    return _A.getNorm() * A_1.getNorm();
}

void Utils::readMatrix(istream &from, Matrix &to) {
    for (size_t i = 0; i < to.m; ++i) {
        for (size_t j = 0; j < to.n; ++j) {
            from >> to[i][j];
        }
    }
}

void Utils::readVector(istream &from, Vector &to) {
    for (size_t i = 0; i < to.n; ++i) {
        from >> to[i];
    }
}
