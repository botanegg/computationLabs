#include "utils.h"

#include <iostream>
#include <cmath>

void Utils::printMatrix(const Matrix &m) {
    size_t mn_d1 = m.n - 1;
    for (size_t i = 0; i < m.m; ++i) {
        for (size_t j = 0; j < m.n; ++j) {
            cout << m[i][j];
            if (j != mn_d1) cout << ' ';
        }
        cout << endl;
    }
}

void Utils::printVector(const Vector &v) {
    size_t vn_d1 = v.n - 1;
    for (size_t i = 0; i < v.n; ++i) {
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
            for (size_t k = 0; k < i; ++k) {
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

    for (size_t i = 0; i < size; ++i) {
        for (size_t j = i + 1; j < size; ++j) {
            double at = A[i][i]; //a11 ... (diagonal element)
            double bt = A[j][i]; //a21 ... (next row element)
            double c = at / sqrt(at * at + bt * bt);
            double s = bt / sqrt(at * at + bt * bt);
            for (size_t k = i; k < size; ++k) {
                double t1 = A[i][k]; //first row element
                double t2 = A[j][k]; //next row element
                A[i][k] = c * t1 + s * t2;
                A[j][k] = -s * t1 + c * t2;
            }

            double t1 = b[i]; //first row element
            double t2 = b[j]; //next row element
            b[i] = c * t1 + s * t2;
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

Vector Utils::solveSOR(const Matrix &_A, const Vector &_b, double w) {
    Vector x = Vector::get0(_A.n);
    Vector x_n = x;
    size_t step = 0;

    do {
        step++;
        x = x_n;
        for (size_t j = 0; j < _A.m; ++j) {
            x_n[j] = w * _b[j] / _A[j][j] + (1 - w) * x[j];

            for (size_t k = 0; k < j; ++k) {
                x_n[j] -= w * _A[j][k] / _A[j][j] * x_n[k];

            }

            for (size_t k = j + 1; k < _A.n; ++k) {
                x_n[j] -= w * _A[j][k] / _A[j][j] * x[k];
            }
        }
    }
    while ((x - x_n).getNorm() > 0.00001 && step < ITTR_MAX_STEPS);

    return x_n;
}

Vector Utils::solveGZ(Matrix const &_A, Vector const &_b) {
    return solveSOR(_A, _b, 1.0);
}

Vector Utils::solveProgon(Matrix const &_A, Vector const &_b) {
    size_t n = _A.n;
    Vector P = Vector::get0(n);
    Vector Q = Vector::get0(n);
    P[0] = _A[0][1] / _A[0][0];
    Q[0] = -_b[0] / _A[0][0];

    for (size_t i = 1; i < n; i++) {
        auto yy = (i + 1) == n ? 0 : _A[i][i + 1];
        P[i] = yy / (-_A[i][i] - _A[i][i - 1] * P[i - 1]);
        Q[i] = (_A[i][i - 1] * Q[i - 1] - _b[i]) / (-_A[i][i] - _A[i][i - 1] * P[i - 1]);
    }

    Vector X = Vector::get0(n);
    X[n - 1] = Q[n - 1];
    for (size_t i = n - 1; i > 1; i--) {
        X[i-1] = X[i] * P[i-1] + Q[i-1];
    }

    return X;
}

double Utils::powerLambdaMethod(const Matrix &_A, double eps) {
    size_t n = _A.n;
    Vector x;
    x.setDimension(n);
    for (size_t j = 0; j < n; ++j) {
        x[j] = 1 / sqrt((double)n);
    }

    Vector y = _A * x;


    double prevLambda;
    double nextLambda = x * y;

    x = y / y.length();

    size_t step = 0;

    do {
        step++;
        prevLambda = nextLambda;
        y = _A * x;
        nextLambda = x * y;
        x = y / y.length();
    } while (abs(nextLambda - prevLambda) > eps && step < ITTR_MAX_STEPS);

    return nextLambda;
}

double Utils::rotationLambdaMethod(const Matrix &_A, double eps) {
    size_t n = _A.n;
    Matrix M = _A;

    double error = 1;
    size_t step = 0;
    while (error > eps && step < ITTR_MAX_STEPS) {
        step++;
        size_t mi = 0, mj = 1;
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j)
                if (i != j) if (abs(M[i][j]) > abs(M[mi][mj])) {
                    mi = i;
                    mj = j;
                }

        double k = sqrt(1 - (4 * M[mi][mj] * M[mi][mj]) / ((M[mi][mi] - M[mj][mj]) * (M[mi][mi] - M[mj][mj]) + (4 * M[mi][mj] * M[mi][mj])));
        double c = sqrt((1 + k) / 2);
        double s = sqrt((1 - k) / 2);
        if (fabs(M[mi][mi] - M[mj][mj]) < _EPS) {
            if (M[mi][mj] < 0) s = -s;
        }
        else if ((M[mi][mi] - M[mj][mj]) / M[mi][mj] < 0) s = -s;
        cout << "mi " << mi << " mj " << mj << endl;
        Matrix u = Matrix::getE(n, n);
        u[mi][mi] = u[mj][mj] = c;
        u[mi][mj] = -s;
        u[mj][mi] = s;
        if (step < 10) printMatrix(u);
        Matrix uT = u.getInverse();
        M = (uT * M) * u;
        if (step < 10) printMatrix(M);

        error = 0;
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j)
                if (i != j) error += abs(M[i][j]);
        cout << endl << error << endl;
        cout << endl;
    }

    printMatrix(M);
	return 0; //!TODO FIX IT
}
