#ifndef SHARE_UTILS_H
#define SHARE_UTILS_H

#include "matrix.h"
#include "vector.h"

#include <istream>

class Utils {
public:
    static void printMatrix(const Matrix &m);

    static void printVector(const Vector &v);

    static Vector solveSystemLU(const Matrix &_A, const Vector &_b);

    static Vector solveSystemRotation(const Matrix &_A, const Vector &_b);

    static Vector computeResidual(const Matrix &_A, const Vector &_x, const Vector &_b);

    static double computeCondition(const Matrix &_A);

    static void readMatrix(istream &from, Matrix &to);

    static void readVector(istream &from, Vector &to);
};

#endif // SHARE_UTILS_H
