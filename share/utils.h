#ifndef SHARE_UTILS_H
#define SHARE_UTILS_H

#include "matrix.h"
#include "vector.h"

class Utils {
public:
    static void printMatrix(const Matrix &m);

    static void printVector(const Vector &v);
};

#endif // SHARE_UTILS_H