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
