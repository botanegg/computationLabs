#include <iostream>
#include <fstream>
#include <share/vector.h>
#include <share/matrix.h>
#include <share/utils.h>

using namespace std;

int main() {
    size_t size = 5;

    Matrix A = Matrix::getE(size, size);

    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < size; ++j) {
            A[i][j] = 1.0 / (i + j + 1);
        }
    }

    Vector fake_x = Vector::get0(size);
    fake_x[0] = 5;
    Vector b = A * fake_x;

    cout << "vector b[]" << endl;
    Utils::printVector(b);
    cout << endl;

    Matrix A_Inverse = A.getInverse();

    cout << "matrix A[][]" << endl;
    Utils::printMatrix(A);
    cout << endl;

    cout << "matrix A^(-1)[][]" << endl;
    Utils::printMatrix(A_Inverse);
    cout << endl;

    Vector x_LU = Utils::solveSystemLU(A, b);

    cout << "vector x[] LU method" << endl;
    Utils::printVector(x_LU);
    cout << endl;

    cout << "residual[] LU method" << endl;
    Utils::printVector(Utils::computeResidual(A, x_LU, b));
    cout << endl;

    Vector x_Rotation = Utils::solveSystemRotation(A, b);

    cout << "vector x[] Rotation method" << endl;
    Utils::printVector(x_Rotation);
    cout << endl;

    cout << "residual[] Rotation method" << endl;
    Utils::printVector(Utils::computeResidual(A, x_Rotation, b));
    cout << endl;

    cout << "Condition num is " << Utils::computeCondition(A);
    cout << endl;

    fstream file("input.txt");

    Matrix lab01 = Matrix::getE(4, 4);

    Utils::readMatrix(file, lab01);

    Utils::printMatrix(lab01);

    cout << "Condition num lab01 is " << Utils::computeCondition(lab01);
    cout << endl;

    return 0;
}
