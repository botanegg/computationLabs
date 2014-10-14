#include <iostream>
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

    Matrix A_Inverse = A.getInverse();

    cout << "matrix A[][]" << endl;
    Utils::printMatrix(A);
    cout << endl;

    cout << "matrix A^(-1)[][]" << endl;
    Utils::printMatrix(A_Inverse);
    cout << endl;

    Vector x = Utils::solveSystemLU(A, b);

    cout << "vector x[]" << endl;
    Utils::printVector(x);
    cout << endl;

    cout << "residual[]" << endl;
    Utils::printVector(Utils::computeResidual(A, x, b));
    cout << endl;

    return 0;
}
