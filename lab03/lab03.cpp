#include <iostream>
#include <fstream>
#include <share/vector.h>
#include <share/matrix.h>
#include <share/utils.h>

using namespace std;

int main() {
    fstream file("input.txt");

    Matrix lab03 = Matrix::getE(4, 4);

    Utils::readMatrix(file, lab03);

    Utils::printMatrix(lab03);
    Vector xx = Vector::get0(4);
    xx[0] = 5;

    Vector bb = lab03 * xx;
    Utils::printVector(bb);
    cout << endl;
    Utils::printVector(Utils::solveSOR(lab03, bb));
    cout << endl;
    Utils::printVector(Utils::solveGZ(lab03, bb));
    cout << endl;
    cout << "Condition num lab03 is " << Utils::computeCondition(lab03);
    cout << endl;

    return 0;
}
