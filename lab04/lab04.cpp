#include <iostream>
#include <fstream>
#include <share/matrix.h>
#include <share/utils.h>

using namespace std;

int main() {
    fstream file("input.txt");

    Matrix lab04 = Matrix::getE(4, 4);

    Utils::readMatrix(file, lab04);

    cout << "powerLambdaMethod lab04 is ";
    cout << endl;
    cout << Utils::powerLambdaMethod(lab04, 0.000001);
    cout << endl;

    cout << endl;

    cout << "rotationLambdaMethod lab04 is ";
    cout << endl;
    cout << Utils::rotationLambdaMethod(lab04, 0.001);
    cout << endl;

    return 0;
}
