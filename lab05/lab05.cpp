#include <iostream>
#include <fstream>
#include <cmath>

#include <share/derivable.h>
#include <share/utils.h>


using namespace std;


Derivable linear_f(double x) {
    Derivable xd = Derivable::IndependendVariable(x);
    return cos(xd * xd) + 10 * xd;
};


Derivable system_f1x(double x) {
    Derivable xd = Derivable::IndependendVariable(x);
    return sin(xd + 0.5);
}

Derivable system_f1y(double y) {
    Derivable yd = Derivable::IndependendVariable(y);
    return ((-1)*yd) - 1.2;
}


Derivable system_f2x(double x) {
    Derivable xd = Derivable::IndependendVariable(x);
    return xd;
}

Derivable system_f2y(double y) {
    Derivable yd = Derivable::IndependendVariable(y);
    return cos(yd - 2);
}


void linear() {

    cout << "###" << endl << "Linear" << endl;
    double x = 10;
    double EPS = 0.0000000001;
    while (1) {
        Derivable thisF = linear_f(x);
        if (fabs(thisF.getValue()) < EPS) break;
        x = x - thisF.getValue() / thisF.getDerivative();
        std::cout << x << std::endl;
    }

    cout << "###" << endl << endl;
}


void system() {
    cout << "###" << endl << "System" << endl;
    Vector xk = Vector::get0(2);
    Vector fk = Vector::get0(2);
    double EPS = 0.0000000001;
    while (1) {

        fk[0] = system_f1x(xk[0]).getValue() + system_f1y(xk[1]).getValue();
        fk[1] = system_f2x(xk[0]).getValue() + system_f2y(xk[1]).getValue();


        Matrix A = Matrix::getE(2, 2);
        A[0][0] = system_f1x(xk[0]).getDerivative();
        A[0][1] = system_f1y(xk[1]).getDerivative();
        A[1][0] = system_f2x(xk[0]).getDerivative();
        A[1][1] = system_f2y(xk[1]).getDerivative();

        Vector dx = Utils::solveSystemLU(A, fk);
        xk = xk - dx;

        if (dx.getNorm() < EPS) break;

        Utils::printVector(xk);
    }

    cout << "FK = "; Utils::printVector(fk);
    cout << "XK = "; Utils::printVector(xk);

    cout << "###" << endl << endl;
}

int main() {
    linear();
    system();

    return 0;
}
