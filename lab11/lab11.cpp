#include <iostream>
#include <iomanip>
#include <cmath>

#include "share/matrix.h"

using namespace std;

double f(double x, double t, double b) {
    return (2 * (b * b - 1)) / ((x + b * t + 1) * (x + b * t + 1) * (x + b * t + 1));
}

double u0(double x) {
    return 1 / (1 + x);
}

double ur(double t, double b) {
    return 1 / (2 + b * t);
}

double u1(double t, double b) {
    return 1 / (b * t + 1);
}

double atx(double x, double t) {
    return 1;
}

double diff(double x, double b) {
    return -b / ((1 + x) * (1 + x));
}

double real(double x, double t, double b) {
    return 1 / (1 + x + b * t);
}

/*double f(double x, double y,double alpha)
{
	return (2 / (1 + y)) * (1 + ((x+alpha) / (1 + y)) * ((x+alpha) / (1 + y)) );
}*/

void solution(double x, double x0, double T, double T0, double h, double tau, double b) {
    size_t n = (size_t) ((T - T0) / tau + 1);
    size_t m = (size_t) ((x - x0) / h + 1);
    double res = 0;

    Matrix solutionMatrix = Matrix::getE(n, m);
    for (size_t i = 0; i < m; i++) {
        solutionMatrix[0][i] = u0(x0 + i * h);
        solutionMatrix[1][i] = tau * diff(x0 + i * h, b) + u0(x0 + i * h);
    }

    for (size_t i = 2; i < n; i++) {
        solutionMatrix[i][0] = u1(tau * i, b);
        for (size_t j = 1; j < m - 1; j++) {
            res = (tau * tau * atx(x0 + h * j, tau * i)) / (h * h);
            solutionMatrix[i][j] = res * solutionMatrix[i - 1][j + 1] + 2 * (1 - res) * solutionMatrix[i - 1][j] +
                                   res * solutionMatrix[i - 1][j - 1] - solutionMatrix[i - 2][j] +
                                   tau * tau * f(x0 + h * j, tau * i, b);
        }
        solutionMatrix[i][m - 1] = ur(tau * i, b);
    }
    for (size_t i = 0; i < 10; i++) {
        for (size_t j = 0; j < 10; j++)
            cout << setw(5) << solutionMatrix[i][j] << " ";
        cout << endl;
    }
}

int main() {
    double x0 = 0;
    double x = 1;
    double T0 = 0;
    double T = 1;
    double h = 0.1;
    double tau = 0.01;
    double b = 0.5;
    cout << "Full matrix of solutions: " << endl;
    solution(x, x0, T, T0, h, tau, b);

    return 0;
}
